# Calculating averages

The final result that is extracted from a molecular dynamics simulation is an average.  In other words, molecular dynamics (like Monte Carlo) is a method for generating multiple samples
of a random variable so that one can take an average.  This average then provides an estimate for a physical observable, which should really be calculated by computing the ensemble average.

When we run molecular dynamics simulations there are three methods that we can use to generate the multiple samples that are used to calculate these averages:

1. We can calculate observables for each of the microstates in the trajectory and take a time average.
2. We can have multiple indistinguishable representations of the quantity of interest in each microstate we generate.  We can thus perform a spatial average over these indistinguishable entities.
3. We can run multiple trajectories in parallel and then average over these distinct replicas of the whole physical system.

The various methods that are implemented in PLUMED combine these three types of averaging in a myriad of different ways. To make the code more transparent I have thus tried to separate these various
ways of collecting data for averaging.  My hope that these separate actions can be used when implementing new methods for calculating averages and that the resulting code will provide a more transparent
explanation of the similarities and differences between all these methods.  Obviously, shortcut actions can be used to enusre that the complex inputs that are actually used to run the various methods are generated
correctly.  My hope, however, is that if folks are willing to use shortcuts when implementing new methods we can generate rich documentation that illustrates how the averages that are used in the simulations
are generated from the trajectory.  This documentation will, I hope, help users to better understand the theory behind the methods that they are employing and encourage cross-fertilisation of ideas between 
developers of different methods.

With all that in mind, allow me to explain here how the three forms of averaging described above are implemented in PLUMED.

## Unbiased time averages

Lets start with simple PLUMED input. The input below calculates and prints the time average of the distance between atom 1 and 2 to a file called average.

```plumed
d1: DISTANCE ATOMS=1,2
a1: AVERAGE ARG=d1 STRIDE=1
PRINT ARG=a1 FILE=average STRIDE=1
```

As the stride is set equal to one on the AVERAGE line samples are taken from every microstate that is visited in the trajectory.  If we had used the following input:

```plumed
d1: DISTANCE ATOMS=1,2
a1: AVERAGE ARG=d1 STRIDE=10
PRINT ARG=a1 FILE=average STRIDE=1
```

We would only accumulate averages on every 10th step.  We are still printing a1 on every step, however.  We would thus find that the values of a1 that are printed out on steps
10,11,12,13,14,15,16,17,18 and 19 were all the same.  No additional samples were added on these steps so the average is unchanged.

The input for the AVERAGE action is the same as it was before we merged all the new code I have been working on.  This action is now a shortcut, however.  The full input for the first of the two 
inputs above is:

```plumed
d1: DISTANCE ATOMS=1,2
a1_weight: CONSTANT VALUES=1
a1_prod: CUSTOM ARG=d1,a1_weight FUNC=x*y PERIODIC=NO
a1_denom: ACCUMULATE ARG=a1_weight STRIDE=1
a1_numer: ACCUMULATE ARG=a1_prod STRIDE=1
a1: CUSTOM ARG=a1_numer,a1_denom FUNC=x/y PERIODIC=NO
PRINT ARG=a1 FILE=average STRIDE=1
```

Notice that the non-standard action that is used here is the one called ACCUMULATE.  This action simply adds together the instantaneous values that are passed to it.  You can use this action to 
do a sum over a time series.  The rest of the averaging operation is done with CUSTOM actions.  You can also hopefully see how the input above can easily be modified so that you can calculate 
an average from a biased simulation.  In that case you simply need to modify the input for the a1_weight command in a way that accounts for the reweighting that must be done on the biased simulation.

Notice that it is also easier to explain how block averages are computed with this input.  When describing an input like this one:

```plumed
d1: DISTANCE ATOMS=1,2
a1_weight: CONSTANT VALUES=1
a1_prod: CUSTOM ARG=d1,a1_weight FUNC=x*y PERIODIC=NO
a1_denom: ACCUMULATE ARG=a1_weight STRIDE=1 CLEAR=100
a1_numer: ACCUMULATE ARG=a1_prod STRIDE=1 CLEAR=100
a1: CUSTOM ARG=a1_numer,a1_denom FUNC=x/y PERIODIC=NO
PRINT ARG=a1 FILE=average STRIDE=1
```

You simply state that the CLEAR=100 keywords ensure that a1_denom and a1_numer are set equal to zero on every hundreth step. 

Lastly, note that, because we are using shortcuts, we can show users what is done for special cases.  To see what I mean consider the following input:

```plumed
t1: TORSION ATOMS=1,2,3,4
a1: AVERAGE ARG=t1 STRIDE=1
PRINT ARG=a1 FILE=average STRIDE=1
```

When this is expanded you have the following:

```plumed
t1: TORSION ATOMS=1,2,3,4
a1_weight: CONSTANT VALUES=1
a1_denom ACCUMULATE ARG=a1_weight 
a1_sin: CUSTOM ARG=t1,a1_weight FUNC=y*sin((x+pi)/1) PERIODIC=NO
a1_sin: CUSTOM ARG=t1,a1_weight FUNC=y*cos((x+pi)/1) PERIODIC=NO
t1_sinsum: ACCUMULATE ARG=a1_sin STRIDE=1
t1_cossum: ACCUMULATE ARG=a1_cos STRIDE=1
a1: CUSTOM ARG=a1_sinsum,a1_cossum,a1_denom FUNC=-pi+1*atan2(x/z,y/z) PERIODIC=-pi,pi
PRINT ARG=a1 FILE=average STRIDE=1
```

Readers of the manual can thus see how the average is calculated in a way that takes the periodicity of the variable into account without having to look inside the PLUMED source code.

## Spatial averages

Suppose that you have an input file that calculates the coordination numbers of 100 atoms like the one shown below:

```plumed
c1: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1}
a1: ACCUMULATE ARG=c1 STRIDE=1
```

Should `a1`, the value that is output by the ACCUMULATE action here be a scalar or a vector?  In other words, should PLUMED assume that the user wants to add all 100 numbers calculated by the 
COORDINATIONNUMBER command all together when ACCUMULATE is called or should it instead assume that the user wants to accumulate individual averages for each of the 100 coordination numbers?

To me the answer to this question is obvious.  The ACCUMULATE action should output a vector.  In other words, the first element of the vector `a1` should be the sum of all the values the coordination
number of atom 1 took, the second element should be the sum of the vlaues the coordination number of atom 2 took and so on.  If the user had wanted us to convert the vector computed by the 
COORINDATIONNUMBER action to a scalar they would have told us they were doing this by using an input like this:

```plumed
c1: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1}
c1_mean: MEAN ARG=c1 PERIODIC=NO
a1: ACCUMULATE ARG=c1_mean STRIDE=1
``` 

I think this keeps things transparent for users. I wouldn't even add shortcuts around for this sort of functionality. I think that the flexibility provided here for folks doing a mixture of time and spatial 
averaging makes things pretty straightforward and transparent.

## Averages over replicas

The example input file below shows how you can gather all the data from all replicas in one place using the GATHER_REPLICAS command:

```plumed
d1: DISTANCE ATOMS=1,2
r1: RESTRAINT ARG=d1 AT=@replicas:1.0,1.2,1.4,1.6 KAPPA=10
r1g: GATHER_REPLICAS ARG=r1.bias
PRINT ARG=r1g.rep-1,,r1g.rep-2,r1g.rep-3,r1g.rep-4 FILE=colvar
```

The GATHER_REPLICAS command creates one PLUMED value object for each replica. The $N$ values output by GATHER_REPLICA contain the values of the quantities on each of the $N$ replicas.  By using GATHER_REPLICAS you 
can meke the values that you have on other replicas accessible on all replicas.  Consequently, if you can run a four replica simulation to compute the average distance between atom 1 and atom 2 you can get the average by doing:


```plumed
d1: DISTANCE ATOMS=1,2
d1g: GATHER_REPLICAS ARG=d1
d1s: COMBINE ARG=d1.* PERIODIC=NO
a1_weight: CONSTANT VALUES=4
a1_denom: ACCUMULATE ARG=d1s STRIDE=1
a1_numer: ACCUMULATE ARG=a1_prod STRIDE=1
a1: CUSTOM ARG=a1_numer,a1_denom FUNC=x/y PERIODIC=NO
PRINT ARG=a1 FILE=average STRIDE=1
```

Flexibility is again offered to users here by separating the act of gathering data from replicas and adding it all together in the input file.  Also notice that if the quantity input to GATHER_REPLICAS is a vector/matrix or grid
then the output actions will be vectors matrices or grids.  In short, the syntax here agsin offers users a great deal of flexibility by separating out these three ways of collecting samples for averaging.  This flexibility hopefully also
helps users to understand how the averages that they are computing are being constructed.
