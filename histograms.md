# Histograms 

Whenver we compute a histogram we are computing multiple averages at the same time.  Consequently, I think that the input syntax syntax that we use for computing histograms should respect the 
distinction between:

1. Computing a time average
2. Computing a spatial average over multiple indistinguishable instances of the same quantity.
3. Computing an average over replicas

that I discussed in the post on [averaging](averaging).  In this post I will explain how I have implemented histograms in a way that respects this distinction.

## Time-average histograms

The code for calculating a histogram is the same as it was in previous versions of PLUMED.  To calculate the distribution of values that the distance between atoms 1 and 2 took during your simulation
you use the following input:

```plumed
x: DISTANCE ATOMS=1,2
hA: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 STRIDE=1
DUMPGRID ARG=hA FILE=histo STRIDE=1 
```

If you look in more detail at this input though you see that the HISTOGRAM action is a shortcut that expands to the following much longer and more complicated looking input:

```plumed
# Compute the distance
x: DISTANCE ATOMS=1,2
# Setup the weight of the kernel that will be added to our histogram
hA_weight: CONSTANT VALUES=1
# Convert the distance into the single kernel that is added to the grid
hA_kde: KDE ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
# Multiply the kernel we created by the weight of this frame
hA_kdep: CUSTOM ARG=hA_kde,hA_weight FUNC=x*y PERIODIC=NO
# Accumulate all the kernels into the final grid
hA_u: ACCUMULATE ARG=hA_kdep STRIDE=1
# Accumulate the sum of all the weights of all frames
hA_nsum: ACCUMULATE ARG=hA_weight STRIDE=1
# And compute the average histogram
hA: CUSTOM ARG=hA_u,hA_nsum FUNC=x/y PERIODIC=NO
DUMPGRID ARG=hA FILE=histo STRIDE=1
```

One can, in fact, go even further and rewrite this input like this:

```plumed
x: DISTANCE ATOMS=1,2
# This bit sets up values that hold the shape of the kernel we are adding.
hA_sigma: CONSTANT VALUES=0.1
hA_cov: CUSTOM ARG=hA_sigma FUNC=x*x PERIODIC=NO
hA_icov: CUSTOM ARG=hA_cov FUNC=1/x PERIODIC=NO
# This bit sets up the height of the Gaussian in a way that ensures that the 
# volume of the kernel is one.
hA_h: CUSTOM ARG=hA_sigma FUNC=1/(x*sqrt(2*pi)) PERIODIC=NO
hA_weight: CONSTANT VALUES=1
hA_kde: KDE ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 HEIGHTS=hA_h METRIC=hA_icov
hA_kdep: CUSTOM ARG=hA_kde,hA_weight FUNC=x*y PERIODIC=NO
hA_u: ACCUMULATE ARG=hA_kdep STRIDE=1
hA_nsum: ACCUMULATE ARG=hA_weight STRIDE=1
hA: CUSTOM ARG=hA_u,hA_nsum FUNC=x/y PERIODIC=NO
DUMPGRID ARG=hA FILE=histo STRIDE=1
```

In other words, the input syntax for these actions is flexibile enough that you can pass the shape of the Gaussians that you want to use when constructing your histogram in the input and the heights of the Gaussians to use.  
You can even pass a matrix in the input for the METRIC keyword if you want to use non-diagonal covariance matrices when adding kernels.  In short, users have an exquisite control over the way their histograms are constructed
if they move away from using the short cuts.  Importantly, they can control how the histogram is normalised or not normalised.  In collaborations I have been a part of I have found that folks want to exploit this flexibilty.  
I also hope that the syntax I will describe in these pages offers this flexibilty while still being quite transparant.  

## Time and spatial average histograms

In the inputs above the DISTANCE command outputs a single scalar.  Creating and storing a grid to hold the single Gaussian centered on this distance value using a KDE command and only then adding it to the histogram that is being 
accumulated by the ACCUMULATE command thus perhaps seems perverse.  When one has multiple indistinguishable representations of the same quantity in each trajectory frame, however, the reason for separating these two actions becomes more clear.  Consider the 
following input that is calculating the distribution of coordination numbers:

```plumed
c1: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=1}
hA_weight: CONSTANT VALUES=1
hA_kde: KDE ARG=c1 GRID_MIN=0.0 GRID_MAX=12.0 GRID_BIN=120 BANDWIDTH=0.1
hA_kdep: CUSTOM ARG=hA_kde,hA_weight FUNC=x*y PERIODIC=NO
hA_u: ACCUMULATE ARG=hA_kdep STRIDE=1
hA_nsum: ACCUMULATE ARG=hA_weight STRIDE=1
hA: CUSTOM ARG=hA_u,hA_nsum FUNC=x/y PERIODIC=NO
DUMPGRID ARG=hA FILE=histo STRIDE=1
```

In the input above the final histogram is caclulated by averaging over time and over the multiple indistinguishable realisations of the coordination number in each trajectory frame.  The KDE action computes a histogram from the instantaneous
values of all the coordination numbers, while the ACCUMULATE then adds all these instantaneous distributions together and hence does the averaging over time.  Notice, furthermore, that we can multiply each instaneous KDE estimate by the 
weight of the corresponding trajectory frame if we are reweighting a biased simulation using a CUSTOM action before we pass the instaneous estimate of the distribution to the ACCUMULATE action.  The input for this type of calculation is pretty 
transparent and flexible as the user does everything with CUSTOM actions (i.e. mathematical equations).  They are not forced to learn an exotic syntax for doing reweighting.

The real reason that I have introduced this separation between the KDE and HISTOGRAM actions is because of how KDE has been used for biasing in papers such as [this one](https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b01027).  In that paper the instantaenous distribution of an order parameter is computed using 
kernel density estimation.  The Kulback Leibler divergence between the instantaneous distribution and some reference distribution is then used as a CV in a biased simulation.  To use a CV like this in PLUMED you would use an input like the one shown below:

```plumed
d1: DISTANCE ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5 ATOMS5=2,3 ATOMS6=2,4 ATOMS7=2,5 ATOMS8=3,4 ATOMS9=3,5 ATOMS10=4,5
d1lt: LESS_THAN ARG=d1 SWITCH={RATIONAL D_0=2.0 R_0=0.5 D_MAX=5.0}

d1c: DISTANCE ATOMS1=2,1 ATOMS2=3,1 ATOMS3=4,1 ATOMS4=5,1 ATOMS5=3,2 ATOMS6=4,2 ATOMS7=5,2 ATOMS8=4,3 ATOMS9=5,3 ATOMS10=5,4 COMPONENTS
d2: COMBINE ARG=d1c.x,d1c.y,d1c.z POWERS=2,2,2 PERIODIC=NO
aa: CUSTOM ARG=d1c.z,d2 FUNC=acos(x/sqrt(y)) PERIODIC=NO

dd0: FIXEDATOM AT=0,0,0
ddx: FIXEDATOM AT=1,0,0
ddz: FIXEDATOM AT=0,0,1

tt: TORSION ...
   VECTORA1=2,1 VECTORB1=ddx,dd0 AXIS1=ddz,dd0
   VECTORA2=3,1 VECTORB2=ddx,dd0 AXIS2=ddz,dd0
   VECTORA3=4,1 VECTORB3=ddx,dd0 AXIS3=ddz,dd0
   VECTORA4=5,1 VECTORB4=ddx,dd0 AXIS4=ddz,dd0
   VECTORA5=3,2 VECTORB5=ddx,dd0 AXIS5=ddz,dd0
   VECTORA6=4,2 VECTORB6=ddx,dd0 AXIS6=ddz,dd0
   VECTORA7=5,2 VECTORB7=ddx,dd0 AXIS7=ddz,dd0
   VECTORA8=4,3 VECTORB8=ddx,dd0 AXIS8=ddz,dd0
   VECTORA9=5,3 VECTORB9=ddx,dd0 AXIS9=ddz,dd0
   VECTORA10=5,4 VECTORB10=ddx,dd0 AXIS10=ddz,dd0
...

hu: KDE VOLUMES=d1lt ARG=aa,tt GRID_BIN=20,20 GRID_MIN=0,-pi GRID_MAX=pi,pi BANDWIDTH=0.2,0.2
de: SUM ARG=d1lt PERIODIC=NO
h: CUSTOM ARG=hu,de FUNC=x/y PERIODIC=NO
h_ref: REFERENCE_GRID FUNC=1 GRID_BIN=20,20 GRID_MIN=0,-pi GRID_MAX=pi,pi PERIODIC=NO,NO
klg: CUSTOM ARG=h,h_ref FUNC=y*log(y/0.5*(x+y)) PERIODIC=NO
kl: INTEGRATE_GRID ARG=klg PERIODIC=NO

RESTRAINT ARG=kl AT=1.0 KAPPA=10
```

This input looks at the instantaneous distribution of $\theta$ and $\phi$ angles for the bonds between atoms and their neighbours in the first coordination sphere.  $\phi$ is defined as the angle between the vector connecting the two atoms and the positive $z$ direction.
$\theta$ is then the torsional angle between the bond and the positive $x$ direction.  In other words, $\theta$ and $\phi$ are the second two spherical polar coordinates.  You can see that we construct the instantaneous distribution of bonds in this spheircal polar direction
and then compute the Kulbeck-Leibler divergence between this instantaneous distribution and a refernece configuration that is read from a file called "reference.grid."  Notice, finally, that we can add a restraint on this very complicated CV and that PLUMED will work out the 
derivatives on the atoms using the chain rule for us.  

## The radial distribution function 

The grid and histogram functionality that I have described in the previous sections can be used to calculate the radial distribution function (RDF).  I have thus implemented a shortcut action to calculate the RDF that works as follows.  If you click on the RDF command you will be given the option
to expand the shortcut and see how this command is built from simpler PLUMED actions.

```plumed
rdf: RDF GROUP=1-1000 GRID_BIN=100 MAXR=1.0 BANDWIDTH=0.01
DUMPGRID ARG=rdf STRIDE=10
```

Notice that what is output here is the instantaneous RDF.  If you want to compute a time average of the RDF you need to use suitable accumulate actions as described in previous sections.

The RDF was implemented in order to add the [entropy CV](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.015701) that was developed by Piaggi and Parrinello.  This CV computes a Kulbeck-Leibler divergence and thus operates like the CVs described at the end of the previous section.  There is a shortcut to this action to which can be called using the following input:

```plumed
pp: PAIRENTROPY GROUP=1-108 MAXR=2.0 GRID_BIN=20 CUTOFF=1.5 BANDWIDTH=0.13
```

If you implement this (and the related PAIRENTROPIES) shortcuts in a single action you may get better performance.  I believe these implementations that use shortcuts provide a useful reference values for these CVs and our easy to maintain.  They are a useful start point for anyone building code to use such a CV in a production calculation.

## Histograms and replicas

If you are running a simulation with four replicas you can construct a histogram using all the data from the replicas by using an input like the one shown below:

```plumed
#SETTINGS NREPLICAS=4
d1: DISTANCE ATOMS=1,2
d1c: GATHER_REPLICAS ARG=d1
d1v: CONCATENATE ARG=d1c.*
wh: CONSTANT VALUES=0.25,0.25,0.25,0.25
hA_weight: CONSTANT VALUES=1
hA_kde: KDE ARG=d1v GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 VOLUMES=wh BANDWIDTH=0.1
hA_kdep: CUSTOM ARG=hA_kde,hA_weight FUNC=x*y PERIODIC=NO
hA_u: ACCUMULATE ARG=hA_kdep STRIDE=1
hA_nsum: ACCUMULATE ARG=hA_weight STRIDE=1
hA: CUSTOM ARG=hA_u,hA_nsum FUNC=x/y PERIODIC=NO
DUMPGRID ARG=hA FILE=histo STRIDE=1
```

Notice how the GATHER_REPLICAS and CONCATENATE actions are used here to transfer the data on each replica to a four dimensional vector that can be used in the input for KDE.  Notice, furthermore, that by using the VOLUMES/HEIGHTS keyword of KDE we can ascribe weights to the data from each replica.
These weights can be used to do any reweighting to take account of different conditions that might be used with different replicas.

## Spatial averages

Over the last few years researchers from Michele's group and I wrote functionality in PLUMED for [this paper]() on the capiliary fluctuation method and [this one]() on determining the shape of a nucleating solid.  These papers made use of an action called MULTICOLVARDENS that built a continuous representation for the average
value of some atom-based order parameter at each point in the simulation cell.  With the new version of PLUMED there is no need to implement a second MULTICOLVARDENS action as you can calculate these fields by using the KDE command directly as shown below:

```plumed
f: FIXEDATOM AT=0,0,0
c1: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=0.1}
d1: DISTANCES COMPONENTS ORIGIN=f ATOMS=1-100
p_numer: KDE VOLUMES=c1 ARG=d1.x,d1.y,d1.z GRID_BIN=100,100,100 BANDWIDTH=0.1,0.1,0.1
p_denom: KDE ARG=d1.x,d1.y,d1.z GRID_BIN=100,100,100 BANDWIDTH=0.1,0.1,0.1
p: CUSTOM ARG=p_numer,p_denom FUNC=x/y PERIODIC=NO
DUMPGRID ARG=p FILE=dens STRIDE=1
```

I hope what is beign calculated here is clear from the syntax above.  Notice, furthermore, that if one wants to compute the average over multiple frames you can accumulate the grids with labels p_numer and p_denom.  On steps where the average of the function is desired you can then compute the ratio of the two grids with CUSTOM.

Notice lastly that this input can also be used to calculate the conditional average of a CV as a function of another CV.  In other words, you can use inputs like the one above to calculate the average value of $s$ when some other CV has a fixed value as is discussed in [this paper]().

## Weights, volumes and normalisation

I hope the examples in the previous sections have demonstrated the flexibility of this new implementation of the histogram in PLUMED.  When designing the syntax for this implementation I have tried to ensure that users can control all aspects of making the histogram. In other words, 
I am trying not to force them to make decisions about the way in which the histogram is constructed.  Furthermore, when I have made decisions I have tried to make the log provide information on the decisions that have been made and there is a way for the user to control those decisions.  My reason for 
desining the code in this way is that I have learned from experience that users will often want histograms that are normalised in different ways.  I thus wanted to avoid writing code that forced the histogram to have a particular normalisation.  Really I want the user to decide how (or if) to normalise
their histogram.

One issue that I think has caused some confusion, and that is worth explaining, is whether the Gaussians that are added when you do KDE have a height of 1 or a volume of one.  The default in PLUMED if you use an input like this one:

```plumed
x: DISTANCE ATOMS=1,2
hA_kde: KDE ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
```

Is that the volume of the added Gaussian is one.  If you want the height of the Gaussian to be one you would write:

```plumed
one: CONSTANT VALUE=1
x: DISTANCE ATOMS=1,2
hA_kde: KDE ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 HEIGHTS=one
```

You can also set the volumes of the added kernels to a value other than one by doing:

```plumed
f: CONSTANT VALUE=4
x: DISTANCE ATOMS=1,2
hA_kde: KDE ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 VOLUMES=f
```

This input is, however, just a shortcut that is equivalent to the following longer input:

```plumed
f: CONSTANT VALUE=4
fh: CUSTOM ARG=f FUNC=x/(sqrt(2*pi)*0.1) PERIODIC=NO
x: DISTANCE ATOMS=1,2
hA_kde: KDE ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 HEIGHTS=fh
```

In conclusion, and if in doubt, you can always control the height of the maximum in the added Gaussian by using the HEIGHTS keyword. 
