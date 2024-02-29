# RMSD and PCAVARS

One of the things that has bothered me for some time RMSD. RMSD is used in many different methods.  Furthermore, in many of these methods the RMSD can be replaced with a simple 
Eudlidean distance that is calculated by taking the difference between two sets of argument values. For these methods I always wanted to write one impelementation that would work with 
both RMSD and Euclidean distances. However, I feel like the impelementations of these things that I wrote in the past were always an overcomplicated mess.  I now think that I have written 
something simpler that still provides the same flexibility.  I would like to explain how this works by discussing how PCAVARS are implemented.

Before getting into the code lets first briefly recap the theory. When you calculate a PCA variable you are essentially projecting a vector connecting the instantaneous configuration to some 
reference configuration on to a reference vector.  In other words, a PCA variable is computed by taking the a dot product like this

$$
s = \sum_{i=1}^n v_i (x_i - y_i)
$$

where $v_i$ represents the components of the reference vector, $x_i$ is the instantenous position and $y_i$ is the reference position.  When this is done with RMSD you do an alignment (i.e. a removal
of translation and rotation of atoms) before calculating the $(x_i - y_i)$ values by calculating how far each position has moved relative to the reference structure.

You can implement this type of CV in plumed by using the following plumed input.

```plumed
v: CONSTANT VALUES=0.123253,0.077840,0.018043,0.027632,0.290613,-0.095447,-0.279804,-0.530408,-0.026324,-0.124997,0.027980,0.068251,-0.054555,-0.094179,-0.041060,-0.033738,-0.392797,-0.306944,-0.303786,-0.401124,-0.000431,0.010767,-0.083125,-0.055704 NCOLS=12 NROWS=2
rmsd: RMSD REFERENCE=reference.pdb TYPE=OPTIMAL DISPLACEMENT SQUARED
pca: MATRIX_VECTOR_PRODUCT ARG=v,rmsd.disp
PRINT ARG=pca FILE=colvar
``` 

The reference.pdb file for the above input looks as follows:

````
ATOM      2  CH3 ACE     1      12.932 -14.718  -6.016  1.00  1.00
ATOM      5  C   ACE     1      21.312  -9.928  -5.946  1.00  1.00
ATOM      9  CA  ALA     2      19.462 -11.088  -8.986  1.00  1.00
ATOM     13  HB2 ALA     2      21.112 -10.688 -12.476  1.00  1.00
END
````

The CV is thus based on the displacement in the positions of these four atoms from the reference structure.  The command RMSD calculates the displacements of these four atoms, which are the $(x_i - y_i)$ values from the equation above.  I then calculate the projection on the 
vectors of interest by doing matrix multiplication.  In the case above I calculate the projection on two vectors as my matrix has two rows.

A similar input can be used to calculate these PCA projections in a CV space as shown below.  For example in the input file below I am doing a projection of a vector between the instantaneous configuration and a reference in a space defined by three distances on a reference vector.

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,2
d3: DISTANCE ATOMS=1,2

d1_ref: CONSTANT VALUE=0.1221 
d2_ref: CONSTANT VALUE=0.0979
d3_ref: CONSTANT VALUE=0.1079

v: CONSTANT VALUES=0.078811,-0.945732,-0.315244 NCOLS=3 NROWS=1
disp: DISPLACEMENT ARG1=d1,d2,d3 ARG2=d1_ref,d2_ref,d3_ref
pca: MATRIX_VECTOR_PRODUCT ARG=v,disp
PRINT ARG=pca
```

The only difference between this input and the previous that used the RMSD displacement is the line that calculates the displacement.  In the previous input we used RMSD whereas here we use DISPLACEMENT instead. For the user the fact that there is only a small difference in the two 
methods is hopefully clear. In addition, when things are implemented in this way it is hopefully easy for users to work with exotic methods of calculating the displacement vector (e.g. combinations of RMSD displacements and differences in arguments).  They can implement these new metrics
without writing new code by working in the input file directly.

If users want a simpler syntax, however, there is a shortcut provided that allows them to get use the older syntax for PCAVARS.  For example, the input below will compute the same quantity as the second input above:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6

pca: PCAVARS ARG=d1,d2,d3 REFERENCE=epath.pdb
PRINT ARG=pca_eig FILE=colvar
```

You can also use PCAVARS to compute the first input above as follows:

```plumed
pca: PCAVARS REFERENCE=reference.pdb TYPE=OPTIMAL
PRINT ARG=pca.* FILE=colvar 
```

## Inside the RMSD command

The way the RMSD command works had to be changed slightly to make inputs like the one above work.  To be clear, however, if you use the following input:

```plumed
rmsd0: RMSD TYPE=OPTIMAL REFERENCE=reference.pdb 
```

Then what PLUMED does during the calclulate step is unchanged.  The only difference is that the command RMSD is a shortcut action.  If there is a single frame in the reference
PDB file and if the DISPLACEMENT option is not used then an RMSD_SCALAR action object is created.  This RMSD_SCALAR object is the old and familiar implementation of the RMSD colvar.

If, however, you use call RMSD with the DISPLACEMENT option as shown below:

```plumed
rmsd: RMSD REFERENCE=reference.pdb TYPE=OPTIMAL DISPLACEMENT SQUARED
```

Then the following input is created by the shortcut:

```plumed
# Converts the input pdb file to a constant
rmsd_ref: PDB2CONSTANT REFERENCE=reference.pdb
# Fix any broken bonds
WHOLEMOLECULES ENTITY0=2,5,9,13
# Get the instaneous positions of the atoms of interest
rmsd_cpos: POSITION NOPBC ATOMS=2,5,9,13
# Concatenate all the instaneous positions into a single vector 
rmsd_pos: CONCATENATE ARG=rmsd_cpos.x,rmsd_cpos.y,rmsd_cpos.z
# Compute the RMSD
rmsd: RMSD_VECTOR ARG=rmsd_pos,rmsd_ref DISPLACEMENT ALIGN=1,1,1,1 DISPLACE=1,1,1,1
```

This input is necessary because the RMSD_VECTOR command is an ActionWithArguments rather than a ActionAtomistic.  the RMSD_VECTOR command takes two vectors with equal lengths as input.
These vectors contain the instantaneous and reference configurations.  

Having RMSD as an ActionWithArguments may initially seem odd.  This change allows you to use this action in a number of different ways, however.  For example, you can read in two reference 
configurations usign the PDB2CONSTANT command and then calculate the RMSD vector connecting these two reference configurations using the RMSD command.  Better still you can pass arguments containing structures 
to plumed and then calculate the RMSD RMSD distance between the two input structures using RMSD directly.  I have used this functionality to reimplement the code in pathtools and adaptive path.  It can also be used 
when implementing dimensionality reduction algorithms such as MDS.  
