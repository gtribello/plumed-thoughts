# Paths

Paths like PCA coordinates are another one of these general methods that can be setup in a variety of different ways.  In other words, paths can be setup in a space of 
collective variables or they can be performed in RMSD space directly.  I thus wanted the new version of the PATH command to a shortcut that generates a more complicated 
input files.  Experienced users can then use their understanding of this more complex input to develop new versions of PATH CVs that use different metrics for measuring
the distances between the waymarks and the reference configurations.

The following input shows an example of the simple input that you can use for the PATH command:

```plumedfile
p: PATH REFERENCE=all.pdb LAMBDA=69087
```

The actual input that plumed uses to run this calculation is as follows:

```plumed
# Calculate the rmsd distance between all the instantaneous configuration
# and all the waypoints that are in the file all.pdb.  This action returns a 
# vector of squared distances from these waypoints
p_data: RMSD SQUARED REFENCE=all.pdb
# Find the shortest distance between the path and the reference configurations
p_mindist: LOWEST ARG=p_data 
# Now calcluate the weight for each point in the path
p_weights: CUSTOM ARG=p_data,p_mindist FUNC=exp(-(x-y)*69087) PERIODIC=NO
# This is the denominator in the expression for s
p_denom: SUM ARG=p_weights PERIODIC=NO
# This computes the z collective variable
p_z: CUSTOM ARG=p_denom,p_mindist FUNC=y-log(x)/69087 PERIODIC=NO
# These are the positions of the waypoints for the path 
p_ind: CONSTANT VALUES=1,2,3,4
# And now compute s
p_numer_prod: CUSTOM ARG=p_weights,p_ind FUNC=x*y PERIODIC=NO
p_numer: SUM ARG=p_numer_prod PERIODIC=NO
p_s: CUSTOM ARG=p_numer,p_denom FUNC=x/y PERIODIC=NO 
```

If you want to use a CV space rather than using RMSD distances you just need to swap out the RMSD command in the above.  Every subsequent step in the method for calculating the 
path CV is the same.  This fact is recognised by the PATH command that will also work if you use a path that is defined in a space of collective variables and does indeed just 
swap out the RMSD action for an action that calculates EUCLIDEAN_DISTANCES between the instantaneous and reference values of the various CVs.

## Geometric and adaptive paths

The method for calculating the position on the path is the original version of this method that was developed by Branduardi, Gervasio and Parrinello.  There is another method 
for calculating the position along and the distance from a path that was developed by Leines and Ensing.  A similarly flexible implementation of this method is available in PLUMED.
Once again shortcut actions have been used in the implementation as you can see by exploring the input below:

```plumed
t1: TORSION ATOMS=5,7,9,15
t2: TORSION ATOMS=7,9,15,17

# Path based on Eculidean distances
pp: GPATH ARG=t1,t2 REFERENCE=epath.pdb

# Path based on rmsd distances 
p2: GPATH TYPE=OPTIMAL-FAST REFERENCE=all.pdb
```

This type of path is used in the adaptive path method that was developed by Leines and Ensing.  This too is available in PLUMED as a shortcut action.  The appropriate command is
ADAPTIVE_PATH 
