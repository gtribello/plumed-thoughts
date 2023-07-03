# Implementing SPRINT collective variables

SPRINT collective variables were developed by Fabio Pietrucci and Wanda Andreoni and introduced in this paper in 2011.
I implemented them in PLUMED shortly after.  I don't think many people have found a use for these variables so there is no 
pressing need to develop a very fast implementation of these CVs.  I think a better approach is thus to implement these variables
in a way that helps potential future users understand what was done in Peitrucci and Andreoni's original paper.  For this reason,
the action SPRINT is implemented as a shortcut.  This means calculating and printing the SPRINT CVs is straightforward.  You 
just use the following input file:

```plumed
s1: SPRINT GROUP1=1-7 SWITCH11={RATIONAL R_0=2.6 NN=6 MM=12}   
PRINT ARG=s1.* FILE=colvar
```

When PLUMED reads this it creates the more complicated input file shown below for calculating the SPRINT CVs.

```plumed
# Calculate a contact matrix that tells us about the bond between these 7 atoms
s1_mat: CONTACT_MATRIX GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12} 
# Calculate the principle eigenvalue of the contact matrix and its corresponding eigenvector
s1_diag: DIAGONALIZE ARG=s1_mat VEC=1
# Scale the eigenvector by multiplying it by the square root of the number of atoms and the corresponding eigenvalue 
s1_sp: CUSTOM ARG=d1_diag.vals-1,s1_diag.vecs-1 FUNC=sqrt(7)*x*y PERIODIC=NO
# Sort the compnents of scaled eigenvector 
s1: SORT ARG=s1_sp
```

A user looking in the log for PLUMED can thus see the individual steps involved in calculating this CV. In other words,
they can use PLUMED to get a sense of how these legacy CVs are calculated.  In fact, if you look at the first example input
above you should see the option to expand the shortcut in the tooltip that appears when you hover above SPRINT.  If you click
on this option something similar to the second input above appears.  Users thus do not even need to look in the log to see how 
these complicated CVs are calculated.  They can instead look through examples where such CVs have been used in the manual, tutorials 
and nest and use this option to expand shortcuts.

To get a sense of just how much PLUMED input can be hidden in these short inputs consider the following input for calculating 
SPRINT coordinates:

```plumed
ss: SPRINT ...
   GROUP1=1-7 GROUP2=8-14 
   SWITCH11={RATIONAL R_0=2.6 NN=6 MM=12}
   SWITCH12={RATIONAL R_0=2.2 NN=6 MM=12}
   SWITCH22={RATIONAL R_0=2.2 NN=6 MM=12}
...
PRINT ARG=ss1.*,ss2.* FILE=colvar
```

In this atom there are two distinct types of atom (A and B).  We thus have three different switching functions for describing 
whether the AA, AB and BB bonds.  If you look at the expanded version of the input above by clicking on the expand shortcut option
for sprint you can see how the contact matrices for AA, AB and BB are constructed separately and then concatentated into a single 
object.  You can also see how the process of sorting the components of the eigenvector has to be modified to take account of the fact 
that all the atoms are not indistinguishable.  In conclusion, my hope is that readers can make sense of the way CVs are constructed from
shortcut in PLUMED by reading the original paper and the PLUMED input files together.



