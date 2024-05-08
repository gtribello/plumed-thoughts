# Clustering

When studying nucleation it is often useful to use a clustering atoms to determine how many atoms are in the largest crystalline nucleus.
The implementation of this approach in PLUMED is detailed in [this paper](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b01073).  A typical 
input in that paper for calcluating the number of atoms in the largest cluster is shown below:

```plumed
# Ccalculate the coordination numbers of the atoms
lq: COORDINATIONNUMBER SPECIES=1-100 SWITCH={CUBIC D_0=0.45 D_MAX=0.55}
# Calculate the contact matrix for the atoms for which we calculated the coordinaion numbers
cm: CONTACT_MATRIX GROUP=lq SWITCH={CUBIC D_0=0.45 D_MAX=0.55}
# Do a clustering using the contact matrix above
dfs: DFSCLUSTERING MATRIX=cm
# Sum the coordination numbers for the atoms in the largest cluster
clust1: CLUSTER_PROPERTIES CLUSTERS=dfs ARG=lq CLUSTER=1 SUM
```

This input is fine but it is also somewhat unweildy and a little confusing.  The problem is that you have to calculate the coordination numbers 
of all the atoms in order to do the clustering and (unless you have a deep understanding of the way the code is implemented) it is not clear why.
With the new sytax you can achieve the same result as follows:

```plumed
# Calculate the contact matrix.  This action computes a 100x100 matrix
cm: CONTACT_MATRIX GROUP=1-100 SWITCH={CUBIC D_0=0.45 D_MAX=0.55}
# Do a clustering using the contact matrix that was computed above as input
# This action returns a 100 dimensional vector. If element i of this matrix
# is equal to 5 this means that atom i in the input to the contact matrix above
# is part of the 5th largest cluster.  
dfs: DFSCLUSTERING ARG=cm
# This next action returns a vector with 100 elements. If element i is equal to 1 then atom 
# i is part of the largest cluster.  If it is equal to zero then it is part of some 
# other cluster.
c1: CLUSTER_WEIGHTS CLUSTERS=dfs CLUSTER=1
# Now calculate the coordination numbers using the usual matrix multiplication trick
ones: ONES SIZE=100
coords: MATRIX_VECTOR_PRODUCT ARG=cm,ones
# Multiply the coordination numbers by c1.  We now have a vector where element i is equal to the 
# coordiation number of atom i if atom i is part of the largest cluster and zero otherwise.
fcoords: CUSTOM ARG=coords,c1 FUNC=x*y PERIODIC=NO
# And lastly sum the coordination numbers of the atoms in the largest cluster
coordsum: SUM ARG=fcoords PERIODIC=NO
```

This new syntax is much more clear as the clustering operation is performed on the CONTACT_MATRIX directly.  The vector returned by the DFSCLUSTERING object 
then tells you which cluster each atom belongs to.  You can thus use simple logical operations on this vector to determine the properties for all your clusters. 
Furthermore, you don't even need to use the coordination numbers.  If you simply want to calculate the number of atoms in the largest cluster you can use the following
input:

```plumed
# Calculate the contact matrix.  This action computes a 100x100 matrix
cm: CONTACT_MATRIX GROUP=1-100 SWITCH={CUBIC D_0=0.45 D_MAX=0.55}
# Do the clustering
dfs: DFSCLUSTERING ARG=cm
# Get a 100 element vector that has ones for those atoms that are part of the largest cluster
c1: CLUSTER_WEIGHTS CLUSTERS=dfs CLUSTER=1   
# Sum the vector above to get the number of atoms in the largest cluster
suml: SUM ARG=c1 PERIODIC=NO
```
