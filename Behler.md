# Implementing Behler-Parinello symmetry functions

There are other codes that implement the [Behler-Parinello symmetry functions](https://pubs.aip.org/aip/jcp/article/134/7/074106/954787/Atom-centered-symmetry-functions-for-constructing) and
to my knowledge no one has used these functions as a CV for a metadynamics simulation.  It is, therefore, not unreasonable to ask if an implementation of these CVs in PLUMED is really necessary.
The implementation I will describe below was done relatively quickly and I hope it demonstrates how new functionality can be quickly prototyped in the modified version of PLUMED that I have described
in these pages.

## Implementing $G^1$ symmetry functions

The $G^1$ symmetry function that Behler introduces in the paper that I have linked above is:

$$
G^1_i = \sum_{j \ne i} f_c(R_{ij})
$$

In this expression, $f_c$ is a switching function and $R_{ij}$ is the distance between atom $i$ and $j$.  We have seen this function already, however, it is just the coordination number.  To calculate it 
using PLUMED we can use the following input.

```plumed
cmat: CONTACT_MATRIX GROUP=1-100 SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)}
ones: ONES SIZE=100
g1: MATRIX_VECTOR_PRODUCT ARG=cmat,ones
# This will print out the 100 coordination number values calculated by the input above.
PRINT ARG=g1 FILE=colvar
```

We have thus already implemented the first of Behler's symmetry functions.

## Implmenting radial symmetry functions

The $G^2$ and $G^3$ symmetry functions that Behler introduces in his paper are:

$$
G^2_i = \sum_{j \ne i} e^{-\nu(R_{ij} -R_s)^2} f_c(R_{ij}) \qquad \textrm{and} \qquad G^3_i = \sum_{j \ne i} \cos(\kappa R_{ij}) f_c(R_{ij})
$$

$\nu$, $R_s$ and $\kappa$ here are parameters, while $f_c$ is a switching function that acts on the distance $R_{ij}$ between atom $i$ and $j$.
In calculating these two functions we need to do a sum of a function of all the $R_{ij}$ values.  We thus have everything we need to calculate 
these variables within PLUMED and need to implement nothing new.  We can simply use the following code to compute $G^2$ and $G^3$ with $\nu=1$, $r_s=3$ and $\kappa=1$:

```plumed
# Calculate the contact matrix.  Element i,j of this matrix tells you if atoms i and j are within a certain cutoff of each other
cmat: CONTACT_MATRIX GROUP=1-100 COMPONENTS SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)}
# Now calculate a matrix with all the R_ij values
cmatr: CUSTOM ARG=cmat.x,cmat.y,cmat.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
# Compute the quantity in the summation in the expression for G^2.  The output here is a 100x100 matrix.
g2_f: CUSTOM ARG=cmatr,cmat.w FUNC=y*exp(-(x-3)^2) PERIODIC=NO
# Compute the quantity in the summation in the expression for G^3.  The output here is a 100x100 matrix.
g3_f: CUSTOM ARG=cmatr,cmat.w FUNC=y*cos(x) PERIODIC=NO
# Now multply the matrices above by a vector of all ones to do the summations
ones: ONES SIZE=100
g2: MATRIX_VECTOR_PRODUCT ARG=g2_f,ones
g3: MATRIX_VECTOR_PRODUCT ARG=g3_f,ones
# Print out the 100 values for g2
PRINT ARG=g2 FILE=g2_file
# Print out the 100 values for g3
PRINT ARG=g3 FILE=g3_file
``` 

Notice the flexibility of this implementation for the Behler symmetry function you can basically compute any weighted function of the bond vectors in the first 
coordination sphere using an input like this one.  A similar input (see below) could be used to calculate the FCC cubic parameter:

```plumed
# Calculate the contact matrix.  Element i,j of this matrix tells you if atoms i and j are within a certain cutoff of each other
cmat: CONTACT_MATRIX GROUP=1-100 COMPONENTS SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)}
# Now calculate a matrix with all the R_ij values
cmatr: CUSTOM ARG=cmat.x,cmat.y,cmat.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
# Now calculate the fcc cubic parameter for each bond.  This outputs a matrix
fcc_f: CUSTOM ARG=cmat.x,cmat.y,cmat.z,cmatr FUNC=((x^4)*(y^4)+(x^4)*(z^4)+(y^4)*(z^4))/(r^8)-27*(x^4)*(y^4)*(z^4)/(r^12) VAR=x,y,z,r PERIODIC=NO
# And multiply the above by the weights to get another 100 x 100 matrix
wfcc_f: CUSTOM ARG=fcc_f,cmat.w FUNC=x*y PERIODIC=NO
# Now sum the rows of the matrix above to get a vector with 100 elements.  One symmetry function for each atom. 
ones: ONES SIZE=100
fcc_u: MATRIX_VECTOR_PRODUCT ARG=wfcc_f,ones 
# This is the coordination number
denom: MATRIX_VECTOR_PRODUCT ARG=cmat.w,ones
# We devide fcc_u by the coordination number to get the average value of the function above for the atoms in the first coordination sphere.
fcc: CUSTOM ARG=fcc_u,denom FUNC=x/y PERIODIC=NO
# And finally we print the 100 symmetry function values to a file called colvar.
PRINT ARG=fcc FILE=colvar
```

The code is faster if the calculations that are done with CUSTOM keywords in the input above are implemented directly in C++.  Even so the fact that you can so quickly prototype complicated 
CVs like these is (I hope) useful.

## Implementing angular symmetry functionsa

One of the angular symmetry functions that Behler introduces is:

$$
G^5_i = 2^{1-\zeta} \sum_{j,k\ne i} (1 + \lambda\cos\theta_{ijk})^\zeta e^{-\nu(R_{ij}^2 + R_{ik}^2)} f_c(R_{ij}) f_c(R_{ik})
$$

In this expression $\zeta$, $\nu$ and $\lambda$ are all parameters.  $f_c$ is a switching function which acts upon $R_{ij}$, the distance between atom $i$ and atom $j$, and
$R_{ik}$, the distance between atom $i$ and atom $k$.  $\theta_{ijk}$ is then the angle between the vector that points from atom $i$ to atom $j$ and the vector that points from 
atom $i$ to atom $k$. 

Even this complicated function can be calculated directly in PLUMED.  The input below caluclates the $G^5_i$ value for atom 1 with $\zeta=3$, $\lambda=3 and $\nu=0.1$.

```plumed
# Calculate the distances between atom one and all the other atoms
d1: DISTANCE ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5 ATOMS5=1,7 ATOMS6=1,8
# Now transform the distances above by a switching function
d1lt: LESS_THAN ARG=d1 SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)}
# Calculate the vectors between atom one and all the other atoms
d1c: DISTANCE COMPONENTS ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5 ATOMS5=1,7 ATOMS6=1,8 
# Now calculate the lengths of all the vectors connecting atom 1 to the other atoms
d1r2: COMBINE ARG=d1c.x,d1c.y,d1c.z POWERS=2,2,2 PERIODIC=NO
d1r: CUSTOM ARG=d1r2 FUNC=sqrt(x) PERIODIC=NO
# And find unit vectors that connect atom 1 to all its neighbours
d1ux: CUSTOM ARG=d1c.x,d1r FUNC=x/y PERIODIC=NO
d1uy: CUSTOM ARG=d1c.y,d1r FUNC=x/y PERIODIC=NO
d1uz: CUSTOM ARG=d1c.z,d1r FUNC=x/y PERIODIC=NO 
# And put the three unit vectors into a matrix
stack: VSTACK ARG=d1ux,d1uz,d1uz
# Now create a matrix in which element i,j = f_c(R_ij) f_c(R_ik)
wmat: OUTER_PRODUCT ELEMENTS_ON_DIAGONAL_ARE_ZERO ARG=d1lt,d1lt
# And a matrix containing the cosines of the angles between the two vectors
stackT: TRANSPOSE ARG=stack
cmat: MATRIX_PRODUCT ARG=stack,stackT
# Now compute a matrix containing (1 + \lambda cos(\theta_ijk))^3
pmat: CUSTOM ARG=cmat PERIODIC=NO FUNC=(1+3*cos(x))^3 
# And a matrix containing in which element j,k is R_ij^2 + R_ik^2 
smat: OUTER_PRODUCT ARG=d1r2,d1r2 FUNC=x+y 
# Now take the exponential of this matrix
emat: CUSTOM ARG=smat FUNC=exp(-0.1*x) PERIODIC=NO
# And combine everything to get G^5
g5mat: CUSTOM ARG=cmat,emat,wmat FUNC=x*y*z PERIODIC=NO
# Now sum all the elements of g5mat
g5s: SUM ARG=g5mat PERIODIC=NO
# And multiply by 0.5 and 2^{1-\zeta}.  We need to multiply by 0.5 here to avoid double counting
g5: CUSTOM ARG=g5s FUNC=0.5*0.25*x PERIODIC=NO
# And print the final scalar to a file
PRINT ARG=g5 FILE=colvar
```  

This input computes $G^5$ for a single atom and is rather cumbersome.  The calculation here is not particularly well optimised as we have to calculate the distance between atom 
1 and all the other atoms.  In other words, when an input similar to this one is used PLUMED is not using link cells or neighbour lists to optimse the calculation.  This way of 
implementing $G^5$ is thus not really suited to production applications.  I thus wrote a second implementation of the $G^5$ function.  The input for this alternative implementation is 
as follows:

```plumed
# Calculate the contact matrix and the x,y and z components of the bond vectors
# This action calculates 4 100x100 matrices
cmat: CONTACT_MATRIX GROUP=1-100 SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)} COMPONENTS

# Compute the symmetry function for the 100 atoms from the 4 100x100 matrices output
# by cmat.  The output from this action is a vector with 100 elements
beh3: GSYMFUNC_THREEBODY ...
    WEIGHT=cmat.w ARG=cmat.x,cmat.y,cmat.z
    FUNCTION1={FUNC=0.25*exp(-0.1*(rij+rik))*(1+3*cos(ajik))^3 LABEL=g5}
...

# Print the 100 symmetry function values to a file
PRINT ARG=beh3.g5 FILE=colvar
```

The GSYMFUNC_THREEBODY action sums over all the distinct triples of atoms that are identified in the contact matrix.  This action auses lepton and can thus compute any 
function of the following four quantities:

* `rij` - the distance between atom $i$ and atom $j$
* `rik` - the distance between atom $i$ and atom $k$
* `rjk` - the distance between atom $j$ and atom $k$
* `ajik` - the angle between the vector connecting atom $i$ to atom $j$ and the vector connecting atom $i$ to atom $k$.

Furthermore we can calculate more than one function of these four quantities at a time as illustrated by the input below:

```plumed
# Calculate the contact matrix and the x,y and z components of the bond vectors
# This action calculates 4 100x100 matrices
cmat: CONTACT_MATRIX GROUP=1-100 SWITCH={CUSTOM R_0=4.5 D_MAX=4.5 FUNC=0.5*(cos(pi*x)+1)} COMPONENTS

# Compute the 4 symmetry function below for the 100 atoms from the 4 100x100 matrices output
# by cmat.  The output from this action is a vector with 100 elements
beh3: GSYMFUNC_THREEBODY ...
    WEIGHT=cmat.w ARG=cmat.x,cmat.y,cmat.z
    FUNCTION1={FUNC=0.25*(cos(pi*sqrt(rjk)/4.5)+1)*exp(-0.1*(rij+rik+rjk))*(1+2*cos(ajik))^2 LABEL=g4}
    FUNCTION2={FUNC=0.25*exp(-0.1*(rij+rik))*(1+3.5*cos(ajik))^3 LABEL=g5}
    FUNCTION3={FUNC=0.125*(1+6.6*cos(ajik))^4 LABEL=g6}
    FUNCTION4={FUNC=sin(3.0*(ajik-1)) LABEL=g7}
...

# Print the 4 sets of 100 symmetry function values to a file
PRINT ARG=beh3.g4,beh3.g5,beh3.g6,beh3.g7 FILE=colvar
```

As we saw for the radial symmetry function, we thus again have an implementation of these angular symmetry functions that can be quickly used to prototype new function types.

## Conclusions

I hope the examples above have convinced you that the PLUMED input file now provides the flexibility required to prototype CVs that are really rather complicated.  Using 
CUSTOM actions when implementing these CVs is probably not optimal when it comes to computational performance but I am not convinced that this matters.  I think our objective should always be to have
less code as this makes PLUMED easier to maintain.  

It is also the case that many of the CVs that are implemented in PLUMED are not used by many people.  The understanding of what particular actions do is thus lost with time.  When CVs are implemented in 
the way above the series of action that are performed within them are more transparent.  This transparency allows allows other to understand what has been done in the past and reduces the requirements to 
write large amounts of documentation.  

Lastly, if users need code that runs faster they can always reimplement expensive CVs that are implemented in PLUMED.  Furthermore, in any effort to reimplement CVs they are helped
by the fact that PLUMED offers reference values that they can test their new code against.
