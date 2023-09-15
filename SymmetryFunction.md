# Calculating symmetry functions

In the other blog articles on these pages I have written a great deal about how you can calculate coordination numbers 
using inputs such as the one below:

```plumed
# Calculate the contact matrix.  This is a square 7x7 matrix
c1: CONTACT_MATRIX GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
# Calculate the coordination numbers for the 7 atoms by a vector that contains all ones
ones: ONES SIZE=7
cc: MATRIX_VECTOR_PRODUCT ARG=c1.w,ones
# And output the 7 coordination numbers for the atoms
PRINT ARG=cc FILE=colvar
``` 

In this article I want to explain how we can use PLUMED input files that are similar to the one above to calculate 
other input files.

## CONTACT_MATRIX and the COMPONENTS flag

In the input at the top of this file I used a CONTACT_MATRIX.  This command outputs a matrix, $\mathbf{A}$, in which element $i,j$ is given by:

$$
A_{ij} = \sigma(r_{ij})
$$

where $\sigma$ is a switching function and $r_{ij}$ is the distance between atom $i$ and atom $j$.  

Switching functions like $\sigma(r_{ij})$ are used in the definitions symmetry functions such as those used in the FCCUBIC, TETRAHEDRAL and SIMPLECUBIC actions.  These functions have the following general form:

$$
s_i = \sum_{j\ne i } f(x_{ij},y_{ij},z_{ij})\sigma(r_{ij})
$$

where $f$ is some function of the $x_{ij}$, $y_{ij}$ and $z_{ij}$ components of the vector that connects atom $i$ to atom $j$.  To make it easy to implement such functions I have introduced a 
COMPONENTS flag on the CONTACT_MATRIX command that can be used as shown below:

```plumed
# Calculate the contact matrix and the three matrices that contain the components of the vectors connecting atom i and atom j.  In other words, calculate four 7x7 matrices
c1: CONTACT_MATRIX COMPONENTS GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
# Output the four 7x7 matrices calculated by the command above to a file.
PRINT ARG=c1.w,c1.x,c1.y,c1.z FILE=colvar
```

The COMPONENTS flag for the CONTACT_MATRIX action is specfically designed to make calculating symmetry functions such as $s_i$ striaghforward.  The elements of the matrix `c1.x` are 
thus calculated as:

$$
X_{ij} = \begin{cases}
x_{ij} & \textrm{if} \quad \sigma(r_{ij})>0 \\
0 & \textrm{otherwise}
\end{cases}
$$

The elements of `c1.y` and `c1.z` are also calculated in a similar fashion.  __In other words, the $i,j$ element of these three matrices are only non-zero if atom $j$ is within the first coordination sphere 
of atom $i$.__

## Combining CONTACT_MATRIX and CUSTOM

We can use the CONTACT_MATRIX command with the COMPONENTS keyword to calculate a function similar to $s_i$ above as follows:

```plumed
# Calculate the contact matrix and the three matrices that contain the components of the vectors connecting atom i and atom j.  In other words, calculate four 7x7 matrices
c1: CONTACT_MATRIX COMPONENTS GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
# This applies the function to the three input matrices element-wise and thus outputs a 7x7 matrix
r: CUSTOM ARG=c1.x,c1.y,c1.z FUNC=sqrt(x*x+y*y+z*z) PERIODIC=NO
# Again we are applying the function element wise to the five input matrices and thus outputting a 7x7 matrix
f: CUSTOM ARG=c1.w,c1.x,c1.y,c1.z,r FUNC=w*(x^4+y^4+z^4)/(r^4) VAR=w,x,y,z,r PERIODIC=NO
# We now multiply the 7x7 by a vector of ones to calculate the symmetry function
ones: ONES SIZE=7
s: MATRIX_VECTOR_PRODUCT ARG=f,ones 
# And print the sum of the symmetry function to a file
sums: SUM ARG=s PERIODIC=NO
PRINT ARG=sums FILE=colvar
```

The graph for this input is shown below:

```mermaid
flowchart TB 
MD(positions from MD)
Box("label=Box 
 PBC 
")
Box -- Box --> c1
linkStyle 0 stroke:red,color:red;
MD --> c1
linkStyle 1 stroke:violet,color:violet;
subgraph subc1 [c1]
subgraph subc1_mat [c1]
c1(["label=c1 
 CONTACT_MATRIX 
"])
r(["label=r 
 CUSTOM
FUNC=sqrt(x*x+y*y+z*z) 
"])
f(["label=f 
 CUSTOM
FUNC=w*(x^4+y^4+z^4)/(r^4) 
"])
s(["label=s 
 MATRIX_VECTOR_PRODUCT 
"])
end
style subc1_mat fill:lightblue
sums(["label=sums 
 SUM 
"])
end
c1 -- c1.x --> r
linkStyle 2 stroke:red,color:red;
c1 -- c1.y --> r
linkStyle 3 stroke:red,color:red;
c1 -- c1.z --> r
linkStyle 4 stroke:red,color:red;
c1 -- c1.w --> f
linkStyle 5 stroke:red,color:red;
c1 -- c1.x --> f
linkStyle 6 stroke:red,color:red;
c1 -- c1.y --> f
linkStyle 7 stroke:red,color:red;
c1 -- c1.z --> f
linkStyle 8 stroke:red,color:red;
r -- r --> f
linkStyle 9 stroke:red,color:red;
ones(["label=ones 
 CONSTANT 
"])
f -- f --> s
linkStyle 10 stroke:red,color:red;
ones -- ones --> s
linkStyle 11 stroke:blue,color:blue;
s -- s --> sums
linkStyle 12 stroke:blue,color:blue;
sums -- sums --> 10
10("label=#64;10 
 PRINT
FILE=colvar 
")
```

Notice that all the commands are done in a single chain.  There is no need to store any matrix elements as the functions in `r` and `f` are applied to the elements of the matrices calculated
by `c1` immediately after they are calculated.  Furthermore, if one were to sum the elements of the vector `s` and add a bias upon the sum, the forces on `s` are passed back through the code as
shown below:

```mermaid
flowchart BT 
10(["label=#64;10 
 BIASVALUE 
"])
10 -- sums --> sums
subgraph subc1 [c1]
subgraph subc1_mat [c1]
c1(["label=c1 
 CONTACT_MATRIX 
"])
r(["label=r 
 CUSTOM
FUNC=sqrt(x*x+y*y+z*z) 
"])
f(["label=f 
 CUSTOM
FUNC=w*(x^4+y^4+z^4)/(r^4) 
"])
s(["label=s 
 MATRIX_VECTOR_PRODUCT 
"])
end
style subc1_mat fill:lightblue
c1 -. c1.x .-> r
linkStyle 1 stroke:red,color:red;
c1 -. c1.y .-> r
linkStyle 2 stroke:red,color:red;
c1 -. c1.z .-> r
linkStyle 3 stroke:red,color:red;
c1 -. c1.w .-> f
linkStyle 4 stroke:red,color:red;
c1 -. c1.x .-> f
linkStyle 5 stroke:red,color:red;
c1 -. c1.y .-> f
linkStyle 6 stroke:red,color:red;
c1 -. c1.z .-> f
linkStyle 7 stroke:red,color:red;
r -. r .-> f
linkStyle 8 stroke:red,color:red;
f -. f .-> s
linkStyle 9 stroke:red,color:red;
sums(["label=sums 
 SUM 
"])
s -. s .-> sums
linkStyle 10 stroke:blue,color:blue;
end
sums == sums ==> c1
sums == sums ==> s
subgraph subc1 [c1]
end
s -- ones --> ones
linkStyle 13 stroke:blue,color:blue;
ones(["label=ones 
 CONSTANT 
"])
Box("label=Box 
 PBC 
")
c1 -- Box --> Box
linkStyle 14 stroke:red,color:red;
c1 --> MD
linkStyle 15 stroke:violet,color:violet;
MD(positions from MD)
```

This is possible because the chain of actions the $(i,j)$ matrix elements for `r` and `f` (and their derivatives with respect to the atomic positions) are calculated immediately after $(i,j)$ matrix elements
for `c1.w`, `c1.x`, `c1.y` and `c1.z`.  In other words, in calculating `sums` we do a single loop over $i$ and $j$ and can thus accumulate values and derivatives without storing any matrices or vectors.

## Extending these symmetry functions

There are ways of defining the neighbourhood of atom $i$ that do not involve using the CONTACT_MATRIX keyword.  For example, you can define a matrix in which element $i,j$ is only non-zero if there is a hydrogen
bond between atom $i$ and atom $j$ (HBOND_MATRIX or HBPAMM_MATRIX) or, if you have molecules, you can say that molecule $i$ and $j$ are connected if they are within a cutoff of each other and if the two molecules
have some favourable orientation with respect to each other.  In all these cases you can also use a COMPONENTS keywords to get the direction of the bonds in the first coordination sphere as I have described in the 
article above.  You are thus not forced to use a switching function to define the $\sigma(r_{ij})$ part of the symmetry function that was defined above.  You can use other methods to determine whether atom $i$ and 
atom $j$ are adjacent or not.  Consequently, the implementation of these symmetry function in PLUMED is quite flexible and allows users to try many things without modifying the code.
