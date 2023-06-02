# Representing PLUMED inputs using flowcharts

Ever since we wrote the PLUMED paper in 2014, I have wondered if we can automate the process of generating diagrams
similar to the ones in that paper that show how Values and forces are passed between the various actions in the PLUMED input.
If we could automate the generation of these diagrams, we could have similar charts for all the entries 
in the PLUMED nest. We could also have these diagrams for the example inputs in the tutorials and manuals.

The time to implement these diagrams is now when we have changed the code so that actions can also pass vectors,
scalars and grids between them, as discussed [here](Passing.md). I have thus written a command line tool that allows one to 
generate a graph from a plumed input. Consequentially, if I take the contents of the file plumed.dat below:

```plumed
c1: COM ATOMS=1-10
c2: COM ATOMS=11-20
d1: DISTANCE ATOMS=c1,c2 COMPONENTS
r: RESTRAINT ARG=d1.z AT=1 KAPPA=1 
f1: CUSTOM ARG=d1.x,d1.y FUNC=x*x+y*y PERIODIC=NO
PRINT ARG=d1.x,d1.y,f1,r.bias FILE=colvar
``` 

And run the command:

````
plumed show_graph --plumed plumed.dat --out graph.md
````

I can generate the graph shown below:

```mermaid
flowchart TB
MD{{positions from MD}}
Box{{"`label=Box 
 PBC
`"}}
Box -- Box --> c1
linkStyle 0 stroke:red,color:red;
MD --> c1
linkStyle 1 stroke:violet,color:violet;
c1(["label=c1
 COM"])
Box -- Box --> c2
linkStyle 2 stroke:red,color:red;
MD --> c2
linkStyle 3 stroke:violet,color:violet;
c2(["label=c2
 COM"])
Box -- Box --> d1
linkStyle 4 stroke:red,color:red;
c1 -- c1 --> d1
linkStyle 5 stroke:violet,color:violet;
c2 -- c2 --> d1
linkStyle 6 stroke:violet,color:violet;
d1(["label=d1
 DISTANCE"])
d1 -- d1.z --> r
r(["label=r
 RESTRAINT"])
d1 -- d1.x --> f1
d1 -- d1.y --> f1
f1(["label=f1
 CUSTOM
FUNC=x\*x+y\*y"])
d1 -- d1.x --> 7
d1 -- d1.y --> 7
f1 -- f1 --> 7
r -- r.bias --> 7
7("label=#64;7
 PRINT
FILE=colvar")
```

The file `graph.md` output by the command above is renderable using [mermaid](https://mermaid.js.org/syntax/flowchart.html). You can see the resulting flow chart if you copy and paste the file's contents 
[here](https://mermaid.live/). I used Mermaid to build the charts, as you can insert Mermaid syntax into GitHub markdown. The rendered diagrams then
appear when GitHub shows the rendered markdown online.

Each node in the diagram above represents one of the actions from the PLUMED input file. The arrows then indicate how PLMD::Value
objects are passed between the actions.  

The shape of the node tells you about the type of action:

* Irregular-hexagonal nodes are PUT actions containing data passed from the MD code. These nodes cannot take PLMD::Value objects created in PLUMED as input.
* Rectangular nodes are actions like PRINT that only take PLMD::Value as arguments. These nodes cannot create PLMD::Value objects and pass them to other actions.
* Rounded nodes are actions that can take PLMD::Value objects created within PLUMED as input and pass on such objects as output.

The arrows connecting the actions provide information about the PLMD::Value object being passed.

* Passing of __scalars__ is indicated using __black arrows__ 
* Passing of __vectors__ is indicated using __blue arrows__
* Passing of __matrices__ is indicated using __red arrows__
* Passing of __grids__ is indicated using __greeen arrows__
* Passing of __atomic positions__ is indicated using __violet arrows__.  An atomic position is just five PLMD::Value objects that are all vectors. These five vectors contain the x, y and z positions of the atoms and the masses and charges of the atoms.

You can also show how forces are passed between actions by using the command:

````
plumed show_graph --plumed plumed.dat --out graph.md --force
```

When I run the command above on the plumed input above, I obtain the following flowchart:

```mermaid
flowchart BT
r(["label=r 
 RESTRAINT"])
r -- d1.z --> d1
d1(["label=d1
 DISTANCE"])
c2(["label=c2
 COM"])
c1(["label=c1
 COM"])
Box{{"label=Box
 PBC"}}
c1 -- Box --> Box
linkStyle 1 stroke:red,color:red;
c1 --> MD
linkStyle 2 stroke:violet,color:violet;
c2 -- Box --> Box
linkStyle 3 stroke:red,color:red;
c2 --> MD
linkStyle 4 stroke:violet,color:violet;
d1 -- Box --> Box
linkStyle 5 stroke:red,color:red;
d1 -- c1 --> c1
linkStyle 6 stroke:violet,color:violet;
d1 -- c2 --> c2
linkStyle 7 stroke:violet,color:violet;
MD{{positions from MD}}
```

Notice that fewer actions are shown in this new graph. This is because the graph above only shows actions that play some role in the force calculation.

I have found these diagrammatic representations of PLUMED input files enormously beneficial when dealing with complicated PLUMED input files. I will thus use them extensively in these notes about the work that I have done in refining PLUMED.
