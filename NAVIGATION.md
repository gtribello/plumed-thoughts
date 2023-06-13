# Revising PLUMED 

After the Trieste meeting in 2017, I (Gareth Tribello) started revising some of the code in PLUMED in a separate branch called hack-the-tree.
The code in this new branch quickly diverged from the code in the master branch. This divergence, I now realise, was a good thing. The extensive
code base for a piece of scientific software places innumerable constraints on its development. Sometimes it helps to free yourself of these
constraints and to allow yourself to "break" some things. I think this process has allowed me to improve how some of the CVs
and methods I have written for PLUMED are implemented.

I think it is now time to merge some of these features from the hack-the-tree branch back into the master branch of PLUMED. I decided to write these pages as
I made the code from hack-the-tree amenable to the rest of the code. I hope they explain my thinking and help others use the
new features I have implemented. I also hope they inspire discussions about how we approach code development for a particular research community.
When I started this process I identified the following three themes that have driven a lot of my thinking over the last three years:

* [Passing data between actions](Passing.md) 
* [Reproducibility and extensibility](Reproducibility.md)
* [Hierachy and community](Community.md)

A lot of the posts in the pages that follow expand on these three ideas.


```mermaid
flowchart TB;
  A[Passing data] -.-> B[Reproducibility];
  B -.-> C[Hierachy];
  C -.-> D[Passing data to and from PLUMED I];
  D --> E[Passing data to and from PLUMED II];
  A --> D;
  E --> F[Passing the energy];
  A --> G[Virtual atoms]
  B --> H[Graphs]
  A --> I[multi colvar]
  H --> I
  I --> J[multi colvar shortcuts]
  B --> J
  C --> J
  J --> K[contact matrices]
  I --> K
  click A "Passing.md" "8th May 2023: General thoughts about how data is passed between PLUMED actions";
  click B "Reproducibility.md" "8th May 2023: General thoughts about why we want to do calculations that are reproducibile";
  click C "Community.md" "8th May 2023: General thoughts about how we support communities of scholars"
  click D "MDInterfaceI.md" "8th May 2023: A description of how data is passed to and from PLUMED"
  click E "MDInterfaceII.md" "21st May 2023: A description of how atomic properties are passed to and from PLUMED"
  click F "PassingEnergy.md" "21st May 2023: A description of how potential energy is passed to PLUMED"
  click G "VirtualAtoms.md" "22nd May 2023: A description of how virtual atom positions are passed between Actions in PLUMED"
  click H "Graphs.md" "1st June 2023: A description of the way we can use graphs to illustrate PLUMED input files"
  click I "MultiColvar.md" "1st June 2023: An explaination of how multicolvar is implemented"
  click J "MultiColvarShortcut.md" "9th June 2023: An explanation of how backwards compatibility for MultiColvar has been ensured and an introduction to ActionShortcut"
  click K "contactMatrix.md" "12th June 2023: An explanation of how coordination numbers and contact maps can be used to construct CVs"
```
