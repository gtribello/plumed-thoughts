# Passing atoms to/from PLUMED

As discussed in [this article](MDInterfaceI.md), data is passed from the MD code to PLUMED by creating a PUT action.
These PUT actions take the data from a void pointer that is passed to PLUMED from the MD code and transfer it into a 
PLMD::Value object. Passing a void pointer and using a PUT action to convert the data within it
to a PLMD::Value is also used when the atomic positions, masses and charges are sent to PLUMED. However, 
there are some other subtleties for these quantities because MD codes use a domain decomposition and scatter the properties of atoms across
multiple domains. We thus need to use the action DOMAIN_DECOMPOSITION when passing these quantities to make sense of the 
information in the void pointers that the MD code passes.

## Creating a DOMAIN_DECOMPOSITION action

A DOMAIN_DECOMPOSITION can be created by using a call to plumed.cmd as follows:

```c++
plumed.cmd("readInputLine dd: DOMAIN_DECOMPOSITION NATOMS=20 VALUE1=vv UNIT1=length PERIODIC1=NO CONSTANT1=False");
```

The DOMAIN_DECOMPOSTION command above creates a PUT action with the label vv. The pointer to the data that needs to be transferred to the PLMD::Value
object that is created by the PUT action is then set by using the command below:

```c++
plumed.cmd("setInputValue vv, &val);
```

Meanwhile, the pointer to the forces that should be modified is passed as follows:

```c++
plumed.cmd("setValueForces vv", force);
```

In other words, pointers to values and forces in the MD code are passed to PUT actions that are created by the DOMAIN_DECOMPOSION in 
[the way you pass data to other PUT actions](MDInterfaceI.md). 

The PLMD::Value objects created by a DOMAIN_DECOMPOSITION action are always vectors with the same number of components as atoms in the system. Furthermore, you can create multiple PUT
actions from a single DOMAIN_DECOMPOSITION action. To see why this is useful, consider the following DOMAIN_DECOMPOSITION command:

```plumed
gromacs: DOMAIN_DECOMPOSITION ...
   NATOMS=2000
   VALUE1=posx UNIT1=length PERIODIC1=NO CONSTANT1=False ROLE1=x
   VALUE2=posy UNIT2=length PERIODIC2=NO CONSTANT2=False ROLE2=y
   VALUE3=posz UNIT3=length PERIODIC3=NO CONSTANT3=False ROLE3=z
   VALUE4=Masses UNIT4=mass PERIODIC4=NO CONSTANT4=True ROLE4=m
   VALUE5=Charges UNIT5=charge PERIODIC5=NO CONSTANT5=True ROLE5=q
...
```

This action is created when you call `plumed.cmd("setNatoms",&natoms)` from gromacs. It makes 5 PLMD::Value called posx, posy, posz, Masses and Charges. 
These PLMD::Value objects then hold the x, y and z positions of the atoms and the masses and charges of the atoms. It is important to note that this command will 
also, create a PBC_ACTION to hold the cell.

The ROLE keywords above are only needed because the five quantities passed by the command above play a unique role within PLUMED. If you pass 
some other quantities, this instruction is not required. PLMD::ActionAtomistic searches for atomic positions, masses and charges by looking for PUT Actions
that have these five particular roles and for ActionWithVirtualAtom objects.

## Differences from regular PUT actions

PUT actions created from a DOMAIN_DECOMPOSITION action behave differently from other PUT actions. In particular:

* If a DOMAIN_DECOMPOSITION action creates a PUT action, then the PUT action depends on the DOMAIN_DECOMPOSITION action. ActionToPutData::apply thus does nothing for these PUT actions.
* Similarly, when DOMAIN_DECOMPOSITION actions create PUT actions, data is transferred from the input pointer to the PLMD::Value object by the DOMAIN_DECOMPOSITION action. When ActionToPutData::wait is called for these PUT Actions `wasset=true`, ActionToPutData::wait does nothing.
* Lastly, if a constant PUT action is created by DOMAIN_DECOMPOSITION, the values in the vector are set during the first step of MD rather than during initialisation. 

These differences are necessary because PUT actions cannot understand the information in the pointers that are passed from the MD code. These pointers are only understandable if you know 
which atoms are on each processor. This information is only passed to the DOMAIN_DECOMPOSITION action. DOMAIN_DECOMPOSITION must translate the information passed from the MD code before it is 
passed back on through PLMD::Value objects created by the PUT actions. DOMAIN_DECOMPOSITION thus keeps pointers to all the PUT actions that it creates. It sets the data in these action's PLMD::Value objects
within `DomainDecomposition::share(const std::vector<AtomNumber>& unique)` and `DomainDecomposition::wait().`  The forces on the PUT actions created by the DOMAIN_DECOMPOSITION action are added in `DomainDecomposition::apply()`.

