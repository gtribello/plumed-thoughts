# Backwards compatibility for MultiColvar

As discussed in [this article](MultiColvar.md), I have made substantial changes to the way MultiColvars are implemented. In previous version the PLUMED
input to calculate the number of distances are less than 0.1 nm is as follows:

```plumed
d1: DISTANCES ATOMS=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS=9,10 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan FILE=colvar
```

Now, however, the same operation is achieved by the following code:

```plumed
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10
d1_lt: LESS_THAN ARG=d1 SWITCH={RATIONAL R_0=0.1}
d1_lessthan: SUM ARG=d1l PERIODIC=NO
PRINT ARG=d1_lessthan FILE=colvar
``` 

Obviously, we do not want users to have to re-learn PLUMED because there have been changes to the code.  In other words, we want to ensure that the first
input above continues to work with the new version of the code.  

## Shortcut actions

It is straightforward to write code that allows users to use the old old PLUMED input syntax.  It is a matter of simply using shortcut actions to read 
the input for actions that are no longer implemented and convert them to the new input.  In other words, and for the example above, the command
DISTANCES is connected to a class that inherits from ActionShortcut.  In the constructor for this object the new input could be created as follows:

```c++
std::string inum, atomstr;  
for(unsigned i=1;;++i) {
    std::string atstring; parseNumbered("ATOMS",i,atstring);
    if( atstring.length()==0 ) break;
    Tools::convert( i, inum ); atomstr += " ATOMS" + inum + "=" + atstring;
}
// Create the action to compute all the distances
readInputLine( getShortcutLabel() + ": DISTANCE" + atomsstr );
// Create the LESS_THAN line
std::string lt_string; parse("LESS_THAN", lt_string );
readInputLine( getShortcutLabel() + "_lt: LESS_THAN ARG=" + getShortcutLabel() + " SWITCH={" + lt_string + "}");
// And finally the sum
readInputLine( getShortcutLabel() + "_lessthan: SUM ARG=" + getShortcutLabel() + "_lt PERIODIC=NO");
```  

This is the essense of what is done in `multicolvar/Distances.cpp` which is where the DISTANCES command is implemented within PLUMED.  The actual code 
is slightly more complicated as DISTANCES allowed for more options than simply LESS_THAN.  However, the key point I am trying to explain here is how 
we can use ActionShortcut objects to implement a single-line version of the three-line input required to calculate the number of distances that are less than 0.1 nm.

## Components of Shortcut actions

Before concluding there is one final thing I need to explain for you to fully understand how these ActionShortcut objects work. There is a wrinkle that is best 
understood by considering the following (old-style) PLUMED input again:

```plumed
# This is a shortcut
d1: DISTANCES ATOMS=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS=9,10 LESS_THAN={RATIONAL R_0=0.1}
# This is not a shortcut
PRINT ARG=d1.lessthan FILE=colvar
```
 
As explained above, once the shortcuts have worked their magic PLUMED will actually read in:

```plumed
# This is the expanded shortcut
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10
d1_lt: LESS_THAN ARG=d1 SWITCH={RATIONAL R_0=0.1}
d1_lessthan: SUM ARG=d1l PERIODIC=NO
# This is not part of the expanded shortcut
PRINT ARG=d1.lessthan FILE=colvar
``` 

The wrinkle that you need to notice here is that the quantity that we are asking PLUMED to print out `d1.lessthan` does not appear to be defined in the above input.  
The quantity that we would like to print out is called `d1_lessthan`.  Obviously, given that the input above works, there is a workaround for this problem and it is this
that still remains to be explained.  

The workaround is that you can define output components for shortcuts just as you would define output components in other ActionWithValue objects.  For our Distances action
above we can thus write a registerKeywords method as follows:

```c++
void Distances::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","ATOMS","the pairs of atoms that you would like to calculate the angles for");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("optional","LESS_THAN","calculate the number of variables that are less than a certain target value.");
  keys.addOutputComponent("lessthan","LESS_THAN","the number of colvars that have a value less than a threshold");
```

When `ActionWithArguments` tries to interpret a Value name such as `d1.lessthan` it now does two things.

* It checks if there is an ActionWithValue object with the label `d1` that has a component called `d1.lessthan`.
* It checks if there is an ActionShortcut object with the label `d1`.  If there is such an ActionShortcut object, PLUMED checks if `lessthan` is amongst the list of registered componets for that action.  If that has indeed been registered then PLUMED interprets `d1.lessthan` as `d1_lessthan` and searches for an action with the label `d1_lessthan`

Notice that this behaviour is also triggered if you use wildcards like this:

```plumed
d1: DISTANCES ATOMS=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS=9,10 LOWEST LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.* FILE=colvar
```

PRINT here would ouptut `d1_lessthan` and `d1_lowest` as these are the two actions that are registerd and created by the DISTANCES shortcut command.

## Conclusion

The new syntax for MultiColvar is (I think) clearer than the old one.  It also more flexible and allows one to do many more things with PLUMED.  I thus want to encourage 
users to learn it and to discard the old-style inputs in future.  To help with this aim, I have ensured that you have the option to see the full input whenever a shortcut
action appears in an example input on the tutorials page and in the nest.  You can try this functionality above.  If you hover over any of the DISTANCES commands in the example 
inputs above you will see an option to expand the shortcut.  If you click this option the full three-line input for `DISTANCES LESS_THAN` is shown.

I think that allowing users to see how code has been implemented in this way is enormously useful. Users can start understanding how the code works by interacting with the rich
set of (automatically-annotated) examples in the nest and tutorials that are likely related to the calculations they wish to perform.  They thus no longer need to battle through 
material in a manual that is full of information that is not relevant to the work they are doing. 

