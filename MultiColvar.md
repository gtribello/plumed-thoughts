# Reimplementing MultiColvar

The original idea for the MultiColvar class came from a conversation I had with Fabio Pietrucci.
The line of thought went something like this:

* When calculating many collective variables you calculate the same function $f$ for multiple sets of atomic positions $\{X\}_i$.
* The sum of the functions is then computed. The final CV is thus:

$$
s = \sum_{i=1}^N f(\{X\}_i)
$$

* As it is straightforward to parallelise the sum in the above expression, why should we not write a base class for calculating CVs like this?
* If this base class is implemented, developers can inherit from it and only need to implement the $f(\{X\})$ part of the expression above.

The MultiColvarBase class I implemented after this conversation was designed with this line of thinking in mind. With this class, you can thus use
a command like:

```plumed
d1: DISTANCES ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan FILE=colvar
```
to calculate and print the number of the five distances above that are less than 0.1 nm.

A problem with the code I wrote quickly emerged. Many users and developers wanted access to the full vector of $f(\{X\})_i$ values rather than access to sums of these values. 
I had to add code into the base class to give them access to this vector of values. This tinkering complicated the MultiColvarBase class, making using or developing features 
where the MultiColvarBase class was involved difficult. Much of the rewriting I have done has
 aimed to simplify the MultiColvarBase class and resolve these issues. This is why I have added functionality
to pass vectors between actions. Giving users and developers direct access to the vectors reduces the amount of code that needs to be in the 
MultiColvarBase class. This functionality can be moved to other actions, and complicated CVs can be implemented directly from the input file (or using shortcuts).
For example, The number of distances less than 0.1 nm that I calculated above is now calculated using the following input:

```plumed
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10
d1l: LESS_THAN ARG=d1 SWITCH={RATIONAL R_0=0.1}
d1s: SUM ARG=d1l PERIODIC=NO
r: RESTRAINT ARG=d1s KAPPA=1 AT=3
```  

If you look at the flowchart representation for this input, you can understand how it works more clearly:

```plumed
#MERMAID=value
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10
d1l: LESS_THAN ARG=d1 SWITCH={RATIONAL R_0=0.1}
d1s: SUM ARG=d1l PERIODIC=NO
r: RESTRAINT ARG=d1s KAPPA=1 AT=3
```

The first action above calculates the five distances and passes a vector with five elements to the LESS_THAN action that follows it. This LESS_THAN action
then applies a function elementwise to the five components of the vector. The output PLMD::Value object from this action is thus another vector. This vector 
is then converted to a scalar by the SUM action, which adds all the elements of the vector together.

The derivatives present a problem when implementing CVs using the method described above. The derivative for a vector of CVs is a matrix. For complicated 
CVs that depend on the positions of many atoms, this matrix quickly gets too large to be stored. I have resolved this problem by 
calculating the CV once during the forward (calculate) loop and a second time during the backwards (apply) loop.  This turns out to be much faster than what was
done in the older version of PLUMED where the derivative of the sum with respect to the atoms is accumulated in the forward loop by using the chain rule.  The flowchart 
representations that describe how forces are applied are thus similar to those that show how the forces are applied on scalar quantities.

```plumed
#MERMAID=force 
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10
d1l: LESS_THAN ARG=d1 SWITCH={RATIONAL R_0=0.1}
d1s: SUM ARG=d1l PERIODIC=NO
r: RESTRAINT ARG=d1s KAPPA=1 AT=3
```

Exposing the vector of CV values calculated by a MultiColvar by putting them in the output  PLMD::Value object has dramatically simplified this base class. I have been able 
to write the base class of a MultiColvar as a template. If you have written a method that inherits from Colvar, you can thus also write an action to calculate a vector that contains multiple 
instances of your CV as well as the usual scalar-valued CV by adding the following lines to your Colvar's cpp file:

```c++
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"

// Create a shortcut action that determines if the output PLMD::Value will be a scalar or a vector
typedef ColvarShortcut<Distance> DistanceShortcut;
// Ensure that the shortcut is created when the keyword DISTANCE appears
PLUMED_REGISTER_ACTION(DistanceShortcut,"DISTANCE")
// Register your scalar-valued CV as <name>_SCALAR
PLUMED_REGISTER_ACTION(Distance,"DISTANCE_SCALAR")
// Create a class to compute multiple instances of your CV
typedef MultiColvarTemplate<Distance> DistanceMulti;
// Reister the vector of CVs as <name>_VECTOR
PLUMED_REGISTER_ACTION(DistanceMulti,"DISTANCE_VECTOR")
```

You must implement Colvars slightly differently if you want to use the above functionality. In particular, you need to write three static methods:

* `void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa )` reads the atoms that are used to calculate the CV by parsing the keywords that have type `atom`.
* `unsigned getModeAndSetupValues( ActionWithValue* av )` creates the named components passed between PLUMED actions.
* `void calculateCV( const unsigned& mode, const std::vector<double>& masses, const std::vector<double>& charges, const std::vector<Vector>& pos, std::vector<double>& vals, std::vector<std::vector<Vector> >& derivs, std::vector<Tensor>& virial, const ActionAtomistic* aa )` calculates the CV.

As you can see by looking at what was done in `colvar/Distance.cpp`, for example, this is not particularly hard. Implementing these methods is usually a matter of moving code that you would have to write 
in the constructor of calculate methods of the Colvar to the functions above. Furthermore, you can still inherit from Colvar and implement a single, scalar-valued CV as you did in the past.

Similar changes were also required to deal with Function actions, as these actions can now have scalar or vector input and output. For functions, I implemented two template classes 
`FunctionOfScalar` and FunctionOfVector`.  The template parameters for these classes should inherit from `FunctionTemplateBase` and an action that applies a switching function act 
on a scalar or vector can then be written as follows:

```c++
#include "FunctionTemplateBase.h"
#include "tools/SwitchingFunction.h"
#include "FunctionShortcut.h"
#include "FunctionOfVector.h"
#include "ActionRegister.h"

namespace PLMD {
namespace function {

class LessThan : public FunctionTemplateBase {
  bool squared;
  SwitchingFunction switchingFunction;
public:
  void registerKeywords( Keywords& keys ) override;
  void read( ActionWithArguments* action ) override;
  bool getDerivativeZeroIfValueIsZero() const override { return true; }
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};

typedef FunctionShortcut<LessThan> LessThanShortcut;
PLUMED_REGISTER_ACTION(LessThanShortcut,"LESS_THAN")
typedef FunctionOfVector<LessThan> VectorLessThan;
PLUMED_REGISTER_ACTION(VectorLessThan,"LESS_THAN_VECTOR")
  
void LessThan::registerKeywords(Keywords& keys) {
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present, you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.addFlag("SQUARED",false,"is the input quantity the square of the value that you would like to apply the switching function to");
}

void LessThan::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) action->error("should only be one argument to less_than actions");
  if( action->getPntrToArgument(0)->isPeriodic() ) action->error("cannot use this function on periodic functions");
  
  string sw,errors;
  action->parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
    if( errors.length()!=0 ) action->error("problem reading SWITCH keyword : " + errors );
  } else {
    int nn=6; int mm=0; double d0=0.0; double r0=0.0; action->parse("R_0",r0);
    if(r0<=0.0) action->error("R_0 should be explicitly specified and positive");
    action->parse("D_0",d0); action->parse("NN",nn); action->parse("MM",mm);
    switchingFunction.set(nn,mm,r0,d0);
  }
  action->log<<"  using switching function with cutoff "<<switchingFunction.description()<<"\n";
  action->parseFlag("SQUARED",squared);
  if( squared ) action->log<<"  input quantity is square of quantity that switching function acts upon\n";
}

void LessThan::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  plumed_dbg_assert( args.size()==1 );
  if( squared ) vals[0] = switchingFunction.calculateSqr( args[0], derivatives(0,0) );
  else vals[0] = switchingFunction.calculate( args[0], derivatives(0,0) );
  derivatives(0,0) = args[0]*derivatives(0,0);
}

}
}
```

As you can see, the key methods here are `read`, which reads the parameters of the function from the input line, and `calc`, which calculates the value of the 
function and the derivative at the point specified in `args`.
 
As with Colvars, the code above is just a reordering of what you would have done in the old version of PLUMED. Furthermore, you don't need to implement 
functions using the new method outlined above. You can still inherit from `Function` and implement a scalar-valued function 
that takes scalar arguments only in the way that you always did
