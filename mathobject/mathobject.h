#ifndef MATHOBJECT
#define MATHOBJECT



/*! 
  Base class of all classes that wish to use
  basic numerical methods. The reason for 
  this class to exist is the fact that if an object
  wants to call a function that is not part of the 
  object (say some integration routine) which
  in turn needs to be given some function belonging
  to the object and using data of this object,
  then both the object itself and the integration
  routine need to be member of one and the
  same class. This is similar to the QObject 
  strategy of the Qt library. 
  
  For convenience, two typedefs are made:
  moDerivs is a shorthand notation for a function
  of the type Miscmath::odeint() needs.
  moSingle is just a function of a double value
  returning double.

  For a nice example, see the Spline class itself,
  where the integrate() method calls Miscmath::odeint().

  Of course, you always have to give a reference
  to the object that calls some Miscmath routine,
  if you want the Miscmath routine to use a function
  of your object which you specify.
  If not, how should the Miscmath routine know 
  which object to ask for the data / functions ?  
*/
class  Mathobject {};

typedef void (Mathobject::*moDerivs) (const double, const double*, double*) const;  //! common derivs for Nrecipes ODE
typedef double (Mathobject::*moSingle) (const double) const; //! Simple f(x)   

#endif
