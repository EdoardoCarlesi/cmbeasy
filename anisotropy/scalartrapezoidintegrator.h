#ifndef SCALARTRAPEZOIDINTEGRATOR_H
#define SCALARTRAPEZOIDINTEGRATOR_H

#include "scalarintegrator.h"

/*! 
  Integrator for scalar fluctuations. It folds the sources with
  the bessel functions. It does so in quite the same manner
  as the original CMBFAST routine (except that it monitors
  the sum of the integration to stop earlier -- if possible, 
  uses different  precisions for different k-values and uses Splines). 

  However, it should be emphasized that one may come up
  with an entirely different integration strategy. Unfortunately,
  the [IMHO clever] things I tried turned out not to beat the original 
  CMBFAST strategy, which is why this ScalarIntegrator sticks
  with the slightly modified original one.

  So if you come up with some way to do this job in a better fashion 
  -- maybe even with error - estimates -- please feel 
  encouraged to write a substitute.
*/

class ScalarTrapezoidIntegrator : public ScalarIntegrator {
 public:
  ScalarTrapezoidIntegrator(Cosmos *c,double,double,SplineWeb* d, SplineWeb *dp);
  ~ScalarTrapezoidIntegrator();

  virtual ClReturn integrate(int l, Spline * bessel,const ControlPanel&);   
};

#endif
