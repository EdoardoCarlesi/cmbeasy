#ifndef TENSORTRAPEZOIDINTEGRATOR_H
#define TENSORTRAPEZOIDINTEGRATOR_H

#include "tensorintegrator.h"

/*! 
  Integrator for tensor C_l's. Convolutes the sources
  with the bessel functions. Then performs the final k-integration.
*/

class TensorTrapezoidIntegrator : public TensorIntegrator {
 public:
  TensorTrapezoidIntegrator(Cosmos *c,double,double,double,SplineWeb*,SplineWeb*,SplineWeb*); 
   virtual ClReturn integrate(int l, Spline *bessel,const ControlPanel&);
};

#endif
