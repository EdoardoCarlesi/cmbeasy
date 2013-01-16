#ifndef TENSORINTEGRATOR_H
#define TENSORINTEGRATOR_H

#include "integrator.h"
#include <map>

class Cosmos;
class Spline;
class SplineWeb; 
/*! 
  Integrator for tensor C_l's. Convolutes the sources
  with the bessel functions. Then performs the final k-integration.
*/

class TensorIntegrator : public Integrator {
 protected:
  SplineWeb *Dt, *Dte, *Dtb;   //!< tensor sources web
  double stpt; // from k * tau > stpt on, the sources are identical 0
  vector<Spline*> source,sourceE,sourceB; // the (interpolated in k) sources coming from the dtWeb
  
 public:
  TensorIntegrator(Cosmos *c,double,double,double,SplineWeb*,SplineWeb*,SplineWeb*); 
  ClReturn integrate(int l, Spline *bessel,const ControlPanel&);
};

#endif
