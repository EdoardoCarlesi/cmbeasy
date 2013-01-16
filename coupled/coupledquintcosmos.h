#ifndef coupledquintcosmos_H
#define coupledquintcosmos_H 

#include "global.h"
#include "spline.h"
#include <iostream>
#include "quintessence.h"
#include "anchor.h"
#include "safevector.h"
#include "controlpanel.h"
#include "cl.h"
#include "coupling.h"

#include "cosmos.h"
#include "quintcosmos.h"

class CoupledQuintCosmos : public QuintCosmos {

 public: 
  CoupledQuintCosmos(Quintessence::Type = Quintessence::none);
  ~CoupledQuintCosmos();
  
  void normExtraSplines();

  Coupling* couple;
  void setCoupling(Coupling* c){couple=c;}
  Coupling* coupling() const { return couple; }
  
  virtual void propagateHistoryInTau(const double , const double*,double *); 
  void history(bool inform=false); 
  virtual void initializeDensities(double* y, double* ydot); 
};
#endif
