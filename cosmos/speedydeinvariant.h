#ifndef SPEEDYDEINVARIANT_H
#define SPEEDYDEINVARIANT_H 

#include "speedyinvariant.h"

/*!
  
  Gauge invariant perturbation class implementing the speed ups of
  astro-ph/0503277

  The dark energy is modelled as a scalar fluid with rest-frame
  speed of sound mu^2 which you can change in fderivs().

  The propagateScalar() routine is a complete rewrite and re-implentation
  of the one from the Perturbation parent class.

*/

class QuintCosmos;

class SpeedyDEInvariant : public SpeedyInvariant {
 protected:
  QuintCosmos *quintcosmos;
  int qidx; // the index position of the delta phi and delta phidot variable in the y[] array
  int xidx; //!< the index of the field perturbation derivatives
  bool IsPhantomCrossing;

 public:
  SpeedyDEInvariant(QuintCosmos* c);
  virtual ~SpeedyDEInvariant() { if (false) cout << "called: " << called << " thereof called2: " << called2 << " called3: " << called3 << endl; }

  virtual void fderivs(const double tau, const double *y, double *yprime);
  virtual void initialScalarPerturbations(const ControlPanel &control, const double tau);

  virtual void fillPhiDot(double tau); //!< calculate d/dtau of phi and psi
 
 
  virtual void calcPerturbationNr(const ControlPanel &); //! calculate the number of scalar and tensor equations 
  //void foutputt_(const ControlPanel& control,int *n, double *y, double *ypr,  double tau, double *dt, double *dte, double *dtb);
 
  virtual void getReady(const ControlPanel&);
};

#endif 
