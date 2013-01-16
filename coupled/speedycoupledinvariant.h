#ifndef SPEEDYCOUPLEDINVARIANT_H
#define SPEEDYCOUPLEDINVARIANT_H 

#include "speedydeinvariant.h"
#include "coupledinvariant.h"
#include "coupledquintcosmos.h"

/*!
  Gauge invariant perturbation class implementing the speed ups of
  astro-ph/0503277, modiefied for coupling between cold dark matter
  and dark energy.

  The dark energy is modelled as a scalar fluid with rest-frame
  speed of sound mu^2 which you can change in fderivs().

  The propagateScalar() routine is a complete rewrite and re-implentation
  of the one from the Perturbation parent class.
*/

class QuintCosmos;

class SpeedyCoupledInvariant : public SpeedyDEInvariant {
 protected:
  QuintCosmos *quintcosmos;
  int qidx; // the index position of the delta phi and delta phidot variable in the y[] array
  bool IsPhantomCrossing;
 public:
  SpeedyCoupledInvariant(QuintCosmos* c);
  virtual ~SpeedyCoupledInvariant() { 
 if (false) cout << "called: " << called << " thereof called2: " << called2 << " called3: " << called3 << endl; 
  }
	
  virtual void fderivs(const double tau, const double *y, double *yprime);
  virtual void initialScalarPerturbations(const ControlPanel &control, const double tau);
  virtual void fillPhiDot(double tau); //!< calculate d/dtau of phi and psi
  virtual void calcPerturbationNr(const ControlPanel &); //! calculate the number of scalar and tensor equations 
  //void foutputt_(const ControlPanel& control,int *n, double *y, double *ypr,  double tau, double *dt, double *dte, double *dtb);

  virtual void getReady(const ControlPanel&);
};

#endif 
