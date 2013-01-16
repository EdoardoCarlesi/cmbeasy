// class definition for CoupledInvariant
// version 1.0
#ifndef COUPLEDINVARIANT_H
#define COUPLEDINVARIANT_H 

#include "invariant.h"

class QuintCosmos;
class OriginalCoupledQuintCosmos;
class Quintessence;

/*! This class calculates the perturbations for a CoupledQuintCosmos cosmology. It uses invariant variables,
the quintessence field variables have been formed into Dgq and Vq for convenience, since this gives a better
impression in terms of energy density perturbations etc. It is assumed that the field-dependent coupling
splines have been set in CoupledQuintCosmos, which is automatic if history() is called */

class CoupledInvariant : public Invariant {
 public:
  CoupledInvariant(QuintCosmos* c);
  virtual ~CoupledInvariant() {};

  double mona_a; // Mona

  //Georg: Christian and Mona forgot to add 2 for the Nr of Perturbations in calcPerturbationNr ->diff of about 1.5% at l=1400
  static bool UseWrongNumberOfPertEquationLikeOriginal;

  /*! This performs the time steps for each k mode. 
    \param tau time at which we like to evaluate the perturbation equations
    \param y present values
    \param yprime derivates of the respective values*/
  virtual void fderivs(const double tau, const double *y, double *yprime);

  virtual void getReady(const ControlPanel& control);

  virtual void initialScalarPerturbations(const ControlPanel &control, const double tau);

  virtual void fillPhiDot(double tau); 
  virtual void calcPerturbationNr(const ControlPanel &control);

  virtual void scalarSources(double tau, double *d, double *dp, double *dk); // added by Mona!!!!!!!!!
 
 protected:
  OriginalCoupledQuintCosmos *quintcosmos;
  Quintessence* quint;
  int qidx; // the index position of the delta phi and delta phidot variable in the y[] array
};

#endif 
