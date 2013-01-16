#ifndef ARBITRARYOMEGA_H
#define ARBITRARYOMEGA_H

#include "arbitrary.h"

/*!
  Arbitrary Omega_dark_energy(z) class
  Inherit from this class and implement omega(z) to
  get a scalar field dark energy evolution that yields the
  wanted Omega_de(z).
*/

class ArbitraryOmega : public Arbitrary
{
 public:
  ArbitraryOmega(QuintCosmos &c) : Arbitrary(c) {}

  //! Omega_de as a function of redshift
  virtual double omega(double z) {
    throw Bad_Error("ArbitraryOmega::omega - implement omega(z)");
  }
  virtual Type type() { return arbitraryomega; }
  virtual void prepare();

  virtual double w(double a) const  { return W->fastY(-log(a));}  //!< equation of state as function of scale factor
  virtual double dwda(double a) const { return -1.0/a*DW->fastY(-log(a));} //!< first  derivative wrt a:  d w(a)/ da
  virtual double d2wda2(double a) const { return (DDW->fastY(-log(a))  + DW->fastY(-log(a)) )/ (a*a);}  //!< second derivative wrt a:  d^2 w(a)/ da^2

 protected:
  Anchor anchor;
  Spline *W,*DW,*DDW,*K,*Omega,*DOmega;
};


#endif
