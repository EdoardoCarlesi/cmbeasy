#ifndef basecosmos_H
#define basecosmos_H

#include "global.h"
#include "mathobject.h"

/*!
  Abstract base class for cosmology background evaluators.
  The reason is to allow alternative cosmology classes that
  don't need all of the splines and other functions that
  cosmos carries around -- for example, if you aren't including
  CMB calculations.

  Any derived class of basecosmos must contain a method
  to determine its history, and to return various distance
  measures.
*/ 
class baseCosmos : public Mathobject {
 protected:

  bool ValidHistory; //!< true if history() has filled interpolation tables and calculated sucessfully the history of the universe
  bool isRebounding; //!< true if the Universe did not experience a big bang, but is instead rebounding from a previous bout of cosmic contraction.  Only valid if history() has been called.

 public:

  baseCosmos() : ValidHistory(false) { }
    virtual ~baseCosmos() {}
 
  //Plank Mass           
  static double M_p() {return 3.809173e+56;} //!< Reduced Plank Mass in Mpc^-1

  virtual void history(bool inform=false) = 0; //!< Calculates history of cosmos
  bool validHistory() const { return ValidHistory; }  //!< history() sets this flag  

  bool isUniverseRebounding() const { return isRebounding; }; //!< Is Universe rebounding from cosmic contraction? Only valid if history() has been called

  //Scale factor stuff
  virtual double a_0() const = 0; //!< scale factor today
  double a2z(double a) const { return a_0() / a  - 1;}  //!< redshift as a function of a
  double z2a(double z) const { return a_0()/(1+z);} //!< scale factor as a function of redshift

  //Hubble constant as a function of redshift
  virtual double h() const = 0; //!< Small Hubble h today
  virtual double Z2iH(double z) = 0; //!< Translates between redshift and inverse Hubble parameter
  virtual double Z2H(double z) { return 1.0/Z2iH(z); }  //!< Translates from redshift to Hubble parameter

  //Basic Omega stuff
  virtual double omega_m() = 0; //!< Energy density in non-relativistic matter
  virtual double omega_0() = 0; //!< Total omega today (i.e., everything but curvature)
  double omega_k() { return 1 - omega_0(); } //!< Omega curvature
  double omega_curv() { return omega_k(); } //!< Omega curvature ( synonym for omega_k() )

  //Distances
  virtual double propermotionDistance(double) const = 0; //!< returns proper motion distance
  virtual double angulardiameterDistance(double) const = 0; //!< return angular diameter distance
  virtual double luminosityDistance(double) const = 0; //!< return luminosity distance, ignoring difference between heliocentric and CMB frame redshifts
  virtual double luminosityDistance(double,double) const = 0; //!< return luminosity distance

  //Distance derivatives
  virtual double propermotionDistanceDeriv(double) const = 0; //!< return derivative of propermotion distance with respect to redshift
  virtual double angulardiameterDistanceDeriv(double) const = 0; //!< return derivative of angular diameter distance with respect to redshift
  virtual double luminosityDistanceDeriv(double) const = 0; //!< return derivative of luminosity distance with respect to redshift

};

#endif
