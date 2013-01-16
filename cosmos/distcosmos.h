#ifndef dist_cosmos_H
#define dist_cosmos_H

#include <cmath>
#include <iostream>

#include "basecosmos.h"

#include "spline.h"
#include "anchor.h"
#include "miscmath.h"

/*!
  Limited implementation of cosmology class which
  only contains enough information to calculate
  the distance measures, and not CMB or power spectra.

  Like cosmos, all quantities are in powers of Mpc.

  A distCosmos has a pretty basic set of cosmological parameters.
  The only ones considered are Omega_matter, Omega_vacuum,
  and Omega_relativistic particles.  It isn't smart enought to
  handle anything that is relativistic over only part of the redshift
  interval.

  The background evolution is calculated by calling history().

  Access to the background quantities is particularly simple by  using 
  functions that are named  
  \code
  y = x2y(x);
  \endcode

  where some quantity y is calculated as a function of x. 
*/

class distCosmos : public baseCosmos {
protected:
  static const double flattol; //!< Tolerance for flatness.  If abs(Omega_k) is less than this, the Universe is treated as flat
  static const int ndist; //!< Number of luminostiy distance steps
  double zmax; //!< Maxiumum redshift calculated

  double h_h; //!< Current small Hubble h

  double Omega_m; //!< Omega_matter (non-relativistic)
  double Omega_v; //!< Omega_Vacuum
  double Omega_rel; //!< Omega_relativistic
  double Omega_k; //!< Omega_curvature

  Miscmath::keyword odetype; //!< Type of ode integral to do.  Options are rkutta and bstoer

  Anchor SplineAnchor; //!< this is the distCosmos' Anchor -- it takes care of the splines and is asked by reset() to get rid of the all splines to avoid memory leaks
  Spline* distSpline; //!< Stores distance measure
  Spline* dDdzSpline; //!< Derivative of distance measure with respect to z
  Spline* Z2iHSpline; //!< Map from z to inverse Hubble parameter

  virtual bool calc_isrebounding(); //!< Does calculation to determine if Universe is rebounding

  void iEode(double, const double*, double*) const; //!< Integration wrapper for 1/E

public :
  
  distCosmos(); //!< standard constructor
  virtual ~distCosmos() {};

  void reset(); //!< resets() distCosmos. Frees all Splines and calls initSplines()
  virtual void initSplines(); //!< Allocates Splines

  //Scale factor stuff
  double a_0() const {return 1.0;}  //!< scale factor today
  
  //Hubble constant stuff
  double h() const { return h_h; } //!< Small Hubble h today
  double h2() const { return h_h * h_h; } //!< h^2 today

  virtual double E(double z) const; //!< H/H0 at given redshift

  double ZMax() const { return zmax; } //!< Maximum z to integrate out to

  //Omegas
  double omega_m() { return Omega_m; } //!< Omega matter
  double omega_v() { return Omega_v; } //!< Omega Vacuum
  double omega_relativistic() { return Omega_rel; } //!< Omega relativistic
  virtual double omega_0() { return omega_m()+omega_v()+omega_relativistic(); } //!< Return total omega today

  //Current Densities
  double rho_0() const { return  M_p()*M_p()*h2()* 3.3379e-7;  } ; //!< today's critical energy density in Mpc^-4
  double rho_m() { return rho_0() * omega_m(); } //!< today's matter energy density in Mpc^-4
  double rho_relativistic() { return rho_0() * omega_relativistic(); } //!< today's relativistic particle density in Mpc^-4
  double rho_v() { return rho_0() * omega_v(); } //!< today's vacuum energy density in Mpc^-4

  //Initial Densities -- or, really, density at a
  double initialRho_m(const double a) { return rho_m()*pow(a/a_0(),-3); } //! Density at scale factor a of matter
  double initialRho_relativistic(const double a) { return rho_relativistic()*pow(a/a_0(),-4); } //! Density at scale factor a of relativistic stuff
  double initialRho_v(const double a) { return rho_v(); } //! Density at scale factor a of vacuum energy

  double Z2iH(double z) { return Z2iHSpline->fastY(z); }  //!< Translates from redshift to inverse Hubble parameter

  //Distances
  double propermotionDistance(double z) const;//!< Returns proper motion distance
  double angulardiameterDistance(double z) const;//!< Returns angular diameter distance
  double luminosityDistance(double z) const;//!< Returns luminosity distance, assuming that zcmb=zhel
  double luminosityDistance(double zhel, double zcmb) const;//!< Returns luminosity distance at heliocentric redshift zhel and CMB frame redshift zcmb
  
  //Derivatives with distance
  double propermotionDistanceDeriv(double z) const;//!< Returns derivative of proper motion distance with respect to z
  double angulardiameterDistanceDeriv(double z) const;//!< Returns derivative of proper motion distance with respect to z
  double luminosityDistanceDeriv(double z) const;//!< Returns derivative of proper motion distance with respect to z

  //Set various things
  void seth(double h) { h_h = h; setOmegaH2_rel( 2.471e-5 ); }

  void setOmega_m(double x) { Omega_m = x; Omega_k=1.0-omega_0();} //!< Set Omega matter
  void setOmega_v(double x) { Omega_v = x; Omega_k=1.0-omega_0();} //!< Set Omega vacuum
  void setOmega_relativistic(double x) { Omega_rel = x; Omega_k=1.0-omega_0(); } //!< Set Omega relativistic
  void setOmegaH2_rel(double x) { Omega_rel = x / h2(); Omega_k=1.0-omega_0(); } //!< Set Omega*h^2 of relativistic matter

  void setZMax(double x) { zmax = x; if (validHistory()) reset(); } //!< Set maximum z 

  virtual void history(bool inform=false); //!< Evolve Universe, filling in the luminosity distance

};

#endif
