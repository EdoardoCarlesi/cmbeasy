#ifndef DISTCELESTINE_H
#define DISTCELESTINE_H

#include "distcosmos.h"

/*!
  distCosmos derived object to represent anything with dark energy equation 
  of state parameter 

  w(a) = w_0 + w_1*(1-a)

  which means that

  E(z) = sqrt( (1+z)^3 omega_m + (1+z)^2 omega_k + f(z) omega_q )
  
  with

  f(z) = exp( - 3 w_1 z / (1+z) ) * ( 1 + z )^(3*(1+w_0+w_1))

  Note that this class actually supports a combination of the
  quintessence field and a cosmological constant, if one is so
  inclined.

*/

class distCelestine : public distCosmos {
protected:
  
  double w_0; //!< w_0 in w(a) expression
  double w_a; //!< w_a in w(a) expression

  double Omega_q; //!< Omega_quintessence
  
  bool calc_isrebounding(); //!< Partial implementation of calculation of whether or not Universe is rebounding

  double quintEfac(double) const; //!< Factor in front of Omega_q in E(z)

public :
    
  distCelestine(); //!< Default constructor

  double E(double) const; //!< H/H0 at given redshift

  double omega_q() { return Omega_q; } //!< Omega Quintessence

  double omega_0() { return omega_m() + omega_v() + omega_relativistic() +
		       omega_q(); } //!< Return total omega today

  double w0() const { return w_0; } //!< Return w_0
  double wa() const { return w_a; } //!< Return w_a

  double rho_q() { return rho_0() * omega_q(); } //!< Density in quintessence field today in Mpc^-4

  double initialRho_q(double a) { return rho_q() * 
				    quintEfac( a_0() / a - 1 ); } //!< Density at given scale factor of quintessence field

  //Set things
  void setOmega_q(double x) { Omega_q = x; Omega_k = 1.0 - omega_0(); } //!< Set Omega quintessence
  void setW0(double x) { w_0 = x; } //!< Set w_0
  void setWa(double x) { w_a = x; } //!< Set w_a

};

#endif
