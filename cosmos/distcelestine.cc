
#include <cmath>
#include <limits>
#include "distcelestine.h"

using namespace std;

distCelestine::distCelestine() : distCosmos() {

  setOmega_v( 0.0 ); //No cosmological constant
  setOmega_q( 1 - omega_m() - omega_v() - 
	      omega_relativistic() ); //Flat by default

  setW0( -0.8 );
  setWa( 0.0 );

}

/*!
  This function is supposed to determine if the universe is rebounding
  from a previous bout of contraction.  Unfortunately, we don't know
  how to do this in this particular case, so we just assume it isn't --
  unless we have set Omega_quintessence to zero, in which case we know
  what to do.
*/
bool distCelestine::calc_isrebounding() {
  if ( Omega_q == 0 ) return distCosmos::calc_isrebounding();
  if ( w_a == 0 && w_0 == -1 ) {
    //Cosmological constant case
    //Formula from Carrol and Press, ARAA 1993, 30: 499-542
    //Ignoring Omega_relativistic
    double om = omega_m();
    double arg = (1.0 - om)/om;
    if ( om >= 0.5 ) {
      return (omega_v() + omega_q()) >= 
	4*om*pow( cos( acos( arg ) / 3.0 ), 3);
    } else {
      // cmath doesn't include inverse hyperbolic functions, so this
      // looks funny, but cosh^{-1}(z) = ln( z + sqrt(z+1)*sqrt(z-1) )
      return (omega_v() + omega_q()) >= 
	4*om*pow( cosh( log( arg + sqrt(arg*arg-1) )/ 3.0),3);
    }
  }    
  return false;
}

/*!
  This is equal to exp[ 3 * int_0^z (1 + w(z))/(1+z) ] specialized
  to the current case.
  \param z Redshift
*/
double distCelestine::quintEfac(double z) const {
  double opz = 1.0 + z;
  if ( w_a == 0 ) {
    //Simple case
    return pow( opz, 3*(1+w_0) );
  } else {
    return exp( - 3*w_a*z/opz ) * pow( opz, 3 * ( 1 + w_0 + w_a ) );
  }
}

/*!
  \param z Redshift
*/
double distCelestine::E(double z) const {
  double opz = 1.0 + z;
  return sqrt( Omega_v + quintEfac(z) * Omega_q + 
	       opz *opz * ( Omega_k + opz * ( Omega_m + opz*Omega_rel ) ) );
}
