#include <cmath>
#include <iostream>
#include <limits>

#include "distcosmos.h"
#include "controlpanel.h"

using namespace std;

const double distCosmos::flattol = 1e-3; //!< Tolerance for flatness.  If abs(Omega_k) is less than this, the Universe is treated as flat
const int distCosmos::ndist = 500; //!< Number of luminostiy distance steps

distCosmos::distCosmos() : baseCosmos() {
  reset();

  zmax = 2;

  seth(0.65);  //Not very sensitive to this value
  setOmega_m(0.27);
  setOmegaH2_rel( 2.471e-5 );

  //Flat by default
  setOmega_v( 1 - omega_m() - omega_relativistic() );

  //Runge-Kutta by default
  odetype = Miscmath::rkutta;

}


void distCosmos::reset() {
  ValidHistory = false; 
  
  isRebounding = false;
  
  SplineAnchor.kill(); //! Get rid of all splines that are known to splineAnchor
  initSplines();
}


/*
  Init spline that distCosmos needs, which is only Z2iH.
*/
void distCosmos::initSplines() {

  distSpline = new Spline(ndist,"distSpline",&SplineAnchor);
  Z2iHSpline = new Spline(ndist,"z2ih",&SplineAnchor);

}

void distCosmos::history(bool inform) {

  double ciH = 2997.92458 / h_h;  //c/H0 in Mpc units

  if (validHistory()) throw Bad_Error("Cosmos::history() I already have a history! Call reset() first");
  
  double y[2];  //There is only 1 parameter, but the ridiculous NRecipes requires this

  y[1] = 0.0;

  double zstep = zmax / double(ndist - 1);


  double z = 0;
  double hnext = 0.1*zstep;  //Integrator initial guess step

  //Figure out if we are in a rebounding universe.
  //Technically, there is a maximum observable redshift in such
  // cases, so if zmax were less than this limit we would be fine,
  // but this calculation is actually fraught with complications
  // (the standard formulae you will find quoted various places
  //  don't actually work in all cases) so we avoid the situation
  //  by running away.
  isRebounding = calc_isrebounding();
  if (isRebounding) { 
    ValidHistory = false;
    return;
  }

  distSpline->set(0.0,0.0);  //z=0, dist=0
  Z2iHSpline->set(0.0,ciH);  //z=0, Z2iH = 1/H0

  double abs_ok = fabs(Omega_k);
  if (abs_ok < flattol) {
    //Flat Universe case  
    while (z < zmax) {
      hnext = Miscmath::odeint(y,1,z,z+zstep,1e-9,hnext,0,
			       (moDerivs)&distCosmos::iEode,*this, true,
			       odetype);
      z += zstep;
      distSpline->set( z, ciH*y[1] );
      Z2iHSpline->set( z, ciH/E(z) );
    }
  } else if (Omega_k > 0) {
    //Open Universe case
    double sq_abs_ok = sqrt( abs_ok );
    while (z < zmax) {
      hnext = Miscmath::odeint(y,1,z,z+zstep,1e-9,hnext,0,
			       (moDerivs)&distCosmos::iEode,*this, true,
			       odetype);
      z += zstep;
      distSpline->set( z, ciH/sq_abs_ok*sinh(sq_abs_ok*y[1]) );
      Z2iHSpline->set( z, ciH/E(z) );
    }
  } else {
    //Closed Universe Case
    double sq_abs_ok = sqrt( abs_ok );
    while (z < zmax) {
      hnext = Miscmath::odeint(y,1,z,z+zstep,1e-9,hnext,0,
			       (moDerivs)&distCosmos::iEode,*this, true,
			       odetype);
      z += zstep;
      distSpline->set( z, ciH/sq_abs_ok*sin(sq_abs_ok*y[1]) );
      Z2iHSpline->set( z, ciH/E(z) );
    }
  }

  //Last check on valid distances
  if ( isnan( y[1] ) ) {
    ValidHistory = false;
    return;
  }

  distSpline->arm(Spline::thoseReady); //Frees overhead memory
  Z2iHSpline->arm(Spline::thoseReady);

  //Fill in derivative spline
  dDdzSpline = new Spline(distSpline,"dDdzSpline",&SplineAnchor);
  distSpline->derive(*dDdzSpline);
  dDdzSpline->arm();

  ValidHistory = true;

  //Dump to file?
  if (inform) {
    string base = ControlPanel::cmbeasyDir("/output/");
    distSpline->dump(base+"dist");
    dDdzSpline->dump(base+"ddistdz");
    Z2iHSpline->dump(base+"z2iH");
  }

}

/*!
  This is a good function to overload when considering other
  cosmologies (Quintessence, w != -1, etc)
  \param z redshift
  \returns The value E(z), which is the Hubble parameter divided
    by it's current value -- i.e., H/H0.
*/
double distCosmos::E(double z) const {
  double opz = 1.0 + z;
  return sqrt( Omega_v + 
	       opz*opz*( Omega_k + opz*( Omega_m + opz*Omega_rel ) ) );
}

/*!
  \param z Redshift
  \param y Not used.  Needed for integrator compatability
  \param dy The value E(z) is returned in dy[1]
*/
void distCosmos::iEode(double z, const double* y,
		       double *dy) const {
  dy[1] = 1.0 / E(z);
}

//Figure out if Universe is rebounding from previous bout of contraction
bool distCosmos::calc_isrebounding() {
  //Formula from Carrol and Press, ARAA 1993, 30: 499-542
  //Ignoring Omega_relativistic
  double om = omega_m();
  double arg = (1.0 - om)/om;
  if ( om >= 0.5 ) {
    return omega_v() >= 4*om*pow( cos( acos( arg ) / 3.0 ), 3);
  } else {
    // cmath doesn't include inverse hyperbolic functions, so this
    // looks funny, but cosh^{-1}(z) = ln( z + sqrt(z+1)*sqrt(z-1) )
    return omega_v() >= 4*om*pow( cosh( log( arg + sqrt(arg*arg-1) )/ 3.0),3);
  }
}

double distCosmos::propermotionDistance(double z) const {
  if (! validHistory() )  throw Bad_Error("distCosmos::propermotionDistance() The history has not yet been cacluated -- call history");
  if ( z > zmax ) throw Bad_Error("distCosmos::propermotionDistance() Requested redshift outside of calculated range.");
  return distSpline->fastY(z);
} 

double distCosmos::propermotionDistanceDeriv(double z) const {
  if (! validHistory() )  throw Bad_Error("distCosmos::propermotionDistanceDeriv() The history has not yet been cacluated -- call history");
  if ( z > zmax ) throw Bad_Error("distCosmos::propermotionDistanceDeriv() Requested redshift outside of calculated range.");
  return dDdzSpline->fastY(z);
} 

double distCosmos::angulardiameterDistance(double z) const {
  if (! validHistory() )  throw Bad_Error("distCosmos::angulardiameterDistance() The history has not yet been cacluated -- call history");
  if ( z > zmax ) throw Bad_Error("distCosmos::angulardiameterDistance() Requested redshift outside of calculated range.");
  return distSpline->fastY(z)/(1.0+z);
}

double distCosmos::angulardiameterDistanceDeriv(double z) const {
  if (! validHistory() )  throw Bad_Error("distCosmos::angularDistanceDeriv() The history has not yet been cacluated -- call history");
  if ( z > zmax ) throw Bad_Error("distCosmos::angulardiameterDistanceDeriv() Requested redshift outside of calculated range.");
  double opz = 1.0+z;
  return dDdzSpline->fastY(z)/opz - distSpline->fastY(z)/(opz*opz);
}

double distCosmos::luminosityDistance(double z) const { 
  if (! validHistory() )  throw Bad_Error("distCosmos::luminosityDistance() The history has not yet been cacluated -- call history");
  if ( z > zmax ) throw Bad_Error("distCosmos::luminosityDistance() Requested redshift outside of calculated range.");
  return (1+z)*distSpline->fastY(z);
}

double distCosmos::luminosityDistance(double zhel, double zcmb) const { 
  if (! validHistory() )  throw Bad_Error("distCosmos::luminosityDistance() The history has not yet been cacluated -- call history");
  if ( zcmb > zmax ) throw Bad_Error("distCosmos::luminosityDistance() Requested redshift outside of calculated range.");
  return (1+zhel)*distSpline->fastY(zcmb);
}

double distCosmos::luminosityDistanceDeriv(double z) const { 
  if (! validHistory() )  throw Bad_Error("distCosmos::luminosityDistance() The history has not yet been cacluated -- call history");
  if ( z > zmax ) throw Bad_Error("distCosmos::luminosityDistanceDeriv() Requested redshift outside of calculated range.");
  return (1+z)*dDdzSpline->fastY(z) + distSpline->fastY(z);
}

