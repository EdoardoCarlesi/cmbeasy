#include "exponential.h"
#include "quintcosmos.h"

#include <cmath>
#include <iostream> // for debugging

Exponential::Exponential(QuintCosmos &c)
   : Quintessence(c) , imitateConstant(false) ,  mV0(cosmos.M_p(4)),
     mHasInitialQ(false), mHasInitialQDot(false), mInitialQ(-1e100),
     mInitialQDot(-1e100), mHasLambda(false), mLambda(-1e100)
{
  reset();
  param.resize(4);
}

void Exponential::setParameters(const vector<double> &param)
{
  Quintessence::setParameters(param);

  if (!isinf(param[0])) {
    setLambda(param[0]);
  } else {
    cout << "No lambda set, using default parameters!" << endl;
  }

  if (!isinf(param[1])) {
    setInitialQ(param[1]);
  } else {
    cout << "No initial Q set, using default parameters!" << endl;
  }

  if (!isinf(param[2])) {
    setInitialQ(param[2]);
  } else {
    cout << "No initial Q_Dot set, using default parameters!" << endl;
  }

  if (!isinf(param[3])) {
    setV0(param[3]);
  } else {
    cout << "No V0 set, using default parameters!" << endl;
  }
}

/*!
  V0*exp(-lambda * phi/M_p() ) potential
  scale factor is not used, V0 is Mp^4 by default
*/
double Exponential::V(const double q, const double a) const
{
	//cout << " v0  " << V0() << " l:  " << lambda() << " q: " << q << " Mp: " << cosmos.M_p() << endl;
  return V0()*exp(-lambda()*q/cosmos.M_p());
}

double Exponential::Vprime(const double q, const double a,const double adotoa) const
{
  return -lambda()*V0()/cosmos.M_p()*exp(-lambda()*q/cosmos.M_p());
}

double Exponential::Vprime2(const double q, const double a,const double adotoa,const double tau) const
{
  return lambda()*lambda()*V0()/cosmos.M_p(2)*exp(-lambda()*q/cosmos.M_p());
}

//! set the initial field value
void Exponential::setInitialQ(const double q)
{
  if(isinf(q) || isnan(q)) {
    throw Bad_Error("Exponential::setInitialQ() - q is inf or nan.");
  }

  mHasInitialQ=true;
  mInitialQ=q*cosmos.M_p();
}

void Exponential::setInitialQDot(const double qDot)
{
  mHasInitialQDot=true;
  mInitialQDot=qDot*cosmos.M_p();
}

void Exponential::setLambda (const double lambda)
{
  mHasLambda=true;
  mLambda=lambda; /*cout << " ---- lambda: " << mLambda << endl;*/
}

void Exponential::setV0 (const double vzero)
{
  //mV0=vzero*cosmos.M_p(2); 
  mV0=vzero*cosmos.M_p(4); 
  //cout << " ---- vzero: " << mV0  << endl;
}

/*! returns either the field value set by setInitialQ(), or
 *  puts it  on the radiation attractor
 */
double Exponential::initialQ(const double a) const
{
  if(mHasInitialQ) {
    return mInitialQ;
  }

  if (imitateConstant) {
    return  -cosmos.M_p()/lambda() * log(cosmos.omega2rho(cosmos.omega_v()) / cosmos.M_p(4));
  }

  return initialQFromRadiationAttractor(a);
}

double Exponential::initialQFromRadiationAttractor(const double a) const
{
  double rhocrit = (cosmos.initialRho_g(a) + cosmos.initialRho_nu(a) + cosmos.initialRho_nuNR(a) ) / (1. - 4.0/(lambda()*lambda()));
  double pot = 4.0/(3.0*lambda()*lambda())*rhocrit;
  double phiInitial = -cosmos.M_p()/lambda() * log(pot/V0());
  cout << "Exponential::initialQFromRadiationAttractor returning: " << (phiInitial/cosmos.M_p()) << endl;
  return phiInitial;
}

/*! dphi/dtau initial as set by setInitialQDot,
 * or as function of a (actually, the initial photon density).
*/
double Exponential::initialQDot(const double a)  const
{
  if(mHasInitialQDot) {
    return mInitialQDot;
  }
  return sqrt(V(initialQ(a)))*2.0*a;
}

double Exponential::lambda() const
{
  if(mHasLambda) {
    return mLambda;
  }

  if (cosmos.omega_q(false) > 0) {
    return sqrt(3.0 / cosmos.omega_q(false));
  }
  return 8.0;
}

void Exponential::printStatus()
{
  cout << "-------------------" << endl;
  cout << "Exponential Quintessence " << endl;
  cout << "lambda: " << lambda() << endl;
  cout << "V0: " << ( V0()/cosmos.M_p(4) ) << " *M_p()^4" << endl;
  cout << "omega_q() " << cosmos.omega_q() << endl;
  Quintessence::printStatus();
}
