#include "coupledexp.h"
#include "coupledquintcosmos.h"

CoupledExp::CoupledExp(CoupledQuintCosmos& c) : Exponential(c) , imitateConstant(false), B(0.1)
{
  mQinitial = -1;
  mHasQInitial = false;
  mQinitialFactor = 1.;
  mHasQInitialFactor = false;
}

//! Critical density on the radiation attractor
static double rhoCritRDE(const double rho_relativistic, const double lambda, const double beta)
{
  return rho_relativistic/(1.-(4.-beta*(lambda-beta))/(lambda-beta)/(lambda-beta));
}

double CoupledExp::initialQ(const double a) const {
  const double beta = B; // the coupling constant
  double phi_initial;

  if (mHasQInitial) {
    phi_initial = mQinitial;
    //cout << "have initial Q set so I am returning phi_initial= " << phi_initial
    //<< "(i.e. phi_initial/M_p=" << (phi_initial/M_p) << ")" << endl;
    return mQinitial;
  }

//X   phi_initial= Exponential::initialQ(a);
//X   return phi_initial;
//X   cout << "would be returning phi_initial= " << phi_initial << "(i.e. phi_initial/M_p=" << (phi_initial/M_p) << ")" << endl;

  if (cosmos.rho_nuNR0() != 0)
    throw Bad_Error("CoupledExp::initialQ() - take massive nu into account");
  const double rhoRDECrit = rhoCritRDE(cosmos.initialRho_nu(a)+cosmos.initialRho_g(a), lambda(), beta);
  const double Omega_potRDE = (4.-3.*beta*(lambda()-beta))/(3.*pow(lambda()-beta, 2));
  phi_initial = -M_p/lambda()*log(fabs(Omega_potRDE*rhoRDECrit)/cosmos.M_p(4));
  if (mHasQInitialFactor)
    phi_initial *= mQinitialFactor;
//X   cout << "am returning: " << phi_initial << "(i.e. phi_initial/M_p=" << (phi_initial/M_p) << ")" << endl;
  return phi_initial;
}

double CoupledExp::initialQDot(const double a)  const
{
  const double beta = B; // the coupling
  const double qDot = 4.*a*sqrt(fabs(V(initialQ(a),a)/(4.-3.*beta*(lambda()-beta))));
  return qDot;
}

double CoupledExp::lambda() const { 
  const double beta = B; // the coupling
  const double om_q = cosmos.omega_q(false);
  return beta*(1.-1./(2.*om_q))+sqrt(beta*beta/(4.*om_q*om_q)+3./om_q);
}

void CoupledExp::printStatus() {
  cout << "CoupledExponential: Coupling is beta = B = " << B <<  endl;
  Exponential::printStatus();
}
