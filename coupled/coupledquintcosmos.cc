#include "coupledquintcosmos.h"
#include "quintcosmos.h"
#include "cosmos.h"
#include "global.h"
#include "exponential.h"
#include "coupledexp.h"
#include "coupledleaping.h"
#include "coupling.h"
#include "exponentialcoupling.h"
#include "arthur.h"
#include "ratra.h"
#include "arbitrary.h"
#include "celestine.h"
#include "voidQuintessence.h"
#include "cleanvector.h"

#include <fstream>
#include <cmath>
#include <map>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "crossover.h"
#include "crossoverfield.h"

#include "splinetools.h"
#include "massiveneutrinos.h"

#include <fstream>
#include <map>
#include <math.h>
#include <cmath>
#include <limits>

CoupledQuintCosmos::CoupledQuintCosmos(Quintessence::Type quinttype)
              : QuintCosmos(quinttype)
{
}

CoupledQuintCosmos::~CoupledQuintCosmos() {
//  delete quint;
}

void CoupledQuintCosmos::normExtraSplines(){
cout << "CoupledQuintCosmos::normExtraSplines()" << endl;
	double z, a, t, tau;
	double rho_tot, rho_cdm, hub;
	double m_factor, hub_factor;

  A2Hubble_factor = new Spline(splineA2Tau,"a2hubble_factor",&SplineAnchor);
  A2Mass_factor = new Spline(splineA2Tau,"a2mass_factor",&SplineAnchor);
  Z2Hubble_factor = new Spline(Z2iH_spline,"z2hubble_factor",&SplineAnchor);

	int TOT=Tau2A->size();
	double H_0 = tau2Hubble(tau_0())*M_p();
	double rhocrit = 3*H_0*H_0;
	double rho_cdm0 = Tau2Rho_m->fastY(tau_0())-Tau2Rho_b->fastY(tau_0());

	for(int k=0; k<TOT; k++) {
	tau = Tau2A->start(k);
	a = Tau2A->front(k);
	z = 1./a - 1.;

	rho_cdm = Tau2Rho_m->fastY(tau)-Tau2Rho_b->fastY(tau);
  	rho_tot = Tau2Rho->fastY(tau);

  	hub_factor = rho_tot/rhocrit;
	m_factor = rho_cdm/(rho_cdm0*pow(a,-3));

	A2Hubble_factor->set(hub_factor); 
	Z2Hubble_factor->set(hub_factor); 
	A2Mass_factor->set(m_factor); 
}
	
	A2Hubble_factor->arm();
	A2Hubble_factor->dump("hubble_factor");
	A2Mass_factor->arm();
	A2Mass_factor->dump("mass_factor");
}


void CoupledQuintCosmos::propagateHistoryInTau(const double tau,const double* y,double* dy) {
  if (couple==0) cout << "Error in CoupledQuintCosmos::propagateHistoryInTau: coupling not set!" << endl;
  double a = y[1];
//if(a < a_0()*1.5){
//cout << "PropagateHistoryInTau()\n" ;
  double rho_g = y[2];
  double rho_b = y[3];
  double rho_c = y[4];
  double rho_nu = y[6];
  double q = y[8];
  double qprime = y[9];
  quint->setQ(q);
  quint->setQDot(qprime);  // clash of conventions: here both Dot and prime mean "dphi/dtau"
// for debugging, so that overwriting phi in the quint class is enough
  q = quint->q(a);
  qprime = quint->qDot(a);

  const double rho_q = quint->rho(a);

//if(a<1.1e-12) cout << "a: " << a  << " a_0: " << a_0() << " rho_c: " << rho_c << " rho_q:  " << rho_q << endl; 
//for(int k=0; k<N_VAR; k++) cout << "y[" << k << "]: " << y[k] << endl;

  double rho_nuNR = 0;
  double p_nuNR = 0.;
  if (amnu != 0.0) {
    double rnu, pnu;
    MassiveNeutrinos::nu1(a*amnu, &rnu, &pnu);
    rho_nuNR = rho_nu0() * nuNR() * rnu/(a*a*a*a);
    p_nuNR = rho_nu0() * nuNR() * pnu/(a*a*a*a);
  }


  const double totalRho = rho_g + rho_c + rho_b + rho_nu  + rho_nuNR + rho__v + rho_q;
  const double H = a* sqrt(1.0/3.0*totalRho)/M_p();

  if (isnan(H)) 
{
    cout << "cqc.cc: PropagateHistoryInTau():" << endl;
    cout << "isnan H: tau is " << tau << "  a: " << a << endl;
    cout << "quintrho: " << quint->rho(a) << endl;
    cout << "quint->rho_kin(a): " << quint->rho_kin(a) << " " << quint->q(a) <<  " " << quint->V(quint->q(a),a) << endl;
    cout << "qdot/Mp: " <<  (quint->qDot(a)/M_p()) << endl;
    throw Bad_Error("H is nan.");
  }

  dy[1] = H*a;      // da/dtau is H*a
  dy[2] = -4.0*H*rho_g;  // gammas
  dy[3] = -3.0*H*rho_b;  // baryons
  dy[4] = -3.0*H*rho_c - couple->phi2Q_c_0(q, qprime, rho_c);  // CDM - modified for coupling
  dy[5] = a;   // dt/dtau
  dy[6] = -4.0*H*rho_nu;  // massless neutrinos

  if (quint->needsPropagation()) {
    dy[8] = qprime;  // quintessence
    dy[9] = -2.*H*qprime- a*a*quint->Vprime(q,a,H) - a*a*couple->phi2S(q, rho_c);
  } else {
    dy[8] = 0;
    dy[9] = 0;
  }

  double r = 3.0/4.0 * rho_b/rho_g;
  double cs = 1.0/sqrt(3*(1+r));  // sound speed
  dy[7] = cs;  // integrate for sound - horizon
//} // if a < a_0
}


void CoupledQuintCosmos::history(bool inform) {
  if (validHistory()) {
    throw Bad_Error("CoupledQuintCosmos::history() I already have a history! Call reset() first");
  }
  getReady();
  quint->prepare();
  rho__v = omega2rho(omega_v());

  cross_W25 = -4.6;

  double y[10];
  double ydot[10];
  initializeDensities(y, ydot);


  if (quint != 0 && inform)
    quint->printStatus();
  // now we want to get da/dtau in order to estimate the beginning tau....
  propagateHistoryInTau(0,y,ydot);

  double tau = y[1]/ydot[1];

  if (inform) cout << "from initial step, we have: " << tau << endl;

  double dtau = tau*1e-1;     // and dtau is one tenth of this..

  // zeroing tau__equ and tau__ls as they will now be (re)calculated
  tau__equ = 0;
  tau__ls = 0;

  //  ofstream monitor("monit.dat");

  double A=0;
  Background b;
  double hnext = dtau*1e-12; // 1e-3;
  while (y[1] < a_0()) {
    hnext = Miscmath::odeint(y,9,tau,tau+dtau,1e-9,hnext,0,(moDerivs)&CoupledQuintCosmos::propagateHistoryInTauWrapper, *this);
    tau += dtau;

    propagateHistoryInTau(tau,y,ydot);     // We would like to have the tau derivatives for splines
    b = fillHistorySplines(tau,y,ydot);         // this also updates quintessence to the right q and qdot

    A = b.A;

    for (int i=1; i<=7; i++) {
      if (isnan(y[i])) {
        cout << endl << "quintcosmos:: y["<<i<<"] isnan at time: " << tau << "  [ " << A << " ]"<<endl;
        throw Bad_Error("quintcosmos isnan !");
      }
    }

    // determine tau__equ by comparing densities
    // this is rough and should is refined by really narrowing in on
    // the zeros of the equality later
    if (b.rho_b + b.rho_c > b.rho_g && tau__equ == 0) tau__equ = tau;

    // determine t_LS by noting the crossing of a above the a_ls level...
    // rough, as above,  is refined in a later step by looking for
    // the interpolation crossing of a(tau) >= a_ls()
    if (A>= a_ls() && tau__ls == 0) tau__ls = tau;

    if (A < 0.1*a_0()) dtau *= 1.05;
    dtau = min(10.0,dtau);

    if (!mExtraTimeSteps.empty())
      checkForExtraTimestep(tau, dtau);

  }

  finalizeHistorySplines(inform);

  findSpecialMoments();
	ageYr = mpc2year()*tau2t(tau_0());

  if (inform)
    printStatus();
}

void CoupledQuintCosmos::initializeDensities(double* y, double* ydot)
{
  for (int k = 0; k < 9; k++) { y[k]=0; ydot[k] = 0;}
  y[1] = initialScaleFactor(); // we start from here (almost, cause we will note only from second step on)
  y[2] = initialRho_g(y[1]);
  y[6] = initialRho_nu(y[1]);
  y[3] = initialRho_b(y[1]);
  y[4] = initialRho_cdm(y[1]);
  y[5] = 0;
  y[7] = 0;
  quint->setInitial(y[1]);  // Inform quintessence to set its conditions
  y[8] = quint->q(y[1]);
  y[9] = quint->qDot(y[1]);

  //  cout << "setting initial q to " << (y[8]/M_p()) << endl;
  //  cout << "setting initial q' to " << (y[9]/M_p()) <<  endl;

  double atmp = y[1];
  double rhonuNR = 0.;
  double pnuNR = 0.;
  if (haveMassiveNeutrinos()) {
    MassiveNeutrinos::initnu1(amnu);
    double rnu, pnu;
    const double massNukbT = mass_nu_kbT();
    MassiveNeutrinos::nu1(atmp*massNukbT, &rnu, &pnu);
    rhonuNR = rho_nu0() * nuNR() * rnu/(atmp*atmp*atmp*atmp);
    pnuNR = rho_nu0() * nuNR() * pnu/(atmp*atmp*atmp*atmp);
  }
}
