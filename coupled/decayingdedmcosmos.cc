#include "decayingdedmcosmos.h"
#include "quintcosmos.h"
#include "global.h"
#include "coupledexp.h"
#include "coupledleaping.h"
#include "coupling.h"
#include "exponentialcoupling.h"
#include "anchor.h"
#include "cleanvector.h"

// included only for printOSF
#include "cninvariant.h"

//backward debug
#include "splinetools.h"

#include <fstream>
#include <sstream> // debugging only
#include <iostream> // debugging only
#include <iomanip>
#include <map>
#include <cmath>
#include <limits>

struct BackgroundDEDM;

DecayingDEDMCosmos::DecayingDEDMCosmos(Quintessence::Type quinttype)
              : QuintCosmos(quinttype), bprimespline(0), bprimeprimespline(0)
{
  bspline = 0; //new Spline(1000,"bspline",&anchor);
}

DecayingDEDMCosmos::~DecayingDEDMCosmos() {
}
 
/*
 * Coupled models stuff - to be added e.g. effective eos, variable coupling...
 * */
void DecayingDEDMCosmos::initSplines(){
	
	QuintCosmos::initSplines();

A2Rho_cdm_1 = new Spline(splineA2Tau, "a2_rho_cdm_1", &SplineAnchor);
A2Rho_cdm_2 = new Spline(splineA2Tau, "a2_rho_cdm_2", &SplineAnchor);
A2Rho_cdm_tot = new Spline(splineA2Tau, "a2_rho_cdm_tot", &SplineAnchor);
A2Rho_q = new Spline(splineA2Tau, "a2_rho_q", &SplineAnchor);
A2Rho_b = new Spline(splineA2Tau, "a2_rho_b", &SplineAnchor);
A2Rho_g = new Spline(splineA2Tau, "a2_rho_g", &SplineAnchor);
}

void DecayingDEDMCosmos::reset(){
	QuintCosmos::reset();
}

void DecayingDEDMCosmos::seth(double h)
{
  QuintCosmos::seth(h);
}

void DecayingDEDMCosmos::propagateHistoryInTau(const double tau,const double* y,double* dy) {
  if (couple==0) cout << "Error in DecayingDEDMCosmos::propagateHistoryInTau: coupling not set!" << endl;

  const double a = y[1];

  const double rho_g = y[2];
  const double rho_b = y[3];
  const double rho_c1 = y[4];
  const double rho_c2 = y[10];
  const double rho_nu = y[6];
  
  //cout << "DM1: " << rho_c1 << " DM2: " << rho_c2 << endl; 

  double q = y[8];
  double qprime = y[9];
  quint->setQ(q);
  quint->setQDot(qprime);  // clash of conventions: here both Dot and prime mean "dphi/dtau"
			   // for debugging, so that overwriting phi in the quint class is enough
  q = quint->q(a);
  qprime = quint->qDot(a);
  const double rho_q = quint->rho(a);

  double rho_nuNR = 0;

  const double totalRho = rho_g + rho_c1 + rho_c2 + rho_b + rho_nu  + rho_nuNR + rho__v + rho_q;
  const double H = a* sqrt(1.0/3.0*totalRho)/M_p();
//cout << "rho_g, rho_c1, rho_c2, rho_q" << rho_g << rho_c1 << rho_c2 << rho_q << endl;
  if (isnan(H)) {
    cout << "cqc.cc: PropagateHistoryInTau():" << endl;
    cout << "isnan H: tau is " << tau << "  a: " << a << endl;
    cout << "quintrho: " << quint->rho(a) << endl;
    cout << "quint->rho_kin(a): " << quint->rho_kin(a) << " " 
	 << quint->q(a) <<  " " << quint->V(quint->q(a),a) << endl;
    cout << "qdot/Mp: " <<  (quint->qDot(a)/M_p()) << endl;
    throw Bad_Error("H is nan.");
  }
  
double fact = 1.;//0.8;

  dy[1] = H*a;      // da/dtau is H*a
  dy[2] = -4.0*H*rho_g;  // gammas
  dy[3] = -3.0*H*rho_b;  // baryons
  dy[4] = -3.0*H*rho_c1 - (1./tau_dm)*rho_c1; 
  dy[10]= -3.0*H*rho_c2 - couple->phi2Q_c_0(q, qprime, rho_c2) + fact*(1./tau_dm)*rho_c1;  

  dy[5] = a;   // dt/dtau
  dy[6] = -4.0*H*rho_nu;  // massless neutrinos

    dy[8] = qprime;  // quintessence
    dy[9] = -2.*H*qprime- a*a*quint->Vprime(q,a,H) + a*a*couple->phi2S(q, rho_c2);

  double r = 3.0/4.0 * rho_b/rho_g;
  double cs = 1.0/sqrt(3*(1+r));  // sound speed
  dy[7] = cs;  // integrate for sound - horizon
}

void DecayingDEDMCosmos::findSpecialMoments() {
	QuintCosmos::findSpecialMoments();	
}

void DecayingDEDMCosmos::history(bool inform) {
	cout << " history() called" << endl;

  if (validHistory()) throw Bad_Error("QuintCosmos::history() I already have a history! Call reset() first");

	getReady();
  quint->prepare();
  rho__v = omega2rho(omega_v());
 
  double y[11];
  double ydot[11];

  y[1] = 1e-12;           // we start from here (almost, cause we will note only from second step on)
  y[2] = initialRho_g(y[1]);
  y[3] = initialRho_b(y[1]);
  y[4] = initialRho_cdm(y[1]);
  //y[10] = 0.001*y[4]; // The second DM specie starts with null density
  y[10] = 0; //.001*y[4]; // The second DM specie starts with null density

  //cout << " initialRho_c: " << y[4]<< endl;

  y[6] = initialRho_nu(y[1]);
  y[5] = 0;
  y[7] = 0;

  quint->setInitial(y[1]);  // Inform quintessence to set its conditions

  double q = quint->q(y[1]);
  double totalRhotmp = y[2]+y[6]+y[3]+y[4]+y[10];
  double atmp = y[1];
  const double H = atmp* sqrt(1.0/3.0*totalRhotmp)/M_p();

  y[8] = quint->q(y[1]);
  y[9] = quint->qDot(y[1]);

  if (quint != 0 && inform) quint->printStatus();
  // now we want to get da/dtau in order to estimate the beginning tau....
  
	propagateHistoryInTau(0,y,ydot);

  double tau = y[1]/ydot[1];
  double lasttau = 0.;

  if (inform) cout << " from initial step, we have tau=: " << tau 
                   << " and adotoa=" << (ydot[1]/y[1])    << endl;

  double dtau = tau*1e-1;     // and dtau is one tenth of this..
  // zeroing tau__equ and tau__ls as they will now be (re)calculated
  tau__equ = 0;
  tau__ls = 0;

  double A=0;
  double totalRho,totalPressure;
  BackgroundDEDM b;
  double hnext = dtau*1e-12; // 1e-3;

  while (y[1] < a_0()) {
	  
hnext = Miscmath::odeint(y,10,tau,tau+dtau,1e-9,hnext,0,(moDerivs)&QuintCosmos::propagateHistoryInTauWrapper, *this);

    tau += dtau;

    bool writeSplines = false; 
	//cout << " pre-call propagate..." << endl;
    propagateHistoryInTau(tau,y,ydot);     // We would like to have the tau derivatives for splines
	//cout << " post-call propagate..." << endl;
    b = fillHistorySplinesDEDM(tau,y,ydot, writeSplines);
//cout << "fill history" << endl;

    if ((lasttau != 0) && ((tau-lasttau)<1e-22))    // make sure that we are not sampling 
	    					    // our splines too densly (otherwise
    			                     // we get overflows in the spline interpolation)
    if ( !(y[1] < a_0()) || ((A/y[1]) < 0.99)) writeSplines = false; // force to write splines (otherwise we
                                            // might end up with the last a being 0.2 or something
                                            //we'll probably get spline interpolation errors,
                                            //though...

//cout << " b.rho_c1 : " << b.rho_c1 << " b.rho_c2 : " << b.rho_c2 << endl; 
    for (int i=1; i<=7; i++) {
      if (isnan(y[i])) {
        cout << endl << "quintcosmos:: y["<<i<<"] isnan at time: " << tau << "  [ " << A << " ]"<<endl;
        throw Bad_Error("quintcosmos isnan !");
      }
    }

    // determine tau__equ by comparing densities
    // this is rough and should is refined by really narrowing in on
    // the zeros of the equality later
    if (b.rho_b + b.rho_c1 + b.rho_c2 > b.rho_g && tau__equ == 0) {
      tau__equ = tau;
//X       cout << "Equality at roughly: tau= " << tau << ", a=" << A << endl;
    } else if (tau__equ == 0) {
//X       cout << "No eq. yet; " << endl;
    }

    // determine t_LS by noting the crossing of a above the a_ls level...
    // rough, as above,  is refined in a later step by looking for
    // the interpolation crossing of a(tau) >= a_ls()
    if (A>= a_ls() && tau__ls == 0) tau__ls = tau;

    if (A < 0.1*a_0()) dtau *= 1.05;
    dtau = min(10.0,dtau);
  }
  if (inform)  printStatus();

}

BackgroundDEDM  DecayingDEDMCosmos::fillHistorySplinesDEDM(double tau, const double *y, const double *ydot, bool writeSplines) {
  BackgroundDEDM b;
  double A  = b.A = y[1]; // added normalization factor
  b.rho_g   = y[2];
  b.rho_nu  = y[6];
  b.rho_b   = y[3];
  b.rho_c1  = y[4];
  b.rho_c2  = y[10];
  b.rho_c   = y[4] + y[10];
  b.t = y[5];
  b.rs = y[7];
  b.quint = y[9];

  double plain =0;
  double rnu=0, pnu=0;

  double totalRho = b.rho_g + b.rho_nu + b.rho_nuNR + b.rho_c + b.rho_b  + rho__v + quint->rho(A);
  double totalPressure = 1.0/3.0*(b.rho_g + b.rho_nu) + b.p_nuNR - rho__v + quint->p(A);
  if (amnu !=0.0) {
    MassiveNeutrinos::nu1(b.A*amnu,&rnu,&pnu);
    plain = rho_nu0() * nuNR() /(A*A*A*A);
  }

  b.rho_nuNR = rnu * plain;
  b.p_nuNR = pnu*plain;

    double rho_m = b.rho_b + b.rho_c;

    // scale factor and its derivative
    Tau2A->set(tau,b.A);
    Tau2ADot->set(ydot[1]);
    Tau2ADotToA->set(ydot[1]/y[1]);

    // the photons and massless and massive (NonRelativistic) neutrinos
    Tau2Rho_g->set(b.rho_g);
    Tau2Rho_nu->set(b.rho_nu);
    Tau2Rho_nuNR->set(b.rho_nuNR);
    Tau2P_nuNR->set(b.p_nuNR);
    Tau2w_nuNR->set(b.p_nuNR/b.rho_nuNR);

    // baryons, cdm and total matter
    Tau2Rho_b->set(b.rho_b);
    Tau2Rho_m->set(rho_m);
    Tau2Rho_cdm->set(b.rho_c);

    // 4/3 of the ratio of gammas over baryons (important for soundspeed)
    Tau2R->set(4.0/3.0*b.rho_g/b.rho_b);

    // getting from a <-> tau , t <-> tau
    splineA2Tau->set(A,tau);
    splineT2Tau->set(b.t,tau);

    // t(tau), soundhorizon(tau) and soundspeed^2(tau)
    Tau2T->set(b.t);
    Tau2Rs->set(b.rs);
    Tau2Cs2->set(ydot[7]*ydot[7]);
    Tau2Cs->set(ydot[7]);

    Z2iH_spline->set(a2z(A), A*A / ydot[1]);
    //Z2H_spline->set(a2z(A), ydot[1]/(A*A));

 	Tau2Rho_qpot->set(quint->V(quint->q(A), A));
 	Tau2Rho_qkin->set(quint->rho_kin(A));
  	Tau2Rho_q->set(quint->rho(A));
 	Tau2Phi->set(quint->q(A));
  	Tau2PhiDot->set(quint->qDot(A));

      	Tau2W_q->set(quint->w(A));
  Tau2Omega_q->set(quint->rho(A) / totalRho);
  Tau2Omega_qw ->set(quint->rho(A) / totalRho*quint->w(A));

  Tau2Vprime->set(quint->Vprime(quint->q(A),A,ydot[1]/A));
  Tau2VprimeOverRho->set(quint->Vprime(quint->q(A),A,ydot[1]/A)/quint->rho(A));

  // set eqn of state etc as function of logA (for Jan), also note crossing of
  // w_q = -0.25 if a > 0.01...
  LogA2Omega_q->set(log(A),quint->rho(A) / totalRho);
  LogA2Omega_qw ->set(quint->rho(A) / totalRho / quint->w(A));

  Tau2Rho->set(totalRho);
  Tau2P->set(totalPressure);
  Tau2P_q->set(quint->p(A));

  // 8*pi*G*a^2*rho
  Tau2GRho->set(Gpi8()*A*A*totalRho);

  A2Omega_q->set(quint->rho(A) / totalRho);
  A2Omega_qw ->set(quint->rho(A) / totalRho*quint->w(A));

    // Different cdm species
    A2Rho_cdm_1->set(b.rho_c1);
    A2Rho_cdm_2->set(b.rho_c2);
    A2Rho_cdm_tot->set(b.rho_c);
    A2Rho_q->set(quint->rho(A));
    A2Rho_b->set(b.rho_b);
    A2Rho_g->set(b.rho_g);

  Z2W->set(quint->w(A));

  return b;
}

double DecayingDEDMCosmos::expDecay(double tau, double a) {
double fac;
fac = 1.;
//fac = exp(tau*a / tau_dm);
return fac;

}

void DecayingDEDMCosmos::dumpAll(){

A2Rho_cdm_1->arm();
A2Rho_cdm_2->arm();
A2Rho_cdm_tot->arm(); 
A2Rho_q->arm();
A2Rho_b->arm();
A2Rho_g->arm(); 

A2Rho_cdm_1->dump();
A2Rho_cdm_2->dump();
A2Rho_cdm_tot->dump(); 
A2Rho_q->dump();
A2Rho_b->dump();
A2Rho_g->dump();

double rhoTot = A2Rho_cdm_tot->fastY(1.) +  A2Rho_q->fastY(1.) + A2Rho_b->fastY(1.);
double o1, o2, ob, oq;

o1 = A2Rho_cdm_1->fastY(1.) /rhoTot;
o2 = A2Rho_cdm_2->fastY(1.) /rhoTot;
ob = A2Rho_b->fastY(1.) /rhoTot;
oq = A2Rho_q->fastY(1.) /rhoTot;

cout << "=============" << endl;
cout << "Omega_cdm 1: " << o1 << endl;
cout << "Omega_cdm 2: " << o2 << endl;
cout << "Omega_b    : " << ob << endl;
cout << "Omega_q    : " << oq << endl;
cout << "=============" << endl;

}

void DecayingDEDMCosmos::printStatus(ostream &o) {
  cout << "This is a DecayingDEDMCosmos, coupling is ";
  couple->printStatus(o);
  cout << endl;
  QuintCosmos::printStatus(o);
}
