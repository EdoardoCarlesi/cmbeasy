#include "quintcosmos.h"
#include "global.h"
#include "exponential.h"
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

#include "massiveneutrinos.h"

#ifdef INTERACTING
#include "interactingquint.h"
#endif

#ifndef PRERELEASE
#include "corasaniti.h"
#include "tristar.h"
#include "omegastep.h"
#include "binomega.h"
#include "arbitraryomega.h"
#include "constantomega.h"
#include "cuscutan.h"
#endif

#include <fstream>
#include <map>
#include <math.h>
#include <cmath>
#include <limits>

#ifdef MACOSX_PANTHER
extern "C" int isnan(double); 
extern "C" int isinf(double); 
#endif 


QuintCosmos::QuintCosmos(Quintessence::Type quinttype) :   Cosmos() , quint(0) {
  Omega_quintessence = 0;  // initialize with 0 dark energy
  setQuintessence(quinttype);
  setOmega_quintessence(0);
}

QuintCosmos::~QuintCosmos() {
  delete quint;
}

/*!
  Only Quintessence stuff
*/
void QuintCosmos::initSplines() {
  Cosmos::initSplines();

  Tau2Rho_qpot = new Spline(Tau2A, "Tau2rhoqpot",&SplineAnchor);  Tau2Rho_qpotX = -1.0;
  Tau2Rho_qkin =  new Spline(Tau2A, "Tau2rhoqkin",&SplineAnchor);  Tau2Rho_qkinX = -1.0;
  Tau2Rho_q = new Spline(Tau2A, "Tau2rhoq",&SplineAnchor);  Tau2Rho_qX = -1.0;
  Tau2P_q = new Spline(Tau2A, "Tau2pq",&SplineAnchor);  Tau2P_qX = -1.0;

  Phi2Tau = new Spline(Tau2A->size(), "Phi2tau",&SplineAnchor); Phi2TauX = -1.0;
  Tau2Phi = new Spline(Tau2A, "Tau2phi",&SplineAnchor); Tau2PhiX = -1.0;
  Tau2PhiDot = new Spline(Tau2A, "Tau2phidot",&SplineAnchor); Tau2PhiDotX = -1.0;
  Tau2W_q = new Spline(Tau2A, "Tau2wq",&SplineAnchor); Tau2W_qX = -1.0;  
  Tau2WDot_q = new Spline(Tau2A, "Tau2dotwq",&SplineAnchor); Tau2WDot_qX = -1.0;

  Tau2Vprime = new Spline(Tau2A,"Tau2Vprime",&SplineAnchor); Tau2VprimeX= -1.0;
  Tau2DotVprime = new Spline(Tau2A,"Tau2DotVprime",&SplineAnchor); Tau2DotVprimeX= -1.0;

  Tau2VprimeOverRho = new Spline(Tau2A,"Tau2VoverRho",&SplineAnchor); Tau2VprimeOverRhoX= -1.0;
  Tau2DotVprimeOverRho = new Spline(Tau2A,"Tau2VoverRho",&SplineAnchor); Tau2DotVprimeOverRhoX= -1.0;

  Tau2P_qDot = new Spline(Tau2A, "Tau2P_qDot", &SplineAnchor); Tau2P_qDotX = -1.0;
  Tau2Rho_qDot = new Spline(Tau2A, "Tau2Rho_qDot", &SplineAnchor); Tau2Rho_qDotX = -1.0;

  Tau2Omega_q = new Spline(Tau2A, "Tau2OmegaQ",&SplineAnchor);
  Tau2Omega_qw = new Spline(Tau2A, "Tau2OmegaQwQ",&SplineAnchor);
  Tau2Mass_factor = new Spline(Tau2A,"Tau2mass_factor",&SplineAnchor);

  A2Omega_q = new Spline(splineA2Tau,"a2omegaq",&SplineAnchor);
  A2Hubble = new Spline(splineA2Tau,"a2hubble",&SplineAnchor);
  A2Hubble_factor = new Spline(splineA2Tau,"a2hubble_factor",&SplineAnchor);
  A2Mass_factor = new Spline(splineA2Tau,"a2mass_factor",&SplineAnchor);
  A2Omega_qw = new Spline(splineA2Tau,"a2omegaqw",&SplineAnchor);
  A2Omega_qkin =  new Spline(splineA2Tau, "a2omega_qkin",&SplineAnchor);

  Z2W = new Spline(Z2iH_spline,"z2w",&SplineAnchor);
  Z2Hubble_factor = new Spline(Z2iH_spline,"z2hubble_factor",&SplineAnchor);
  Z2Hubble_ratio = new Spline(Z2iH_spline,"z2hubble_ratio",&SplineAnchor);
  Z2dWdZ = new Spline(Z2iH_spline,"z2dwdlnz",&SplineAnchor);
  Z2OmegaWZ = new Spline(Z2iH_spline,"z2Omegawz",&SplineAnchor);
  Z2Omega_qkin =  new Spline(Z2iH_spline, "z2omega_qkin",&SplineAnchor);

  LogA2Omega_q = new Spline(1000,"logA2omegaq",&SplineAnchor);
  LogA2Omega_qw = new Spline(LogA2Omega_q,"logA2omegaqw",&SplineAnchor);


#define RESETSPLINE(a) a##X=-1.23456995233e-245

  RESETSPLINE(Tau2Rho_qpot);
  RESETSPLINE(Tau2Rho_qkin);
  RESETSPLINE(Tau2Rho_q);
  
  RESETSPLINE(Tau2Vprime);
  RESETSPLINE(Tau2DotVprime);

  RESETSPLINE(Tau2VprimeOverRho);
  RESETSPLINE(Tau2DotVprimeOverRho);

  RESETSPLINE(Tau2P_q);
  RESETSPLINE(Phi2Tau);
  RESETSPLINE(Tau2Phi);
  RESETSPLINE(Tau2PhiDot);
  RESETSPLINE(Tau2W_q);
  RESETSPLINE(Tau2WDot_q);

  RESETSPLINE(Tau2P_qDot);
  RESETSPLINE(Tau2Rho_qDot);
#undef RESETSPLINE
}

/*!
  reset quintessence (if any) and call cosmos::reset()
*/
void QuintCosmos::reset() { 
  //cout << "Quintcosmos:: reset() " << endl;
  quint->reset();
  Cosmos::reset();
}





/*!
  Create new quintessence. Usually this is done only once.

  ATTENTION: if omega_q() is zero only your wish will
  be stored in userRequest. Instead of your wish, 
  a voidQuintessence() will be constructed.
  Only when you set omega_quintessence to none-zero
  values, will your wish be fulfilled and automatically
  your quintessence will be constructed.

  This is in fact a very similar routine than what PerturbationFactory
  does. Maybe, one should move this to Quintessence and 
  make it a QuintessenceFactory...

*/
void QuintCosmos::setQuintessence(Quintessence::Type quinttype) {
  userRequest = quinttype;

  if (omega_q(false) == 0)
    quinttype = Quintessence::none;

  delete quint;
  quint =0;

  switch (quinttype) {
  case Quintessence::exponential: 
    quint = new Exponential(*this);
    break;
  case Quintessence::ipl: 
    quint = new Ratra(*this);
    break;
  case Quintessence::leaping:
    quint = new Arthur(*this);
    break;
#ifdef INTERACTING
  case Quintessence::interacting: 
    quint = new InteractingQuint(*this);
    break;
#endif
  case Quintessence::arbitrary:
    quint = new Arbitrary(*this);
    break;
  case Quintessence::celestine:
    quint = new Celestine(*this);
    break;
  case Quintessence::crossover:
    quint = new Crossover(*this);
    break;
  case Quintessence::crossoverfield:
    quint = new CrossoverField(*this);
    break;
#ifndef PRERELEASE
 case Quintessence::corasaniti:
    quint = new Corasaniti(*this);
    break;
  case Quintessence::tristar:
    quint = new Tristar(*this);
    break;
  case Quintessence::omegastep:
    quint = new OmegaStep(*this);
    break;
  case Quintessence::binomega:
    quint = new BinOmega(*this);
    break;
  case Quintessence::arbitraryomega:
    quint = new ArbitraryOmega(*this);
    break;
  case Quintessence::constantomega:
    quint = new ConstantOmega(*this);
    break;
#endif     
  default: 
    quint = new VoidQuintessence(*this);
  }
}


/*!
  Sets omega_quintessence.

  If this is zero, it will destroy the Quintessence quint and
  create a voidQuintessence() instead.

  If it has been zero and now will have some none
  zero value. Your wish Quintessence (which
  you may already have specified using setQuintessence()
  or QuintCosmos(QuintType) will be created
*/
void QuintCosmos::setOmega_quintessence(double x) {  
  if (omega_q() == 0) {
    if (x !=0) {
    
      if (quint) { delete quint; quint=0;}
      Omega_quintessence = x;
  
      setQuintessence(userRequest);
     
    } // else: if it is zero and stays zero do nothing
  } else {
    if (x ==0) {
      if (quint) { delete quint; quint=0;}
      quint = new VoidQuintessence(*this);
    }  // else if it has been none-zero and stazs none-zero, do nothing   
  }
  Omega_quintessence = x;
}  

void QuintCosmos::tuneQuintessence(double Omebar, double Weff ) {
  if (!quint) throw Bad_Error("QuintCosmos::tuneQuintessence() There is no Quintessence yet. Call setQuintessence() first");
  quint->tuneQuintessence(Omebar);
  return;
}

/*! 
  Given the array y[], fill the array dy[] with the 
  derivatives w.r.t conformal time tau. 

  As always, output of tau is in Mpc, i.e. d/dtau is in Mpc^-1 
*/
void QuintCosmos::propagateHistoryInTau(const double tau,const double* y,double* dy) {
  double a = y[1];
  double rho_g = y[2];
  double rho_nu = y[6];
  double rho_b = y[3];
  double rho_c = y[4];

  quint->setQ(y[8]);
  quint->setQDot(y[9]);

  double rho_nuNR = 0;
  double p_nuNR = 0.;
  if (amnu != 0.0) {
    double rnu, pnu;
    MassiveNeutrinos::nu1(a*amnu, &rnu, &pnu);
    rho_nuNR = rho_nu0() * nuNR() * rnu/(a*a*a*a);
    p_nuNR = rho_nu0() * nuNR() * pnu/(a*a*a*a);
  }
  //double rho_nuNR = nuNR() * massiveNeutrinoRho(0, a);

  
  double totalRho = rho_g + rho_c + rho_b + rho_nu  + rho_nuNR + rho__v + quint->rho(a);
  double H = a* sqrt(1.0/3.0*totalRho)/M_p();

  if (isnan(H)) {
    cout << "isnan H: " << tau << "a: " << a << endl;
    cout << "quintrho: " << quint->rho(a) << endl;
    throw Bad_Error("quintcosmos isnan H  !");
  }

  dy[1] = H*a;      // da/dtau is H*a
  dy[2] = -4.0*H*rho_g;  // gammas
  dy[3] = -3.0*H*rho_b;  // baryons
  dy[4] = -3.0*H*rho_c;  // cold dark matter


  if (quint->needsPropagation()) {
    dy[8] = quint->qDot(a);  // quintessence
    dy[9] = -2*H*quint->qDot(a) - a*a*quint->Vprime(quint->q(a),a,H);     // quintessence
  } else {
    dy[8] = 0;
    dy[9] = 0;
  }
  dy[5] = a;   // dt/dtau
  dy[6] = -4.0*H*rho_nu;  // massless neutrinos

  double r = 3.0/4.0 * rho_b/rho_g;
  double cs = 1.0/sqrt(3*(1+r));  // sound speed
  dy[7] = cs;  // integrate for sound - horizon
}

void QuintCosmos::history(bool inform) {
  if (validHistory()) {
    throw Bad_Error("QuintCosmos::history() I already have a history! Call reset() first");
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
    hnext = Miscmath::odeint(y,9,tau,tau+dtau,1e-9,hnext,0,(moDerivs)&QuintCosmos::propagateHistoryInTauWrapper, *this);
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

  findSpecialMoments();
  
finalizeHistorySplines(inform);

  if (inform)
    printStatus();

}

void QuintCosmos::initializeDensities(double* y, double* ydot)
{
  for (int k = 0; k < 9; k++) { y[k]=0; ydot[k] = 0;}
  y[1] = initialScaleFactor();                         // we start from here (almost, cause we will note only from second step on)
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


Background  QuintCosmos::fillHistorySplines(double tau, const double *y, const double *ydot, bool writeSplines)
{
  Background b;
  double A = b.A = y[1];
  b.rho_g = y[2];
  b.rho_nu = y[6];
  b.rho_b = y[3];
  b.rho_c = y[4];
  b.t = y[5];
  b.rs = y[7];

  double plain =0;
  double rnu=0, pnu=0;
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
    Z2H_spline->set(a2z(A), ydot[1]/(A*A));
/*
  Background b = Cosmos::fillHistorySplines(tau, y, ydot, writeSplines);
  double A = b.A;
*/

  double Rho_crit = 3*pow(H_0_cpm()*M_p(),2); // Today's critical density
  omega_cdm_final = omega_cdm(false);
  double omega_b_final = omega_b(false);
  double dm0 = Rho_crit*omega_cdm_final; 
  double b0 = Rho_crit*omega_b_final; 
 // total rho and total pressure
  double totalRho = b.rho_g + b.rho_nu + b.rho_nuNR + b.rho_c + b.rho_b  + rho__v + quint->rho(A);
  double totalPressure = 1.0/3.0*(b.rho_g + b.rho_nu) + b.p_nuNR - rho__v + quint->p(A);
  double h_factor = totalRho/Rho_crit;
  double h_factor_z = (sqrt(1.0/3.0*totalRho)/M_p())*h()*100*h()*100;
  double m_factor = b.rho_c*pow(A,3)/dm0;
  double b_factor = b.rho_b*pow(A,3)/b0;
	
  double h_factor_lcdm = (1-omega_cdm_final-omega_b(false)) + (omega_cdm_final+omega_b(false))*pow(A,-3);

//    Tau2Rho_g->set(b.rho_g);
//    Tau2Rho_b->set(b.rho_b);
//    Tau2Rho_m->set(rho_m);
//    Tau2Rho_cdm->set(b.rho_c);

  Tau2Rho_qpot->set(quint->V(quint->q(A), A));
  Tau2Rho_qkin->set(quint->rho_kin(A));
  Tau2Rho_q->set(quint->rho(A));
  Tau2Phi->set(quint->q(A));
  Tau2PhiDot->set(quint->qDot(A));

 double rho_qkin_units = quint->rho_kin(A)/totalRho;
  Z2Omega_qkin->set(rho_qkin_units);
  A2Omega_qkin->set(rho_qkin_units);
  // A2W.set(A,quint->w(A));
  A2Hubble->set(sqrt(1.0/3.0*totalRho)/M_p());
  Tau2W_q->set(quint->w(A));
  Tau2Omega_q->set(quint->rho(A) / totalRho);
  Tau2Omega_qw ->set(quint->rho(A) / totalRho*quint->w(A));

  Tau2Vprime->set(quint->Vprime(quint->q(A),A,ydot[1]/A));
  Tau2VprimeOverRho->set(quint->Vprime(quint->q(A),A,ydot[1]/A)/quint->rho(A));

  LogA2Omega_q->set(log(A),quint->rho(A) / totalRho);
  LogA2Omega_qw ->set(quint->rho(A) / totalRho / quint->w(A));
  if (quint->w(A) < -0.25 && cross_W25 == -4.6 && log(A) > -4.6)  cross_W25 = log(A);

  Tau2Rho->set(totalRho);
  Tau2P->set(totalPressure);
  Tau2P_q->set(quint->p(A));

  Tau2GRho->set(Gpi8()*A*A*totalRho);

  A2Omega_q->set(quint->rho(A) / totalRho);
  A2Omega_qw ->set(quint->rho(A) / totalRho*quint->w(A));
  A2Hubble_factor->set(h_factor);
  Z2Hubble_factor->set(h_factor_z);
  Z2Hubble_ratio->set(h_factor/h_factor_lcdm);
  Tau2Rho_cdm->set(b.rho_c);
  A2Mass_factor->set(m_factor);
  //A2Mass_factor->set(b_factor);
  Tau2Mass_factor->set(m_factor);
  Z2W->set(quint->w(A));

  return b;
}

void QuintCosmos::finalizeHistorySplines(bool inform)
{
  Tau2A->arm(Spline::thoseReady);
  splineT2Tau->arm();
  splineA2Tau->arm(Spline::thoseReady);
  LogA2Omega_q->arm(Spline::thoseReady);

  Z2W->flip(); // we flip this cause we do not want 10^6 ... 0  but the other way around (right order)
  Z2iH_spline->flip();
    Z2Hubble_factor->flip();  
    Z2Hubble_ratio->flip();

  Z2iH_spline->arm(Spline::thoseReady);

  //calculate distance stuff
  double dzlum = distzmax / double(ndistpoints-1);
  for (double ztemp = 0.0; ztemp < distzmax; ztemp += (ztemp<3.?0.01:dzlum))
    Z2dist->set(ztemp, Z2iH_spline->integrate(0,ztemp));
  Z2dist->arm(Spline::thoseReady);

  // if inform is true, dump various background quantities to files
  if (inform) {

	dumpExtraSplines();

 }

  Tau2Rho->derive(*Tau2RhoDot);
  Tau2RhoDot->arm();
  if (inform) Tau2RhoDot->dump("rhodot");

  Tau2Rho_g->derive(*Tau2RhoDot_g);
  Tau2RhoDot_g->arm();

  Tau2Rho_q->derive(*Tau2Rho_qDot);
  Tau2Rho_qDot->arm();

  Tau2P_q->derive(*Tau2P_qDot);
  Tau2P_qDot->arm();

  Tau2Rho_nu->derive(*Tau2RhoDot_nu);
  Tau2RhoDot_nu->arm();

  Tau2Rho_nuNR->derive(*Tau2RhoDot_nuNR);
  Tau2RhoDot_nuNR->arm();


  Tau2w_nuNR->derive(*Tau2wdot_nuNR);
  Tau2wdot_nuNR->arm();

  Tau2Rho_b->derive(*Tau2RhoDot_b);
  Tau2RhoDot_b->arm();

  Tau2Rho_cdm->derive(*Tau2RhoDot_cdm);
  Tau2RhoDot_cdm->arm();

  Tau2Rho_m->derive(*Tau2RhoDot_m);
  Tau2RhoDot_m->arm();

  Tau2GRho->derive(*Tau2GRhoDot);
  Tau2GRhoDot->arm();

  Tau2Vprime->derive(*Tau2DotVprime);
  Tau2DotVprime->arm();


  Tau2VprimeOverRho->derive(*Tau2DotVprimeOverRho);
  Tau2DotVprimeOverRho->arm();

  Tau2W_q->derive(*Tau2WDot_q);
  Tau2WDot_q->arm();

   Z2W->derive(*Z2dWdZ,Spline::all);

   if (inform) Z2dWdZ->dump("dwdz");

   for (int k = 0; k < Z2dWdZ->size(); k++) {
     double z = Z2dWdZ->x(k);
     double tau = z2tau(z);  // building blocks for the distinguishing of fast and slow quintessence w-change
     Z2OmegaWZ->set(pow(Z2dWdZ->y(k),2) * tau2rho_q(tau)/tau2rho(tau) *z *z);
   }

   Z2OmegaWZ->arm();

   ValidOmega = true; // new
   quint->postPrepare(); // maybe some things need to be done with the knowledge of history ?
   if (quint->needsToUpdatePhi()) { // uups, is Phi not the real one ? (only for Arbitrary, so far)
     //     Tau2Phi->dump("phibefore");
     Tau2Phi->disarm();
     Tau2PhiDot->disarm();
     for (int i=0; i < Tau2Phi->size(); i++) {
       double tau=Tau2Phi->x(i);
       double a = tau2a(tau);
       Tau2Phi->setY(i,quint->q(a));
       Tau2PhiDot->setY(i,quint->qDot(a));
     }
     Tau2Phi->arm();
     Tau2PhiDot->arm();
   }

   for (int i=0; i < Tau2Phi->size(); i++) {
     double tau = Tau2Phi->x(i);
     double phi = Tau2Phi->y(i);
     Phi2Tau->set(phi, tau);
   }
   Phi2Tau->makeProper();
   if (Phi2Tau->size()>10) {
     Phi2Tau->arm();
   }
}

void QuintCosmos::dumpExtraSplines(){
	cout << "dumpExtraSplines()" << endl;
    Z2iH_spline->dump("inverseH",false);
    Tau2A->dump("his");
    Tau2ADot->dump("adot");
    Tau2Rho->dump("rho");
    Tau2Rho_g->dump("rho_g");
    Tau2Rho_b->dump("rho_b");
    Tau2Rho_nu->dump("rho_nu");
    Tau2Rho_nuNR->dump("rho_nuNR");
    Tau2P_nuNR->dump("p_nuNR");
    Tau2w_nuNR->dump("w_nuNR");
    Tau2Rho_cdm->dump("rho_cdm");
    Tau2Rho_m->dump("rho_m");

    Tau2Phi->dump("quint");
    Tau2PhiDot->dump("quintdot");
    Tau2Rho_qpot->dump("rho_qpot");
    Tau2Rho_qkin->dump("rho_qkin");
    Tau2Rho_q->dump("rho_q");
    Tau2P_q->dump("pres_q");
    Tau2P->dump("pres");
    Tau2W_q->dump("wq");
    Tau2Omega_qw->dump("woq");
    Tau2Omega_q->dump("omega_q");

    Z2Omega_qkin->dump("z_omega_qkin");
    A2Omega_qkin->dump("a_omega_qkin");
    A2Hubble->dump("a_hubble");  
    A2Hubble_factor->dump("hubble_factor");  
    Z2Hubble_factor->dump("hubble_factor_z");  
    Z2Hubble_ratio->dump("hubble_ratio_to_lcdm");  
    A2Mass_factor->dump("mass_factor");  
    Tau2Mass_factor->dump("tau_mass_factor");  

}


/*!
  Slight subtelty:
  
  The original exp() Quintessence uses omega_q() to set the parameter lambda in a deterministic way. Now, the field,
  its energy etc. are calculated using the *wish* omega_q() from the very beginning. 

  The fluctuation evolution, however would later on use a (potentially) different omega_q(), thus leading to wrong
  spectra. It is therefore possible to call omega_q(false) which will give the *wish* back. This is used by
  Quintessence Exponential 
*/
double QuintCosmos::omega_q(bool reality) { if (validHistory() && reality) return tau2rho_q(tau_0())/rho_crit();else return Omega_quintessence; }

void QuintCosmos::printStatus(ostream &o) {
  
	double tt1, tt2;
	double tau1, tau2;

	tau1 = z2tau(59.7);
	tau2 = z2tau(60.4);
	//tau1 = 1744.3475538;
	//tau2 = 1744.185439;
	tt1 = 1.e-09*pow(mpc2year(),1)*tau2t(tau1);
	tt2 = 1.e-09*pow(mpc2year(),1)*tau2t(tau2);
	//tt2 = 

 	cout << " ************************************** " << endl;
	cout << " tt1 : " << tt1  << " tt2 : " << tt2  << " tt1-tt2  : " << tt1-tt2 << endl;
 	cout << " tau1: " << tau1 << " tau2: " << tau2 << " tau1-tau2: " << tau1-tau2 << endl;
	cout << " ************************************** " << endl;
	
	if (quint) quint->printStatus();
  o.setf(ios::scientific);
  o << endl << endl;
  o << " =================================================================="  << endl;
  o << " ======================= QuintCosmos Status =======================" << endl;
  o << " =================================================================="  << endl;
  o << " h : " << h() << endl;
  o << " Hubble's constant: " << H_0_cpm() << " Mpc^-1" << endl;
  o << " Hence hubble time is : " << 1.0 / (   cpm2sInv()*tau2Hubble(tau_0()) ) << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " omega_phot: " << omega_g() <<      "    rho_g0:    " << rho_g0() << endl;
  o << " omega_nu  : " << omega_nu() <<     "    rho_nu0:   " << rho_nu0() << endl;
  o << " omega_nuNR: " << omega_nuNR() <<   "    rho_nuNR0: " << rho_nuNR0() << endl;
  o << " omega_cdm : " << omega_cdm() <<    "    rho_cdm0:  " << rho_cdm0()  << endl;
  o << " omega_m   : " << omega_m()   << endl;
  o << " omega_b   : " << omega_b() <<      "    rho_b0:    " << rho_b0() << endl;
  o << " omega_quint: " << omega_q() <<      "    rho_q0:    " << tau2rho_q(tau_0()) << endl;
  o << " omega_kurv.   : " << omega_k() << endl;
  o << " --------------------------" << endl;
  o << " omega_0:      " << omega_0() << "  (the user-defined value)" << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " omega_nuNR h^2: " << omega_nuNR()*h2() << endl;
  o << " omega_cdm  h^2: " << omega_cdm()*h2() << endl;
  o << " omega_b    h^2: " << omega_b()*h2() <<  endl;
  o << " omega_m    h^2: " <<  omega_m()*h2() << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " sum of neutrino masses (in ev) : " << nuNR()*mass_nu_eV() << endl;

if (validHistory())
//if(false) 
 {
    o << " -------------------------------------------------------------"  << endl;
    o << " ---------------conformal times & scale factors  -------------" << endl;
    o << " -------------------------------------------------------------"  << endl;
    o << " equality::        tau = "  << tau_equ() << "   a = " << a_equ() << endl;  
    o << " last scattering:: tau = "  << tau_ls() << "   a = " << a_ls()  << endl;  
    o << " today::           tau = "  << tau_0() << "   a = " << a_0() << endl;  
    o << " z_ls: " << z_ls() << endl;
   
    o << " -------------------------------------------------------------"  << endl;
    o << " optdlss: " << optDistanceLss() << endl; 
    o << " reionizationFraction " << reionizationFraction() << endl;
    o << " -------------------------------------------------------------"  << endl;
    o << " sound horizons: " << endl;
    o << " -------------------------------------------------------------"  << endl;
    o << " equality: " << tau2rs(tau_equ()) << " Mpc" << endl;
    o << " ls      : " << tau2rs(tau_ls()) << " Mpc" << endl;
    o << " -------------------------------------------------------------"  << endl;
 
    o << " time equ:   " << mpc2year()*tau2t(tau_equ())  << " years" << endl;
    o << " time ls:    " << mpc2year()*tau2t(tau_ls()) << " years" << endl;
    o << " time today: " << mpc2year()*tau2t(tau_0()) << " years" << endl;
    o << " ------------------------------------------------------------"  << endl;
   
    o << "::::  w_0: " <<  tau2w_q(tau_0()) << endl;
    o << "::: Omega_q average till ls: " << tau2AvOmega_q(tau_ls()) << endl;
    o << "::: Omega_q average till 0: " << Tau2Omega_q->average(tau_0()) << endl;
    o << "::: Omega_q structure f. : "<<  omesf(1.0/3.0) << endl; 
    o << " ------------------------------------------------------------"  << endl;
    o << "::: RhoCrit(z=0)  : " << 3*pow(tau2Hubble(tau_0())*M_p(),2) << endl;
    o << "::: RhoCrit(a=1.)  : " << 3*pow(tau2Hubble(splineA2Tau->fastY(1.))*M_p(),2) << endl;
    o << "::: RhoCrit_0(z=0): " << 3*pow(H_0_cpm()*M_p(),2) << endl;
    o << "::: RhoTot (z=0)  : " << Tau2Rho->fastY(tau_0()) << endl;

    o << " ------------------------------------------------------------"  << endl;
    o << endl << endl;
  }
}

double QuintCosmos::tau2AvOmega_q(const double tau) {
  return Tau2Omega_q->integrate(Tau2W_q->start(),tau)/(tau - Tau2W_q->start());
}


double QuintCosmos::sigma8Omega(WeffType t, int pn) {
  return sigma8Omega(InitialPower[pn],sigma8[pn],t);
}

double QuintCosmos::sigma8Omega(double n, double sig8, WeffType t) {
  double theta = n -1.0 + h()-0.65;
  double gamma_ = 0.21 - 0.22*weff(t) + 0.33*omega_m() + 0.25*theta;

  cout << "sig8om:: n" << n << " h: " << h() << endl;
  cout << "sig8om:: O_m : " << omega_m() << endl;
  cout << "sig8om:: s8: " << sig8 << endl;

  return sig8 * pow(omega_m(),gamma_);
}

/*!
  Equation (3) from astro-ph 0107525 
*/
double QuintCosmos::estimateSigma8Q(double sig8 , double tau0) {

  if (omega_q(false) == 0.0) return sig8;

  double tautr = a2tau(exp(cross_W25));
  double odsf =  tau2AvOmega_q(tautr);
  

  double e1 = 3.0/5.0 *  odsf;
  double e2 = -0.2 * (1 + 1.0/weff(logAWeff));
  
  double s8 = sig8;

  s8 *= pow( a_equ(), e1);
  s8 *= pow( 1 - omega_q(), e2);
  s8 *= sqrt( tau_0() / tau0);
  
  return s8;
}

/*!
 Return \int_ln(a_equ)^ln(a_tr) Omega_dark(a) dln a  
 divided by  ln(a_tr) -ln(a_equ)
 
 Here a_tr is the transition scale factor at which the
 dark energy becomes dominant. Take your pick :-)
 See astro-ph/0107525 for definition etc.

*/
double QuintCosmos::omesf(double a_tr) {
  const double a_eq = a_equ();
  if (isnan(a_eq)) return std::numeric_limits<double>::quiet_NaN();
  double integral =  LogA2Omega_q->integrate(log(a_eq),log(a_tr));
  //  LogA2Omega_q->dump("logaomega");
  return integral / (log(a_tr) - log(a_equ()));
}


