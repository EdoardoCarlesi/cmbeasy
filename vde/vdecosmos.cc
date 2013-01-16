#include "vdecosmos.h"
#include "vectorde.h"
#include "global.h"
#include "arthur.h"
#include "arbitrary.h"
#include "cleanvector.h"
#include "cosmos.h"

#include <fstream>
#include <cmath>
#include <map>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "massiveneutrinos.h"

#include <fstream>
#include <map>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <cmath>
#include <limits>

#ifdef MACOSX_PANTHER
extern "C" int isnan(double); 
extern "C" int isinf(double); 
#endif 


VdeCosmos::VdeCosmos(): Cosmos(), vde(){
  //Omega_vde = 0;  // initialize with 0 dark energy
  //setOmega_vde(0);
}

VdeCosmos::~VdeCosmos() {
  delete vde;
}

void VdeCosmos::initVde(){
vde = new Vde(*this);
}

void VdeCosmos::initSplines() {
	Cosmos::initSplines();
      	Tau2Rho_v = new Spline(Tau2A, "Tau2rho_vde",&SplineAnchor);  
  Tau2Omega_v = new Spline(Tau2A, "Tau2omega_vde",&SplineAnchor);  
  Tau2Omega_vw = new Spline(Tau2A, "Tau2omega_vde",&SplineAnchor);  
  Tau2P_v = new Spline(Tau2A, "Tau2p_vde",&SplineAnchor);  
  Tau2W_v = new Spline(Tau2A, "Tau2w_vde",&SplineAnchor);  

  a2W = new Spline(1400,"a2w",&SplineAnchor);
  a2Omega_v = new Spline(a2W,"a2_omega_vde",&SplineAnchor);
  a2Omega_c = new Spline(a2W,"a2_omega_cdm",&SplineAnchor);
  a2Omega_b = new Spline(a2W,"a2_omega_baryons",&SplineAnchor);
  a2Omega_g = new Spline(a2W,"a2_omega_gamma",&SplineAnchor);
  a2Omega_l = new Spline(a2W,"a2_omega_lambda",&SplineAnchor);

  a2A_0 = new Spline(a2W,"a2_A_0",&SplineAnchor);
  a2Rho_v = new Spline(a2W,"a2_rho_vde",&SplineAnchor);
  a2Rho_c = new Spline(a2W,"a2_rho_cdm",&SplineAnchor);
  a2Rho_b = new Spline(a2W,"a2_rho_baryons",&SplineAnchor);
  a2Rho_g = new Spline(a2W,"a2_rho_gamma",&SplineAnchor);
  a2Rho_l = new Spline(a2W,"a2_rho_lambda",&SplineAnchor);
  a2dWdZ = new Spline(a2W,"a2dwdlnz",&SplineAnchor);
  a2OmegaWZ = new Spline(a2W,"a2Omegawz",&SplineAnchor);

  LogA2Omega_v = new Spline(1000,"logA2Omega_v",&SplineAnchor);
  LogA2Omega_vw = new Spline(LogA2Omega_v,"logA2Omega_vw",&SplineAnchor);
}

void VdeCosmos::dumpSplines(){
      	Tau2Rho_v->dump();
  Tau2Omega_v->dump(); 
  Tau2Omega_k->dump(); 
  Tau2Omega_vw->dump(); 
  Tau2P_v->dump();
  Tau2W_v->dump(); 

  a2W->dump(); 
  a2A_0->dump(); 
  a2Omega_c->dump(); 
  a2Omega_b->dump(); 
  a2Omega_v->dump(); 
  a2Omega_k->dump(); 
  a2Omega_g->dump(); 

  a2Rho_c->dump(); 
  a2Rho_b->dump(); 
  a2Rho_v->dump(); 
  a2Rho_k->dump(); 
  a2Rho_g->dump(); 
  a2dWdZ->dump(); 
  a2OmegaWZ->dump(); 

  LogA2Omega_v->dump(); 
  LogA2Omega_vw->dump();

}

void VdeCosmos::setOmega_vde(double x) {  
  Omega_vde = x;
}  

void VdeCosmos::propagateHistoryInTau(const double tau,const double* y,double* dy) {
//cout << " VdeCosmos::propagateHistory() " << endl;
  double a = y[1];
  double rho_g = y[2];
  double rho_b = y[3];
  double rho_c = y[4];
  double rho_nu = y[6];
  double rho_vde = y[8];
  double rho_k = y[9];
/* Update rho status in the VDE class */
  vde->set_rho(y[8]);

  double rho_nuNR = 0;
  double p_nuNR = 0.;
  if (amnu != 0.0) {
    double rnu, pnu;
    MassiveNeutrinos::nu1(a*amnu, &rnu, &pnu);
    rho_nuNR = rho_nu0() * nuNR() * rnu/(a*a*a*a);
    p_nuNR = rho_nu0() * nuNR() * pnu/(a*a*a*a);
  }
  
	rho_nuNR = nuNR() * massiveNeutrinoRho(0, a);
  
  double totalRho = rho_k + rho_g + rho_c + rho_b + rho_nu  + rho_nuNR + rho_vde;
  double H = a* sqrt(1.0/3.0*totalRho)/M_p();

  if (isnan(H)) {
    cout << "isnan H: " << tau << "a: " << a << endl;
    //cout << "vderho: " << vde->rho(a) << endl;
    throw Bad_Error("vdecosmos isnan H  !");
  }

  dy[1] = H*a;      // da/dtau is H*a
  dy[2] = -4.0*H*rho_g;  // gammas
  dy[3] = -3.0*H*rho_b;  // baryons
  dy[4] = -3.0*H*rho_c;  // cold dark matter

	/* Write here rho VDE equations*/
  dy[8] = -3.0*H*(1+vde->w(a))*rho_vde;
  dy[5] = a;   // dt/dtau
  dy[6] = -4.0*H*rho_nu;  // massless neutrinos

  double r = 3.0/4.0 * rho_b/rho_g;
  double cs = 1.0/sqrt(3*(1+r));  // sound speed
  dy[7] = cs;  // integrate for sound - horizon
}

void VdeCosmos::history(bool inform) {
if (validHistory()) {
    throw Bad_Error("VdeCosmos::history() I already have a history! Call reset() first");
  }
  getReady();
  rho__v = omega2rho(omega_v());

//inform=true;	
  cross_W25 = -40.6;

  double y[10];
  double ydot[10];
  initializeDensities(y, ydot);
  //if (vde != 0) vde->printStatus();
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
  double hnext = dtau*1e-3; // 1e-3;
//cout << "VdeCosmos::history()" << endl; 
 //for(int k=0; k<9; k++) cout << "y[" << k << "]" << y[k] << endl;


  while (y[1] < a_0()) {
    
//cout << y[1] << "  ---VdeCosmos::history(). Tau:" <<  tau << "dtau: " << dtau << y[1] << "hnext: " << hnext << endl ;
//hnext = Miscmath::odeint(y,11,tau,tau+dtau,1e-9,hnext,0,(moDerivs)&VdeCosmos::propagateHistoryInTauWrapper, *this);
hnext = Miscmath::odeint(y,9,tau,tau+dtau,1e-9,hnext,0,(moDerivs)&VdeCosmos::propagateHistoryInTauWrapper, *this);
//dtau=tau*1e-1;
//cout << y[1] << " VdeCosmos::history(). Tau:" <<  tau << "dtau: " << dtau << y[1] << "hnext: " << hnext << endl ;

    tau += dtau;

    propagateHistoryInTau(tau,y,ydot);   // We would like to have the tau derivatives for splines
    b = fillHistorySplines(tau,y,ydot);  // this also updates vdeessence to the right q and qdot
    A = b.A;

    for (int i=1; i<=7; i++) {
      if (isnan(y[i])) {
        cout << endl << "vdecosmos:: y["<<i<<"] isnan at time: " << tau << "  [ " << A << " ]"<<endl;
        throw Bad_Error("vdecosmos isnan !");
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
  if (inform) printStatus();
}

void VdeCosmos::initializeDensities(double* y, double* ydot)
{
//cout << "VdeCosmos::initializeDensities() " << endl;
  for (int k = 0; k < 9; k++) { y[k]=0; ydot[k] = 0;}
  y[1] = initialScaleFactor();  // we start from here (almost, cause we will note only from second step on)
  y[2] = initialRho_g(y[1]);
  y[3] = initialRho_b(y[1]);
  y[4] = initialRho_cdm(y[1]);
  y[5] = 0;
  y[6] = initialRho_nu(y[1]);
  y[7] = 0;
  y[8] = initialRhoVde(y[1]);
  y[9] = initialRho_k(y[1]);

  double totalRho_i = y[2] + y[3] + y[4] + y[6] + y[8] + y[9];
  double Hi = sqrt(1.0/3.0*totalRho_i)/M_p();

	/* Initialize and print initial scalar field value */ 
	vde->set_H_0(Hi);
	vde->setInitialA0();
double A0 = vde->getInitialA0();
//cout << "VdeCosmos::initializeDensities() " << endl;
cout << "vde->initialA0(): " << A0 << endl;

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
//cout << "VdeCosmos::iD_end() " << endl;
}

Background  VdeCosmos::fillHistorySplines(double tau, const double *y, const double *ydot, bool writeSplines)
{
  Background b = Cosmos::fillHistorySplines(tau, y, ydot, writeSplines);
  double A = b.A;
  // total rho and total pressure
  double totalRho = b.rho_g + b.rho_nu + b.rho_nuNR + b.rho_c + b.rho_b  + rho__v + vde->rho();
  double totalPressure = 1.0/3.0*(b.rho_g + b.rho_nu) + b.p_nuNR - rho__v + vde->p(A);
//cout << "totalRho: " << totalRho << endl;

  Tau2Rho_v->set(vde->rho());
  Tau2P_v->set(vde->p(A));
  Tau2W_v->set(vde->w(A));
  Tau2Omega_v->set(vde->rho() / totalRho);
//cout << "VdeCosmos::fillSplines()" << endl; 
  
  Tau2Omega_vw ->set(vde->rho() / totalRho*vde->w(A));

  // set eqn of state etc as function of logA (for Jan), also note crossing of
  // w_v = -0.25 if a > 0.01...
  LogA2Omega_v->set(log(A),vde->rho() / totalRho);
  LogA2Omega_vw->set(vde->rho() / totalRho / vde->w(A));
  if (vde->w(A) < -0.25 && cross_W25 == -4.6 && log(A) > -4.6)  cross_W25 = log(A);

  Tau2Rho->set(totalRho);
  Tau2P->set(totalPressure);

  // 8*pi*G*a^2*rho
  Tau2GRho->set(Gpi8()*A*A*totalRho);

//cout << "VdeCosmos::fillSplines()" << endl; 

  a2W->set(A,vde->w(A));
  // TODO a2A_0->set();
  a2Omega_v->set(vde->rho()/totalRho);
  a2Omega_c->set(b.rho_c/totalRho);
  a2Omega_b->set(b.rho_b/totalRho);
  a2Omega_g->set(b.rho_g/totalRho);
  a2Rho_v->set(vde->rho());
  a2Rho_c->set(b.rho_c);
  a2Rho_b->set(b.rho_b);
  a2Rho_g->set(b.rho_g);
 
  return b;
}

void VdeCosmos::finalizeHistorySplines(bool inform)
{
  Tau2A->arm(Spline::thoseReady);
  splineT2Tau->arm();
  splineA2Tau->arm(Spline::thoseReady);
  LogA2Omega_v->arm(Spline::thoseReady);
  a2W->arm();
  Z2iH_spline->flip();
  Z2iH_spline->arm(Spline::thoseReady);

  //calculate distance stuff
  double dzlum = distzmax / double(ndistpoints-1);
  for (double ztemp = 0.0; ztemp < distzmax; ztemp += (ztemp<3.?0.01:dzlum))
    Z2dist->set(ztemp, Z2iH_spline->integrate(0,ztemp));
  Z2dist->arm(Spline::thoseReady);

  // if inform is true, dump various background quantities to files
  if (inform) {
    
Z2iH_spline->dump("inverseH",false);
    Tau2A->dump("his");
    Tau2Rho->dump("rho");
    Tau2Rho_g->dump("rho_g");
    Tau2Rho_b->dump("rho_b");
    Tau2Rho_nu->dump("rho_nu");
    Tau2Rho_nuNR->dump("rho_nuNR");
    Tau2P_nuNR->dump("p_nuNR");
    Tau2w_nuNR->dump("w_nuNR");
    Tau2Rho_cdm->dump("rho_cdm");
    Tau2Rho_m->dump("rho_m");
    Tau2Rho_v->dump("rho_v");
    Tau2P->dump("pres");
    Tau2W_v->dump("w_vde");
    Tau2Omega_v->dump("omega_v");
  }

  Tau2Rho->derive(*Tau2RhoDot);
  Tau2RhoDot->arm();
  if (inform) Tau2RhoDot->dump("rhodot");

  Tau2Rho_g->derive(*Tau2RhoDot_g);
  Tau2RhoDot_g->arm();

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

   a2W->derive(*a2dWdZ,Spline::all);

   if (inform) a2dWdZ->dump("dwdz");

   for (int k = 0; k < a2dWdZ->size(); k++) {
     double z = a2dWdZ->x(k);
     double tau = z2tau(z);  // building blocks for the distinguishing of fast and slow vdeessence w-change
     a2OmegaWZ->set(pow(a2dWdZ->y(k),2) * Tau2Rho_v->fastY(tau)/tau2rho(tau) *z *z);
   }

   ValidOmega = true; 
}

double VdeCosmos::omega_vde() { if (validHistory()) return Tau2Rho_v->fastY(tau_0())/rho_crit();else return Omega_vde; }

void VdeCosmos::printStatus(ostream &o) {
  if (vde) vde->printStatus();
  o.setf(ios::scientific);
  o << endl << endl;
  o << " ============================================================="  << endl;
  o << " ======================= VdeCosmos Status =======================" << endl;
  o << " ============================================================="  << endl;
  o << " h : " << h() << endl;
  o << " Hubble's constant: " << H_0_cpm() << " Mpc^-1" << endl;
//  o << " Hence hubble time is : " << 1.0 / (   cpm2sInv()*tau2Hubble(tau_0()) ) << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " omega_phot: " << omega_g() <<      "    rho_g0:    " << rho_g0() << endl;
  o << " omega_nu  : " << omega_nu() <<     "    rho_nu0:   " << rho_nu0() << endl;
  o << " omega_nuNR: " << omega_nuNR() <<   "    rho_nuNR0: " << rho_nuNR0() << endl;
  o << " omega_cdm : " << omega_cdm() <<    "    rho_cdm0:  " << rho_cdm0()  << endl;
  o << " omega_b   : " << omega_b() <<      "    rho_b0:    " << rho_b0() << endl;
  o << " omega_vde : " << omega_vde() <<    "    rho_v0:    " << vde->rho() << endl;
//  o << " omega_kurv.   : " << omega_k() << endl;
  o << " --------------------------" << endl;
  o << " omega_0:      " << omega_0() << "  (the user-defined value)" << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " omega_nuNR h^2: " << omega_nuNR()*h2() << endl;
  o << " omega_cdm  h^2: " << omega_cdm()*h2() << endl;
  o << " omega_b    h^2: " << omega_b()*h2() <<  endl;
  o << " omega_m    h^2: " <<  omega_m()*h2() << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " sum of neutrino masses (in ev) : " << nuNR()*mass_nu_eV() << endl;

  //if (validHistory()) {

  if(true) {
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
   
    o << "::::  w_0: " <<  Tau2W_v->fastY(tau_0()) << endl;
    o << "::: Omega_v average till 0: " << Tau2Omega_v->average(tau_0()) << endl;
    o << "::: Omega_v structure f. : "<<  omesf(1.0/3.0) << endl; 
    o << endl << endl;
  }
}

double VdeCosmos::tau2AvOmega_v(const double tau) {
  return Tau2Omega_v->integrate(Tau2W_v->start(),tau)/(tau - Tau2W_v->start());
}

double VdeCosmos::sigma8Omega(WeffType t, int pn) {
  return sigma8Omega(InitialPower[pn],sigma8[pn],t);
}

double VdeCosmos::sigma8Omega(double n, double sig8, WeffType t) {
  double theta = n -1.0 + h()-0.65;
  double gamma_ = 0.21 - 0.22*weff(t) + 0.33*omega_m() + 0.25*theta;

  cout << "sig8om:: n" << n << " h: " << h() << endl;
  cout << "sig8om:: O_m : " << omega_m() << endl;
  cout << "sig8om:: s8: " << sig8 << endl;

  return sig8 * pow(omega_m(),gamma_);
}

void VdeCosmos::setOmega_vde_flat() {
	//cout << "set_Omega_vde_flat() " << endl; 
     Omega_vde = 1.-omega_cdm()-omega_b()-omega_g();
	//cout << "cdm, b, g: " << omega_cdm() << " " << omega_b() << " " << omega_g() << " " << endl;
     if (nuR()>0)
        Omega_vde-=omega_nu();
     if (nuNR()>0)
       Omega_vde-=omega_nuNR(false);
     if (Omega_vde<0.)
       throw Bad_Error("Cosmos::setOmega_vuint_flat() - would need negative Omega_cdm");
     setOmega_vacuum(0.);
     setOmega_vde(Omega_vde);

     vde->set_OmegaA0(Omega_vde);
     vde->set_OmegaM0(omega_cdm()+omega_b());

    cout << "VdeCosmos::OmegaVde set to: " << Omega_vde << endl;
  } //!< set omega_vuintessence such that a flat universe comes out

double VdeCosmos::initialRhoVde(double a_0){
double rho_A0;
double result=0;
void *p;

int dim = 1000;
double a_1=1;
//a_0=0.00001;

double *steps = new double[dim];
double s = -log(a_0/a_1)/(double (dim-1));
double aa=a_0;

for (int i=0; i<dim; i++){
aa=a_0*exp(s*i);
steps[i]=aa;
}
//cout << "steps[0,1]: " << steps[0] << " " << steps[1] << endl;

double a_m=a_0;

//Cout << "VdeCosmos::initializeRhoVde() " << endl;

//double step = 1./(double (dim-1));// + a_0;
double step; 

//cout << " step: " << step << endl; 

for (int j=1; j<dim; j++){
step = steps[j] - steps[j-1];
a_m = a_0 + steps[j-1] + 0.5*step;
result += step * (1./a_m)*(1+vde->w(a_m)); //*(1+1./3.);//vde->w(a_m));
}

/* The integral calculates the log in base e*/
result*=(1./log(10));

//cout << "VdeCosmos, result : " << result << endl;
delete [] steps; 
double rho_t = rho_0()*vde->get_OmegaA0();
rho_A0 = rho_t*pow(10,3*result);
cout << "VdeCosmos, rho_t, rho_A0: " << rho_t << " " << rho_A0 << endl; 
vde->set_rhoA0(rho_A0);

// TODO !!!! why 20 or 15 or 16??
//return 16*rho_A0;	
return 15.85*rho_A0;	
}

/*!
  Equation (3) from astro-ph 0107525 
*/
double VdeCosmos::estimateSigma8Q(double sig8 , double tau0) {

  if (omega_vde() == 0.0) return sig8;

  double tautr = a2tau(exp(cross_W25));
  double odsf =  tau2AvOmega_v(tautr);
  

  double e1 = 3.0/5.0 *  odsf;
  double e2 = -0.2 * (1 + 1.0/weff(logAWeff));
  
  double s8 = sig8;

  s8 *= pow(a_equ(), e1);
  s8 *= pow(1 - omega_v(), e2);
  s8 *= sqrt(tau_0() / tau0);
  
  return s8;
}

/*!
 Return \int_ln(a_equ)^ln(a_tr) Omega_dark(a) dln a  
 divided by  ln(a_tr) -ln(a_equ)
 
 Here a_tr is the transition scale factor at which the
 dark energy becomes dominant. Take your pick :-)
 See astro-ph/0107525 for definition etc.

*/
double VdeCosmos::omesf(double a_tr) {
  const double a_eq = a_equ();
  if (isnan(a_eq)) return std::numeric_limits<double>::quiet_NaN();
  double integral =  LogA2Omega_v->integrate(log(a_eq),log(a_tr));
  LogA2Omega_v->dump("logaomega");
  return integral / (log(a_tr) - log(a_equ()));
}


