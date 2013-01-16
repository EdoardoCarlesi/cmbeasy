#include "cosmos.h"

#include "global.h"
#include "cleanvector.h"
#include "recfast.h"
#include "recfastalpha.h"
#include "recombinationfactory.h"
#include "massiveneutrinos.h"
#include "perturbation.h"
#include "mass_function.h"

#include <fstream>
#include <map>
#include <limits>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <algorithm>

const double Cosmos::distzmax  = 2.0;// Maximum redshift luminosity distance calculated to
const int Cosmos::ndistpoints = 200;// Number of points in luminosity distance

Cosmos::Cosmos() : baseCosmos(), mInitialScaleFactor(1e-12), ValidThermo(false), TimeStep(0)
{
  // for the sake of clearity set here and not at construction:

  setInitialPower(1.0);  // default only one spectral index, which is 1.0
  reset(); // get all splines new etc

  setReionizationZ(0);
  setOptDistanceLss(0);

  h_h = 0.65;      // small Hubble h
  Y_P = 0.24;      // primordial Helium mass fraction

  rho__v = 0; // initialize vacuum energy
  tau__0 = Rho_crit = std::numeric_limits<double>::quiet_NaN();
  mNrNonRelativisticNu = 0;
  Omega_nuNR = 0;
  NrRelativisticNu = 3.04;
  mMassNu_eV = amnu = 0;

  setOmegaH2_b(0.02);
  setOmegaH2_vacuum(0);
  setOmegaH2_nu(1.695e-5);
  setOmegaH2_g(2.488e-5);
  setOmegaH2_nuNR(0);

  Omega_cdm = 1 - omega_b() - omega_v() - omega_nu() - omega_nuNR()  - omega_g();  // flat universe


  T_cmB = 2.726;   // Temperature of cosmic background radiation
}


/*!
  Init most of the splines that cosmos needs. Most of them are actually children of just
  a few mothers. 
  We give these mothers our SplineAnchor which in turn we will ask on reset() to
  kill those mothers and hence the children too. [sounds cruel, sorry :-]
*/
void Cosmos::initSplines() {
  splineCs2 = new Spline(10000,"cs2",&SplineAnchor);
  Expmmu = new Spline(splineCs2,"emmu",&SplineAnchor);
  Visibility = new  Spline(splineCs2,"visibility",&SplineAnchor);
  DVisibility = new  Spline(splineCs2,"dvisibility",&SplineAnchor);
  DDVisibility = new  Spline(splineCs2,"ddvisibility",&SplineAnchor);
  FreeElectronFraction = new Spline(splineCs2,"emmu",&SplineAnchor);

  //DOpac = new  Spline(splineCs2,"dopac visibility",&SplineAnchor);

  Opac =new Spline(splineCs2,"dotmu",&SplineAnchor);
  D2Kappa = new Spline(splineCs2,"dotmu",&SplineAnchor);

  Tau2A = new Spline(2000,"Tau2A",&SplineAnchor);  Tau2AX = -1.0;  
  Tau2ADot = new Spline(Tau2A, "Tau2ADot",&SplineAnchor); Tau2ADotX = -1.0;
  Tau2Rho = new Spline(Tau2A, "Tau2Rho",&SplineAnchor); Tau2RhoX = -1.0;
  Tau2Rho_g =  new Spline(Tau2A, "Tau2Rho_g",&SplineAnchor); Tau2Rho_gX = -1.0;
  Tau2Rho_nu = new Spline(Tau2A, "Tau2Rho_nu",&SplineAnchor);  Tau2Rho_nuX = -1.0;
  Tau2Rho_nuNR = new Spline(Tau2A, "Tau2Rho_nuNR",&SplineAnchor); Tau2Rho_nuNRX = -1.0;
  Tau2P_nuNR = new Spline(Tau2A, "Tau2pnunr",&SplineAnchor); Tau2P_nuNRX = -1.0;
  Tau2w_nuNR = new Spline(Tau2A, "Tau2wnunr",&SplineAnchor); Tau2w_nuNRX = -1.0;
  Tau2wdot_nuNR = new Spline(Tau2A, "Tau2wnunrDot",&SplineAnchor); Tau2wdot_nuNRX = -1.0;

  Tau2Rho_m = new Spline(Tau2A, "Tau2rhom",&SplineAnchor);  Tau2Rho_mX = -1.0;
  Tau2Rho_b =  new Spline(Tau2A, "Tau2rhob",&SplineAnchor);  Tau2Rho_bX = -1.0;
  Tau2Rho_cdm = new Spline(Tau2A, "Tau2rhocdm",&SplineAnchor);  Tau2Rho_cdmX = -1.0;
 
  Tau2RhoDot =  new Spline(Tau2A, "Tau2rhodot",&SplineAnchor); Tau2RhoDotX = -1.0;
  Tau2RhoDot_g = new Spline(Tau2A, "Tau2rhodotg",&SplineAnchor); Tau2RhoDot_gX = -1.0;
  Tau2RhoDot_m = new Spline(Tau2A, "Tau2rhodotm",&SplineAnchor); Tau2RhoDot_mX = -1.0;
  Tau2RhoDot_cdm = new Spline(Tau2A, "Tau2rhodotcdm",&SplineAnchor); Tau2RhoDot_cdmX = -1.0;
  Tau2RhoDot_b = new Spline(Tau2A, "Tau2rhodotb",&SplineAnchor); Tau2RhoDot_bX = -1.0;
  Tau2RhoDot_nu = new Spline(Tau2A, "Tau2rhodotnu",&SplineAnchor); Tau2RhoDot_nuX = -1.0;
  Tau2RhoDot_nuNR = new Spline(Tau2A, "Tau2rhodotnuNR",&SplineAnchor); Tau2RhoDot_nuNRX = -1.0;

  Tau2GRho =  new Spline(Tau2A, "Tau2grho",&SplineAnchor); Tau2GRhoX = -1.0;
  Tau2GRhoDot = new Spline(Tau2A, "Tau2grhodot",&SplineAnchor); Tau2GRhoDotX = -1.0;

  Tau2P =  new Spline(Tau2A, "Tau2p",&SplineAnchor); Tau2PX = -1.0;

  Tau2R = new Spline(Tau2A, "Tau2r",&SplineAnchor);
  Tau2T = new Spline(Tau2A, "Tau2t",&SplineAnchor);
  Tau2Rs = new Spline(Tau2A, "Tau2rs",&SplineAnchor);
  Tau2Cs2 = new Spline(Tau2A, "Tau2cs2",&SplineAnchor);
  Tau2Cs = new Spline(Tau2A, "Tau2cs",&SplineAnchor);
  Tau2ADotToA =  new Spline(Tau2A, "Tau2ADotToA",&SplineAnchor);

  splineA2Tau = new Spline(1000,"a2tau",&SplineAnchor);
  splineT2Tau =new Spline(1000,"t2tau",&SplineAnchor);
 
  Z2iH_spline = new Spline(1000,"z2iH",&SplineAnchor);
  Z2H_spline = new Spline(1000,"z2H",&SplineAnchor);
  Z2dist = new Spline(ndistpoints,"Z2dist",&SplineAnchor);

  // in the following, we set the caching arguments to some funny
  // value that will never occur, in this world (hopefully). Note that 
  // this is just a precaution, as you would have to call some caching
  // function before and after calling Cosmos::reset() with exactly
  // the same argument. 
  // but, to be on the safe Side

#define RESETSPLINE(a) a##X=-1.23456995233e-245
  RESETSPLINE(Visibility);
  RESETSPLINE(DVisibility);
  RESETSPLINE(DDVisibility);
  RESETSPLINE(Tau2A);
  RESETSPLINE(Tau2ADot);
  RESETSPLINE(Opac);
  RESETSPLINE(D2Kappa);
  RESETSPLINE(Expmmu);

  RESETSPLINE(Tau2Rho);
  RESETSPLINE(Tau2Rho_g);
  RESETSPLINE(Tau2Rho_nu);
  RESETSPLINE(Tau2Rho_nuNR);
  RESETSPLINE(Tau2Rho_m);
  RESETSPLINE(Tau2Rho_b);
  RESETSPLINE(Tau2Rho_cdm);

  RESETSPLINE(Tau2RhoDot);
  RESETSPLINE(Tau2RhoDot_b); 
  RESETSPLINE(Tau2RhoDot_cdm); 
  RESETSPLINE(Tau2RhoDot_m); 
  RESETSPLINE(Tau2RhoDot_g);
  RESETSPLINE(Tau2RhoDot_nu);
  RESETSPLINE(Tau2RhoDot_nuNR);


  RESETSPLINE(Tau2GRho);
  RESETSPLINE(Tau2GRhoDot);

  RESETSPLINE(Tau2w_nuNR);
  RESETSPLINE(Tau2wdot_nuNR);

  RESETSPLINE(Tau2P_nuNR);
  RESETSPLINE(Tau2P);

  RESETSPLINE(Tau2R);
  RESETSPLINE(Tau2T);
  RESETSPLINE(Tau2Rs);
  RESETSPLINE(Tau2Cs2);
  RESETSPLINE(Tau2ADotToA);

#undef RESETSPLINE


  // 
  // Now, let us init the cdm_power splineWeb setPower_cdm() automatically adds this to the SplineAnchor
  //

 
  setPower_cdm("power_cdm");   // this web stores the cdm power spectrum  (with the intial spectrum factored out)
  setPower_baryon("power_baryon");
  setPower_gamma("power_gamma");
  setPower_nu("power_nu");
  setPower_nuNR("power_nuNr");
  setPower_psi("power_psi");
}

void Cosmos::reset() {
  ValidHistory = false; ValidThermo = false; ValidOmega = false;
  
  //The rebounding stuff will have to be changed if non-flat Universes
  // are ever added to Cosmos.
  //Right now the Universe can never be rebounding because this can't
  // happen in the flat case
  isRebounding = false;


  SplineAnchor.kill(); //! Get rid of all splines that are known to splineAnchor
  for (unsigned int i=0; i < InitialPower.size(); i++)  PowerNormalization[i] = 1;   // power normalization back to usual
  initSplines();
  mExtraTimeSteps.clear();
}


double Cosmos::tau(double t) {
  return 0;
}

double Cosmos::cpm2msun(){
// Msun in proton masses 1.189e57
// Proton mass 0.938 Gev
double Gev2Msun=1./(0.938*1.189e+57);
return cpm2Gev(1)*Gev2Msun;
}


double Cosmos::mpc2s(int n) const {
  if (n ==1) return  1.0292454e14; else return 1.0292454e14*mpc2s(n-1); // recursive calls
}

double Cosmos::cpm2sInv(int n) const {
  if (n ==1) return 9.7158559e-15; else return 9.7158559e-15*cpm2sInv(n-1);
}

// From PDG: hbar*c = 197.326 9631 MeV fm
//           1 pc   = 3.085 677 580 7 × 10^16 m
//           --> conversion from GeV to Mpc^-1 is 1.0/(hbar*c*1e-3*1e-15/pc*1e-6)
static const double gev2cmp_const = 1.56373844315248e+38;

/*! Convert Gev^n -> Mpc^-n */
double Cosmos::Gev2cpm(int n) const {
  if (n ==1) return gev2cmp_const; else return  gev2cmp_const*Gev2cpm(n-1);
} 

double Cosmos::cpm2Gev(int n) const  {
   if (n ==1) return 1./gev2cmp_const; else return cpm2Gev(n-1)/gev2cmp_const;
} 

double Cosmos::M_p(int n) const {
  if (n)  {  
    if (n > 0) return M_p() * M_p(n-1);
    return M_p(n+1)/M_p();  
  } 
  return 1.0;
}

/*!
  Given the array y[], fill the array dy[] with the 
  derivatives w.r.t conformal time tau.

  The following indices correspond to the following
  quantities:

  1 -> a
  2 -> rho_gamma
  3 -> rho_baryon
  4 -> rho_cdm

  5 -> t (usual time)
  6->  rho_nu (massless)
  7 -> soundhorizon

  As always, output of tau is in Mpc, i.e. d/dtau is in Mpc^-1
*/
void Cosmos::propagateHistoryInTau(const double tau,const double* y,double* dy)
{
  double a = y[1];
  double rho_g = y[2];
  double rho_nu = y[6];
  double rho_b = y[3];
  double rho_c = y[4];

  double rho_nuNR = 0;
  if (amnu != 0.0) {
    double rnu, pnu;
    MassiveNeutrinos::nu1(a*amnu,&rnu,&pnu);
    rho_nuNR = rho_nu0() * nuNR() * rnu/(a*a*a*a);
  }
//X   double rho_nuNR = nuNR() * massiveNeutrinoRho(0, a);

  double totalRho = rho_g + rho_c + rho_b + rho_nu  + rho_nuNR + rho__v;
  double H = a* sqrt(1.0/3.0*totalRho)/M_p();

  dy[1] = H*a;
  dy[2] = -4.0*H*rho_g;
  dy[3] = -3.0*H*rho_b;
  dy[4] = -3.0*H*rho_c;

  dy[5] = a;   // dt/dtau
  dy[6] = -4.0*H*rho_nu;

  double r = 3.0/4.0 * rho_b/rho_g;
  double cs = 1.0/sqrt(3*(1+r));
  dy[7] = cs;

}

void Cosmos::printOutputList(char* fname, double yr_step, double z0, double z1){

cout << "Cosmos:printOutputList for " << fname << ", step: " << yr_step << " z_min: " << z0 << " z_max: " << z1 << ", total age: " << ageYr << endl;

	double yr;
	double z, a, t, tau;
	int TOT=Tau2A->size();
	//FILE *out=fopen("output_list.txt","w");
	FILE *out=fopen(fname,"w");

	a2year = new Spline(TOT,"a2year");
	year2a = new Spline(TOT,"year2a");

	for(int k=0; k<TOT; k++) {
	tau = Tau2A->start(k);
	a = Tau2A->front(k);
	yr = mpc2year()*tau2t(tau); 
	a2year->set(a,yr);
	year2a->set(yr,a);
}
	//yr_step = 1.2e+8;
	year2a->arm();
	yr=0;

do{
 yr += yr_step;
 a = year2a->fastY(yr);
 z = 1./a - 1;
//	cout << "Year: " << yr << " a: " << a << endl;
if(z<z1 && z>z0) fprintf(out, "%lf\n", a);
if(a>1.) fprintf(out, "%lf\n", 1.);
} while (yr<ageYr);

	fclose(out);
}

void Cosmos::history(bool inform) {
  if (validHistory()) throw Bad_Error("Cosmos::history() I already have a history! Call reset() first");
  getReady();
  rho__v = omega2rho(omega_v());

  if (amnu !=0.0)
    MassiveNeutrinos::initnu1(amnu);

  double y[11];
  double ydot[11];

  y[1] =  1e-12;                         // we start from here (almost, cause we will note only from second step on)
  y[2] = initialRho_g(y[1]);
  y[6] = initialRho_nu(y[1]);
  y[3] = initialRho_b(y[1]);
  y[4] = initialRho_cdm(y[1]);
  y[5] = 0;
  y[7] = 0;

  // now we want to get da/dtau in order to estimate the beginning tau....
  propagateHistoryInTau(0,y,ydot);  
  double tau = y[1]/ydot[1];  // a / adot

  if (inform) cout << "from initial step, we have: " << tau << endl;

  double dtau = tau*1e-1;  //1   // and dtau is one tenth of this..

  // zeroing tau__equ and tau__ls as they will now be (re)calculated
  tau__equ = 0;
  tau__ls = 0;

  double hnext = dtau*1e-3;
  while (y[1] < a_0()) {
    hnext = Miscmath::odeint(y,7,tau,tau+dtau,1e-9,hnext,0,(moDerivs)&Cosmos::propagateHistoryInTauWrapper, *this);
    tau += dtau;

    propagateHistoryInTau(tau,y,ydot);     // We would like to have the tau derivatives for splines
    Background b = fillHistorySplines(tau,y,ydot);  // fill the spline with background information
    double A = b.A;
    // total rho and total pressure
    double totalRho = b.rho_g + b.rho_nu + b.rho_nuNR + b.rho_c + b.rho_b + rho__v ;
    double totalPressure = 1.0/3.0*(b.rho_g + b.rho_nu) + b.p_nuNR - rho__v;
    Tau2Rho->set(totalRho);
    Tau2P->set(totalPressure);
    Tau2GRho->set(Gpi8()*A*A*totalRho);

    // determine tau__equ by comparing densities
    // this is rough and will be refined by findSpecialMoments()
    if (b.rho_b + b.rho_c > b.rho_g && tau__equ == 0) tau__equ = tau;    
    // determine t_LS by noting the crossing of a above the a_ls level...
    // rough, as above, should be refined in findSpecialMoments()
    if (A>= a_ls() && tau__ls == 0) tau__ls = tau;
    if (A < 0.1*a_0()) dtau *= 1.05;  // increase the stepsize 
    dtau = min(10.0,dtau); // level off with stepsize
//    cout << scientific << setw(10) << "computed history to A=" << A
//         << " (rho_nuNR=" << b.rho_nuNR << ")" << '\r';
  }

  Tau2A->arm(Spline::thoseReady);
  splineT2Tau->arm();
  splineA2Tau->arm(Spline::thoseReady);

  Z2iH_spline->flip();  // this is in reverse order, so flip
  Z2iH_spline->arm(Spline::thoseReady);
 
  Z2H_spline->flip();  // this is in reverse order, so flip
  Z2H_spline->arm(Spline::thoseReady);

  //Fill in the distance measure spline
  double dzlum = distzmax / double(ndistpoints-1);
  for (double ztemp = 0.0; ztemp < distzmax; ztemp += (ztemp<3.?0.01:dzlum))
    Z2dist->set(ztemp, Z2iH_spline->integrate(0,ztemp));
  Z2dist->arm(Spline::thoseReady);

  findSpecialMoments();  // determine a_ls, tau_ls etc accurately

	ageYr = mpc2year()*tau2t(tau_0());
  Tau2Rho->derive(*Tau2RhoDot);
  Tau2RhoDot->arm();

  Tau2Rho_g->derive(*Tau2RhoDot_g);
  Tau2RhoDot_g->arm();

  Tau2Rho_b->derive(*Tau2RhoDot_b);
  Tau2RhoDot_b->arm();

  Tau2Rho_cdm->derive(*Tau2RhoDot_cdm);
  Tau2RhoDot_cdm->arm();

  Tau2Rho_m->derive(*Tau2RhoDot_m);
  Tau2RhoDot_m->arm();

  Tau2Rho_nu->derive(*Tau2RhoDot_nu);
  Tau2RhoDot_nu->arm();
  
  Tau2Rho_nuNR->derive(*Tau2RhoDot_nuNR);
  Tau2RhoDot_nuNR->arm();


  Tau2w_nuNR->derive(*Tau2wdot_nuNR);
  Tau2wdot_nuNR->arm();
  
  Tau2GRho->derive(*Tau2GRhoDot);
  Tau2GRhoDot->arm();


  ValidOmega = true; // new

  // if you ever  like to see the behaviour of this background quantities, just
  // call cosmos with inform=true
   if (inform) { 
     string base = ControlPanel::cmbeasyDir("/output/");
     Tau2A->dump(base+"his");
     Tau2ADot->dump(base+"adot");
     Tau2Rho->dump(base+"rho");
     Tau2Rho_g->dump(base+"rho_g");
     Tau2Rho_b->dump(base+"rho_b");
     Tau2Rho_nu->dump(base+"rho_nu");
     Tau2Rho_nuNR->dump(base+"rho_nuNR");
     Tau2P_nuNR->dump(base+"p_nuNR");
     Tau2w_nuNR->dump(base+"w_nuNR");
     Tau2Rho_m->dump(base+"rho_m");     
     Tau2Rho_nu->dump(base+"rho_nu");

     Tau2RhoDot->dump(base+"rhodot");
     Tau2Rho_nu->dump(base+"rho_nu",false);
     Tau2RhoDot_nu->dump(base+"rhodot_nu",false);
     Tau2RhoDot_m->dump(base+"rhodot_m",false);
     Tau2RhoDot->dump(base+"rhodot",false);
     Tau2RhoDot_g->dump(base+"rhodot_g",false);

     Tau2GRhoDot->dump(base+"grhodot");

     Z2iH_spline->dump(base+"Z2iH_spline");
     Z2H_spline->dump(base+"Z2H_spline");
     Z2dist->dump(base+"Z2dist");
  }
}

void Cosmos::setExtraBackgroundTimesteps(const std::vector<double>& extraSteps)
{
  mExtraTimeSteps = extraSteps;
  std::sort(mExtraTimeSteps.rbegin(), mExtraTimeSteps.rend());
}

void Cosmos::checkForExtraTimestep(double tau, double& dtau)
{
  double& nextStep = mExtraTimeSteps.back();
  if (tau+dtau >= nextStep) {
    dtau = nextStep-tau;
    mExtraTimeSteps.pop_back();
  }
}

/*!
  This is split from history(), cause re-implementations in sub-classes, such
  as QuintCosmos also need to set these splines. They just call this
  function, then. 
  A specific structure called Background is returned holding some quantities
  that the sub-class does not need to re-calculate, cause it is already done here
*/
Background  Cosmos::fillHistorySplines(double tau, const double *y, const double *ydot, bool writeSplines) {
  Background b;
  double A = b.A = y[1];
  b.rho_g = y[2];
  b.rho_nu = y[6];
  b.rho_b = y[3];
  b.rho_c = y[4];
  b.t = y[5];
  b.rs = y[7];

//X   b.rho_nuNR = nuNR() * massiveNeutrinoRho(0, A);
//X   b.p_nuNR = nuNR() * massiveNeutrinoPressure(0, A);

//X   for (int i = 1; i <= 7; ++i) {
//X     if (isinf(y[i]) || isnan(y[i])) {
//X         std::stringstream s;
//X         s << "y[" << i << "] is " << (isinf(y[i])?"infinite":"nan") << " at tau= " << tau << " a= " << b.A << " and z=" << (1./b.A - 1.);
//X         throw Bad_Error(s.str());
//X     }
//X   }

  double plain =0;
  double rnu=0, pnu=0;
  if (amnu !=0.0) {
    MassiveNeutrinos::nu1(b.A*amnu,&rnu,&pnu);
    plain = rho_nu0() * nuNR() /(A*A*A*A);
  }

  b.rho_nuNR = rnu * plain;
  b.p_nuNR = pnu*plain;

  if (writeSplines)
  {
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
    Z2H_spline->set(a2z(A), ydot[1]/(A*A));
  }

  return b;
}


/*!
  Finds matter-radiation equality and last scatterting and sets
  the relevant variables
*/
void Cosmos::findSpecialMoments() {
  //cout << "Cosmos::findSpecialMoments()\n";

  tau__0 =  a2tau(a_0());     // tau__0 easy as tau(a_0())

  // here, we set the critical mass density at tau_0(), for we want to
  // access this when calculating Omega's
  Rho_crit = 3.0*pow(tau2Hubble(tau__0)*M_p(),2);  //both agree well
  Rho_crit = tau2rho(tau_0());

  tau__ls = a2tau(a_ls());     // tau__ls easy as tau(a_ls())

  ValidHistory = true;  // we have to set it here and not before, cause omega_b() etc will access the spline from now on

  // tau equ by findEqualtity() returning difference of rho_m and rho_g
  // absolut precision tau_ls()*1e-10 means approx rel. prec. 1e-9 
  tau__equ = Miscmath::zbrent( (moSingle)&Cosmos::findEquality, *this, tau_min(), tau_max(), tau__ls*1e-10);

  t__ls = tau2t(tau__ls);
  t__equ = tau2t(tau__equ);  
  t__0 = a2t(a_0());

  a__equ = tau2a(tau__equ);

  //cout << "Cosmos::~findSpecialMoments()\n";
}

/*!
  Used by Miscmath::zbrent() to find exact conformal time of 
  equality
*/
double Cosmos::findEquality(const double tau) {
  const double relDiff = fabs(tau2rho_m(tau) -  (tau2rho_cdm(tau) + tau2rho_b(tau))) / tau2rho_m(tau);
  if (  relDiff > 1e-7) {
    std::stringstream relDiffString;
    relDiffString << relDiff << " at tau= " << tau << " and z=" << tau2z(tau);
    throw Bad_Error("Cosmos::findEquality() - cdm mismatch; rel difference at " + relDiffString.str());
  }
  return tau2rho_m(tau)-tau2rho_relativistic(tau);
}



/*!
  return rho of all relativistic degrees of freedom (meaning:
  those with mass smaller than kbT()  roughly so...
*/
double Cosmos::tau2rho_relativistic(const double tau) {
  double result = tau2rho_g(tau) + tau2rho_nu(tau);
  if (omega_nuNR()) {
    if (massiveNeutrinosAreRelativistic(tau))
      result += tau2rho_nuNR(tau);
  }
  return result;
}

bool Cosmos::massiveNeutrinosAreRelativistic(const double tau)
{
  return (mass_nu_cpm() < kbT(T_nu()/tau2a(tau)));
}


//All of these distance measures assume flatness currently
// Because that is all that cosmos currently supports
double Cosmos::propermotionDistance(double z) const {
  if ( Z2dist->isWithinBounds(z) ) return Z2dist->fastY(z);
  return Z2iH_spline->integrate(0,z);
}

double Cosmos::propermotionDistance(double z0, double z1) const {
  return Z2iH_spline->integrate(z0,z1);
}

double Cosmos::angulardiameterDistance(double z) const {
  if ( Z2dist->isWithinBounds(z) ) return Z2dist->fastY(z)/(1+z);
  return Z2iH_spline->integrate(0,z)/(1+z);
}

double Cosmos::luminosityDistance(double z) const {
  if ( Z2dist->isWithinBounds(z) ) return (1+z)*Z2dist->fastY(z);
  return (1+z)*Z2iH_spline->integrate(0,z);
}

double Cosmos::comovingVolume(double z0, double z1){
double cV;
cV = (1./h())*3000*H_0_cpm()*Z2iH(z1)*
(propermotionDistance(z0, z1))*(propermotionDistance(z0, z1));
//cout << "ComovingVolume: " << cV*1e-9 << " Gpc^3  "<< endl;
//cout << "H_0/H(z): " << H_0_cpm()*Z2iH(z1)  << endl;
return cV;
}

double Cosmos::integrateComovingVolume(double z0, double z1){
double z, step_z, cv, vol1, vol2; int size = 200;
Spline *integral = new Spline(size, "integral", &SplineAnchor);

//step_z = (z1-z0)/size;
//z=z0;
step_z = (z1-0)/size;
z=0;

for(int i=0; i<size; i++){
z += step_z;
//cv = comovingVolume(z0, z);
cv = comovingVolume(0, z);
integral->set(z,cv);
}
integral->arm();
vol1 = integral->integrate(0,z0);
vol2 = integral->integrate(0,z1);
delete integral;

// TODO make these coordinates modifiable from the exterior
double the1 = 0; 
double the2 = 180;
double phi1 = 0;
double phi2 = 360;

// Return solid angle value
Miscmath *mm;
double angle = mm->solid_angle(the1,the2,phi1,phi2);

cout << "Total ComovingVolume: " << (vol2-vol1)*angle << endl;
return (vol2-vol1)*angle;
}

/*!
  Unlike the angular diameter distance, the luminosity distance 
  cares about time dilation, and hence the factor outside the
  integral is 1 + zhel and not 1 + zcmb.  The integral, however,
  is done in the isotropic frame, and hence is from 0 to zcmb.
 */
double Cosmos::luminosityDistance(double zhel, double zcmb) const {
  if ( Z2dist->isWithinBounds(zcmb) ) return (1+zhel)*Z2dist->fastY(zcmb);
  return (1+zhel)*Z2iH_spline->integrate(0,zcmb);
}


double Cosmos::physical2comovingK(const double k, const double tau) {
  return k * tau2a(tau);
}

double Cosmos::z_equ() {
  return a2z(a_equ());
}

/*!
  Convert rho to Omega at some time tau. If you do not specify tau,
  tau_0() is taken. The feature of tau dependent omega is only
  available, if there is a validHistory()
*/
double Cosmos::rho2omega(const double rho,  double tau) { 
  if (validHistory()) {  // if we have a history, take the calculated values 
    if (tau == -1) tau = tau_0(); 
    return rho/(3.0*pow(tau2Hubble(tau)*M_p(),2)); 
  }
  return rho/(3.0*pow(H_0_cpm()*M_p(),2)); 
}
double Cosmos::omega2rho(double omega) { return omega * 3.0 * pow(H_0_cpm()*M_p(),2); }

double Cosmos::omega_b(bool fromHistory) { if (validOmega() && fromHistory) return tau2rho_b(tau_0())/rho_crit(); else return Omega_b; }
double Cosmos::omega_cdm(bool FromHistory) { if (validOmega() && FromHistory) return tau2rho_cdm(tau_0())/rho_crit(); else return Omega_cdm; }
double Cosmos::omega_v() { if (validOmega()) return tau2rho_v()/rho_crit(); else return Omega_vacuum; }

double Cosmos::omega_g() { if (validOmega()) return  tau2rho_g(tau_0())/rho_crit();else return Omega_g; }
double Cosmos::omega_nu() { if (validOmega()) return tau2rho_nu(tau_0())/rho_crit(); else return Omega_nu; }
double Cosmos::omega_nuNR(bool fromHistory) { if (validOmega() && fromHistory) return tau2rho_nuNR(tau_0())/rho_crit(); else return Omega_nuNR; }


double Cosmos::tau2H(const double tau) { return tau2adot(tau)/tau2a(tau); }
double Cosmos::tau2Hubble(const double tau) { return tau2adot(tau)/pow(tau2a(tau),2); }


double Cosmos::rho_nuNR0() {
  return 3*omega_nuNR()*pow(M_p()*H_0_cpm(),2);
}

//! factor for converting between Omega of massive neutrinos today and mass_nu in eV
static double massNuFactor(const Cosmos& c)
{
  //see e.g. Dodelson, p. 46 (...and exercise 18)
  static const double constFactor=7./180.*M_PI*M_PI*M_PI*M_PI;
  static const double zeta3=1.20205690;
  const double omegaMasslessNu = c.rho_nu0()/c.rho_0();
  const double kbt_conversion = c.Gev2cpm()*1e-9/c.kbT(c.T_nu());
  return constFactor/(zeta3*kbt_conversion*c.h2()*omegaMasslessNu);
}

void Cosmos::setOmegaNuNRFromNeutrinoMass()
{
  const double fac = massNuFactor(*this);
  setOmega_nuNR( nuNR()*mMassNu_eV/(fac*h2()), DoNotAdjustNeutrinoMass);
}

void Cosmos::setNeutrinoMassFromOmegaNuNR()
{
  if(nuNR()==0) {
    if (Omega_nuNR>0) {
      throw Bad_Error("Cosmos::setNeutrinoMassFromOmegaNuNR() - non-zero Omega_nuNR, but no massive neutrino species.");
    } else {
      mMassNu_eV = 0;
    }
    return;
  }
  const double fac = massNuFactor(*this);
  mMassNu_eV = fac*Omega_nuNR*h2()/nuNR();
}


void Cosmos::printStatus(const char* name) {
  ofstream o(name);
  printStatus(o);
}

void Cosmos::printStatus(ostream &o) {

  o.setf(ios::scientific);
  o << endl << endl;
  o << " ============================================================="  << endl;
  o << " ======================= Cosmos Status =======================" << endl;
  o << " ============================================================="  << endl;
  o << " h : " << h() << endl;
  o << " Hubble's constant: " << H_0_cpm() << " Mpc^-1" << endl;
  if (validHistory()) {
    o << " Hence hubble time is : " << 1.0 / (   cpm2sInv()*tau2Hubble(tau_0()) ) << endl;
  }
  o << " -------------------------------------------------------------" << endl;
  o << " omega_phot: " << omega_g() <<      "    rho_g0:    " << rho_g0() << "  Mpc^-4" <<  endl;
  o << " omega_nu  : " << omega_nu() <<     "    rho_nu0:   " << rho_nu0() << "  Mpc^-4" <<endl;
  o << " omega_nuNR: " << omega_nuNR() <<   "    rho_nuNR0: " << rho_nuNR0()  << "  Mpc^-4"<<endl;
  o << " omega_cdm : " << omega_cdm() <<    "    rho_cdm0:  " << rho_cdm0()  << "  Mpc^-4"<<endl;
  o << " omega_b   : " << omega_b() <<      "    rho_b0:    " << rho_b0() << "  Mpc^-4" <<endl;
  o << " omega_Lambda: " << omega_v() <<      "    rho_v0:    " << tau2rho_v() << "  Mpc^-4"<< endl;
  o << " omega_kurv.   : " << omega_k() << endl;
  o << " --------------------------" << endl;
  o << " omega_0:      " << omega_0() << "  (the user-defined value)" << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " omega_cdm  h^2: " << omega_cdm()*h2()  << endl;
  o << " omega_b    h^2: " << omega_b()*h2()    << endl;
  o << " omega_m    h^2: " << omega_m()*h2()    << endl;
  o << " omega_nuNR h^2: " << omega_nuNR()*h2() << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " sum of neutrino masses (in eV) : " << nuNR()*mass_nu_eV() << endl;
  if (validHistory()) {
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
    o << " equality: " << ((tau_equ()!=0)?tau2rs(tau_equ()):-1) << " Mpc" << endl;
    o << " ls      : " << ((tau_ls()!=0)?tau2rs(tau_ls()):-1) << " Mpc" << endl;
    o << " -------------------------------------------------------------"  << endl; 
    double t_eq, t_ls, t_0;
    if (tau_equ()!=0) t_eq = mpc2year()*tau2t(tau_equ()); else t_eq = -1;
    if (tau_ls()!=0) t_ls = mpc2year()*tau2t(tau_ls()); else t_ls = -1;
    if (tau_0()!=0) t_0 = mpc2year()*tau2t(tau_0()); else t_0 = -1;

    o << " time equ:   " << t_eq << " years" << endl;
    o << " time ls:    " << t_ls << " years" << endl;
    o << " time today: " << t_0  << " years" << endl;
    o << " ------------------------------------------------------------"  << endl;
    o << endl << endl;
  }
}


//==============  massive neutrinos ==================================================

//! Calculate the pressure of the massive neutrino species nr at scale factor a
//! Ma & Bertschinger Eq. 52
double Cosmos::massiveNeutrinoPressure(int nr, double a)
{
  static const double points = 100;

  double gsoverhp3 = twooverh3();

  Spline qIntegral(100, "Cosmos::massiveNeutrinoPressure::qIntegral");

  static const double qstart = 3e-2;
  static const double qend = 30.;
  const static double lnstep = (log(qend)-log(qstart))/points;

  const double k_bT = kbT(T_nu());
  const double m_nu = mass_nu_kbT();

  for (double q=qstart; q <= qend; q = exp(log(q)+lnstep))
  {
    const double epsilon = sqrt(1.+a*a*m_nu*m_nu/(q*q));
    const double f_0 = 1./(exp(q)+1.);
    qIntegral.set(q, q*q*q/epsilon * f_0);
  }

  qIntegral.arm();

  return 1./3.*4.*M_PI*pow(k_bT, 4)*pow(a, -4.)*gsoverhp3*qIntegral.integrate();
}

//! 2/h_p^3, where h_p is Planck's constant in units of
double Cosmos::twooverh3() const
{
  static const double h_P_eVs = 4.13566743e-15;  // 2002 CODATA value for Planck's constant in eV s
  const double h_P = Gev2sInv(h_P_eVs*1e-9);

  double gsoverhp3 = 1./pow(h_P, 3.);
  // old factor for comparison, remove
  static const double factor= 4.930625e+00/3.04;
  gsoverhp3 *= factor;

  return gsoverhp3;
}


//! Calculate the energy density of the massive neutrino species nr at scale factor a
//! Ma & Bertschinger Eq. 52
double Cosmos::massiveNeutrinoRho(int nr, double a)
{
  static const double points = 100;

  const double gsoverhp3 = twooverh3();

  Spline qIntegral((int)points, "Cosmos::massiveNeutrinoRho::qIntegral");

  static const double qstart = 3e-2;
  static const double qend = 30.;
  static const double lnstep = (log(qend)-log(qstart))/points;

  const double k_bT = kbT(T_nu()); // T_nu() is T_nu today
  const double m_nu = mass_nu_kbT();

  for (double q=qstart; q <= qend; q = exp(log(q)+lnstep))
  {
    const double epsilon = sqrt(1.+a*a*m_nu*m_nu/(q*q));
    const double f_0 = 1./(exp(q)+1.);
    qIntegral.set(q, q*q*q*epsilon*f_0);
  }


  qIntegral.arm();
//X   double max = qIntegral.maximum();
//X   static double lastmax = 0;
//X   if (max != lastmax)
//X   cout << endl << max << endl;
//X   lastmax = max;

  return 4.*M_PI*pow(k_bT, 4)*pow(a, -4.)*gsoverhp3*qIntegral.integrate();
}

double Cosmos::initialRho_nuNR(const double a, double* p_nuNR)
{
  double rho_nuNR, pNuNR;
  if (amnu != 0.0) {
    static bool isInitialized=false;
    if (!isInitialized) {
      MassiveNeutrinos::initnu1(amnu);
      isInitialized=true;
    }
    double rnu, pnu;
    MassiveNeutrinos::nu1(a*amnu,&rnu,&pnu);
    rho_nuNR = rho_nu0()*nuNR()*rnu/(a*a*a*a);
    pNuNR = rho_nu0()*nuNR()*pnu/(a*a*a*a);
  } else {
    rho_nuNR = 0.;
    pNuNR = 0.;
  }

  if (p_nuNR) {
    *p_nuNR = pNuNR;
  }

  //X   double rho_nuNR = nuNR() * massiveNeutrinoRho(0, a);
  return rho_nuNR;
}

//============== end of massive neutrinos ============================================


//X     //  const=7*pi**4/120.  = 5.68...
#define ZETA3 1.20205690
#define CONST 5.68219698

/*!
  Determine the mass of non-relativistic neutrinos.
  This function is called shortly before CmbCalc
  asks us to evaluate history()
*/
void Cosmos::getReady() {
  if (nuNR() == 0. || mass_nu_eV() == 0.) {
    amnu = 0.;
  } else {
    amnu = mass_nu_kbT();
  }
  //cout << "Set up neutrino mass: " << amnu << " - in eV: " << mass_nu_eV() << endl;
}

//============ CMBFAST ===========


double Cosmos::dtauda(double a) {   return 1.0/tau2adot(a2tau(a)); }


/*! Integrate the ionization fraction xe for hydrogen semi-implicitly 
   from tau to tau+dtau, treating tempb, a, and adot as constants 
   during this interval. 
   Ionization temperature and coefficient. 
   Two-photon decay rate (in 1/Mpc). 
   Switch for implicit (switch=1.0) or semi-implicit (switch=0.5) scheme. */

void Cosmos::ionize(double *tempb, double *a, double * adot, double *dtau, double *xe) {    
    static double beta, crec, bbxe, cpeebles, alpha, b1, alpha0, aa, bb, cp1, cp2, rat, phi2;

    //  Recombination coefficient (in sqrt(K)/Mpc). 
    alpha0 = (1 - Y_he()) * 2.3866e-2 * omega_b() * h2();

    //  Coefficient for correction of radiative decay (dimensionless). 
    crec = (1 - Y_he()) * 8.0138e-22 * omega_b() * h2();

    //  Recombination and ionization rates. 
    phi2 = log(157890. / *tempb) * .448;
    phi2 = max(phi2,0.);
    alpha = alpha0 / sqrt(*tempb) * phi2 / (*a * *a * *a);
    beta = *tempb * phi2 * exp(43.082 - 157890. / *tempb);
    //  Peebles' correction factor. 
    if (*tempb < 200.) cpeebles = 1.;
    else {
      cp1 = crec * 8.468e14 * (1. - *xe) / (*a * *adot);
      cp2 = crec * *tempb * phi2 * exp(43.082 - 39472.5 / *tempb) * (1. - * xe) / (*a * *adot);
      cpeebles = (cp1 + 1.) / (cp1 + 1. + cp2);
    }

    //  Integrate dxe=bb*(1-xe)-aa*xe*xe by averaging rhs at current tau 
    //  (fraction 1-switch) and future tau (fraction switch). 
    aa = *a * *dtau * alpha * cpeebles;
    bb = *a * *dtau * beta * cpeebles;
    b1 = bb * .5 + 1.;
    bbxe = bb + *xe - (bb * *xe + aa * *xe * *xe) * .5;
    rat = aa * .5 * bbxe / (b1 * b1);
    //  Prevent roundoff error. 
    if (rat < 1e-6) {
	*xe = bbxe / b1 * (1. - rat);
    } else {
	*xe = b1 / (aa * 1.) * (sqrt(rat * 4. + 1.) - 1.);
    }
} // ionize_ 

/*!
  Compute the helium ionization fractions using the Saha equation. 
  x0 is the hydrogen ionization fraction n(H+)/n(H) (input), 
  x1 is the helium first ionization fraction n(He+)/n(He) 
  (input and output), and 
  x2 is the helium second ionization fraction n(He++)/n(He) 
  (input and output).  */

void Cosmos::ionhe(double *tempb, double *a, double *x0, double *x1, double *x2) {
    static double b, c, x2new;
    static int niter;
    static double b0, r1, r2, xe, err;


    //  Ionization temperatures. 
    //
    //  Constant for electron partition function per baryon. 
    b0 = 2.15e20 / ((1. - Y_he()) * omega_b() * h2());

    //  Electron partition function per baryon. 
    b = b0 * *a * *a * *a * *tempb * sqrt(*tempb);

    //  Dimensionless right-hand sides in Saha equations. 
    if (fabs(285500. / *tempb) < 100.) r1 = b * 4. * exp(-285500. / *tempb); else r1 = 0.;
    if (fabs(631300. / *tempb) < 150.) r2 = b * exp(-631300. / *tempb); else r2 = 0.;

    //  Solve coupled equations iteratively. 
    c = Y_he() * .25 / (1. - Y_he());
    err = 1.;
    niter = 0;

    while (err > 1e-12) {
      ++niter;    
      xe = *x0 + c * (*x1 + *x2 * 2);
      x2new = r1 * r2 / (r1 * r2 + xe * r1 + xe * xe);
      *x1 = xe * r1 / (r1 * r2 + xe * r1 + xe * xe);
      err = fabs(x2new - *x2);           
      *x2 = x2new;
    }
} // ionhe_ 


/*! 
  This subroutine calculates the redshift of reionization 
  from an optdlss and sets the reionization fraction rif=1. */
void Cosmos::reiopar() {
  setReionizationZ(0);
  setReionizationFraction(0);
  if (optDistanceLss() != 0.) {
    double akthom = (1 - Y_he()) * 2.3048e-5 * omega_b() * h2();
    //*rif = 1.;
    setReionizationFraction(1);

    // Calculating the redshift of reionation the parameter
    // needed by CMBFAST.
    int na = 1;
    double da = 1e-5,a=0;
    double optd = 0.;

    while (optd < optDistanceLss()) {
      a = 1 - na * da;
      if ( a <= a_min()) {
        throw Bad_Error("Cosmos::reiopar - a of reionization <= a_min()");
      }
      optd += da * reionizationFraction() * akthom * dtauda(a) / (a * a);
      na++;
    }

    setReionizationZ(1.0/a - 1);
    setReionizationZStop(0);

  }
} //reiopar


/*!
  Return the fraction of free electrons as a function of redshift.
 
  If you would like to play with differen re-ionizations, this is definitly
  the place to play (but please make sure that determineTimeSteps has
  enough steps at your specific re-ionization redshifts, if you alter this routine here).

  For redshifts that are smaller than the highest one calculated in Recfast, it 
  simply returns this limiting value. 

  For re-ionization, it smoothly (but quickly) turns it on (just like cmbfast) and leaves
  it on from then on.
*/

double Cosmos::x_e(const double Z) {
  //  static double RiseInterval = 0.3;  // redshift interval before z_reionization in which re-ionization smoothly kicks in
  double x;
  double z = Z;
  if (z < 0.0) z = 0.0;  // make sure we don't get negative
  if (z < recombination->maxRedshift() ) {
    x = recombination->xe(-z);  // recombination stores (-z, xe) pairs
  } else { 
    x = recombination->highZFraction();   // else the value of the largest redshift stored in recombination
  }
  
  
  if (reionizationZ() != 0.0) {  // reionization ...
    // smooth tanh matching
    double tau = z2tau(z);
    double ReTau = z2tau(reionizationZ());
    double WouldBe = recombination->xe(-tau2z(0.9*ReTau)); // xe at 0.9 tau
    if (tau > 0.9*ReTau) {
      double leverage = 150*(tau-ReTau) / ReTau; 
      x =  0.5*(reionizationFraction() - WouldBe)*(1.0 + tanh(leverage)) + WouldBe;
    } 
  }
  return x;
}

/*!
  Given a wavenumber k, determine
  at which tau we should start to compute the evolution of perturbations.
*/
double Cosmos::tauStartforKMode(double k)
{
  //     Begin when wave is far outside horizon. 
  //     Conformal time (in Mpc) in the radiation era, for photons 
  //     plus 3 species
  //     of relativistic neutrinos. 
  double taustart = 1e-3 / k;
  //    Make sure to start early in the radiation era. We are not including 
  //    neutrinos as sources of the tensor modes. 
  taustart = min(taustart,1e-3*tau_equ());
  if (isnan(tau_equ()) || (tau_equ()<=0) || (taustart <= 0))
    taustart = tau_min();
  if (taustart <= tau_min()) {
    taustart = tau_min();
    cout << "Warning - using relatively late first tau (taumin= " << taustart << " for k=" << k << endl;
  }
//X   cout << "returning taustart " << taustart << "(taumin is "
//X        << tau_min() << ") for k=" << k << endl;
  return taustart;
}


/*!
  Given the maximum wavenumber in the game, determine
  at which tau the thermal history should start
*/
double Cosmos::determineStartOfThermo(double k) {
  return tauStartforKMode(k);
}

/*!
  The new thermal history routine. It first calls Recfast, which
  calculates the baryon temperature and the fraction of free 
  electrons per hydrogen atom: x_e

  Then, it logarithmically steps from tau_0 backwards in time,
  integrating the differential optical depth d / dtau Kappa.
  From this, it gets exp(kappa(tau)-kappa(tau(0)) and the 
  visibility function. It stores all of the necessary quantities
  (like derivatives of kappa and the visibility) in Splines,
  that the Perturbation classes will acess.
  
  After having done this, it calls determineTimeSteps(), which
  will (given the whole history calculated) determine at which
  point the anisotropy sources will be calculated
*/


void Cosmos::thermo(double MaxK, ControlPanel& control) {
  ThermoAnchor.kill();
  recombination = RecombinationFactory::recombination(*this, control, &ThermoAnchor);
  recombination->compute();
  //  recombination->noteX->dump("ionization");
  // sigma_thomson * number density of baryons today (times  a factor of a^3, which is 1 at intial time) 
  double ThomsonForCanonicalAlpha = (1 - Y_he()) * 2.3048e-5 * h2() * omega_b();   // Thomson cross-section for canonical fine structure constant

  double Kappa[2];
  double hnext = -1e-1;  // suggested initial stepsize
  Kappa[1] = 0;

  double tau = tau_0();  // we go from tau_0 to initial time
  double dtau = -10;
  double taumin = determineStartOfThermo(MaxK);

  double LogStepSpacing = (log(tau) - log(taumin)) / 2000; // two thousand steps
  do {
    if (tau < taumin) tau = taumin;
    double z = max(tau2z(tau), 0.0);
    double a = 1.0 /(z  + 1.0);

    dtau = - tau * LogStepSpacing;
    if (z > 500 && z < 2000) dtau *= 0.2;  // even more dense at likely redshifts of recombination
    if (reionizationZ() !=0 && tau > 0.9*z2tau(reionizationZ()) && tau < 1.1*z2tau(reionizationZ())) {
      dtau *= 0.1; // really dense in the rising and peak region
    } 

    // Here, we factor in a change in the fine structure constant modfying the Thomson cross
    // section. We do cut this off however in a smooth but manner for high redshifts
    double CutOffDelAlpha = recombination->delAlpha(z) * exp(- pow(z / 1e4,4.0)); // steep cutoff at z = 10000
    double Thomson = ThomsonForCanonicalAlpha*(1.0 + 2.0*CutOffDelAlpha);
    double cs2;
    double DotKappa =  Thomson*x_e(z)/(a*a);
    if (z < recombination->maxRedshift()) { // so recfast has something for us
      cs2 = recombination->cs2(-z);
    } else {
      double T_b =  T_cmb() / a; // Baryon temperature = gamma temperature;
      cs2 = 4.0/3.0 * recombination->BoltzmannOverWeightAndC2() *  T_b;  // see Ma & Bertschinger , eqn (68)
    }

    double ExpKappa = exp(Thomson*Kappa[1]);
    double visibility = DotKappa * ExpKappa;

    splineCs2->set(tau, cs2);
    Expmmu->set(ExpKappa);
    Opac->set(DotKappa);
    Visibility->set(visibility);
    FreeElectronFraction->set(x_e(z));
 
    hnext = Miscmath::odeint(Kappa, 1, tau, tau + dtau, 1e-6, hnext,0,
                                    (moDerivs)&Cosmos::integrateDotKappa,*this);
    bool reachedTauMin = false;
    if (tau==taumin) reachedTauMin = true;
    tau += dtau;
    //make sure we have a point for tau=taumin in the spline
    if ((!reachedTauMin) && (tau < taumin)) tau = taumin;
  } while (tau >= taumin);
  //  cout << "TAUMIN: " << taumin << endl;

  // flip the splines, cause we have gone from tau_0 -> small tau
  splineCs2->flip();
  Expmmu->flip();
  Opac->flip();
  Visibility->flip();
  FreeElectronFraction->flip();
  splineCs2->arm(Spline::thoseReady);  // arm those ready
  // if you like to see the splines, just uncomment these
  //   splineCs2->dump("soundspeed",false);
  //    Visibility->dump("visibility_new",false);
  //Opac->dump("opac",false);
  //Expmmu->dump("expmmu",false);
  // FreeElectronFraction->dump("ionization");

  // what is left to do is to calculated  d^2 kappa / dtau^2,
  // d visibility / dtau     and d^2 visibility / dtau^2
  // we do this here
  Anchor convenient;
  Spline *D3Kappa = new Spline(splineCs2, "DDDotKappa",&convenient); // only needed temporarily
  Opac->derive(*D2Kappa);
  D2Kappa->arm();
  D2Kappa->derive(*D3Kappa);
  D3Kappa->arm();
  for (int i = 0; i < Visibility->size(); i++) {
    double tau = Visibility->x(i);
    TimeStep[i] = tau;
    double dkappa = Opac->y(i);
    double d2kappa = D2Kappa->y(i);
    double d3kappa = D3Kappa->y(i);
    double expkappa = Expmmu->y(i);
    DVisibility->set( expkappa*(dkappa*dkappa + d2kappa ));
    DDVisibility->set( expkappa*(dkappa* (dkappa*dkappa + 3*d2kappa) + d3kappa));
  }
  DVisibility->arm();
  DDVisibility->arm();

  // finally, we can use this information to
  // find very definetely the special moments, as for instance last scattering
  ValidThermo = true;
  z__ls = tau2z(Visibility->maximum());
  findSpecialMoments();

  // and now, lo and behold, we calculate the timesteps
  determineTimeSteps();
}

/*!
  return d / dtau of kappa = a *n_e * sigma_thomson. 
  Actually, this returns x_e / (a*a), i.e. one has  to multiply
  this by a factor containing the baryon density and sigma_thomson.
  This is done after each step (no need to multiply in the constant before).
  The factor a^(-2) comes from rho_baryon = const * a^-3 and the factor
  of a from d / dtau kappa = a *n_e .... Together, two powers of a in the
  denominator
*/
void Cosmos::integrateDotKappa(const double tau, double *kappa, double *dotkappa) {
  double a = tau2a(tau);
  double z = 1.0/a - 1.0;
  dotkappa[1] =    x_e(z) / (a*a); // tau2adot(a2tau(a)); 
}

/*!
  Determines the tau-values at which the sources will be calculated.
  Naturally, it will focus on the epoch of recombination and reionization
  It's whole purpose is to fill the vector timeStep with the appropriate
  tau values.
*/
void Cosmos::determineTimeSteps() {  
  double tau_ls = Visibility->maximum();   // the peak of visibility is last scattering 
  double MaxVisibility = Visibility->fastY(tau_ls);  // maximum visibility
  
  // here, we note the crossing of 1e-4 and 1e-2 of the maximum visibility
  vector<double> LowLevel = Visibility->getZeroVector(1e-4*MaxVisibility, 1e-4);  
  vector<double> MidLevel = Visibility->getZeroVector(1e-2*MaxVisibility, 1e-4);

  double EarliestRelevantTime = LowLevel[0]; // before this, just forget the sources
  BeginRecombination = MidLevel[0]; 
  double EndRecomb = MidLevel[1];
  
  // here, we determine a confromal time interval in which the onset
  // of reionization happens
  double BeginReionization = tau_0();
  double EndReionization = tau_0(); 
  if (reionizationZ() != 0) {  
    BeginReionization = 0.95*z2tau(reionizationZ());
    EndReionization = 1.05*z2tau(reionizationZ());
  }
  double TauMax = Visibility->stop();  // this most definitly is tau_0

  double tau = EarliestRelevantTime;
  double dtau = 2; 

  
  //cout << "BeginRecombination: " << BeginRecombination << endl;
  // cout << "End of Recomb: " << EndRecomb << endl;

  SafeVector<double> tmp(10000);  // temporary storage place
  int i =0;
  do {  // in here, we determine the stepsize in tau
    if (tau > BeginRecombination && tau < EndRecomb) dtau = 1;  // during recombination this spacing
    if (tau > EndRecomb) dtau *= 1.05; // slowly increase spacing 
    if (tau < BeginReionization && tau+dtau > BeginReionization) dtau = max(BeginReionization - tau,1.0); // careful approach
    if (tau >= BeginReionization && tau < EndReionization) dtau = tau*1e-3;
    dtau = min(dtau, tau*1e-1); // at most, the logarithmic increase is limited to 1 %
    if (dtau > 50) dtau = 50; // level off at that stepsize

    tmp[i++] = tau; // store
    tau += dtau;    
  } while (tau < TauMax);

  if (tmp[i] < TauMax)  tmp[i++] = TauMax; // last step is today !!
  else throw Bad_Error("Cosmos::determineTimeSteps() tau >= TauMax should never happen");
  // cout << "TAUMAX: "<< TauMax << endl;
  // now copy the temporary array into the timeStep array
  TimeStep.resize(i);
  for (int k = 0; k < i;k ++) TimeStep[k] = tmp[k];

  //  cout << "NUMBER OF TIMESTEPS: " << i << endl;
}



/*! 
  Wrapper around finithermo() for cases where you're really only interested in the
  background. It is by no means faster but requests less parameters and is thus
  a nice and convenient way to get the last scattering redshift right.
void Cosmos::onlyThermo() {
  double taurend,dlntau0=0.01;
  int nstep=0,n1;
  finithermo(0.001,tau_0(), tau_0(), & taurend, &dlntau0,&n1,&nstep,true);
}
*/

  /*
  OBSOLETE, USE replacement thermo() instead !

  finithermo() certainly needs some attention. 

void Cosmos::finithermo(const double taumin, const double taumax, const double tau0, double *taurend, double *dlntau0, int * n1, int *nstep, bool peebles) {
  ThermoAnchor.kill();
  recombination = RecombinationFactory::recombination(Recombination::Recfast,*this,&ThermoAnchor);
 
#define NTHERMO 10000
     // Local variables 
    static double taurend1, taurend2, fact, adot;
    static int j2ri2;
    static double tau01, taurist1;
    static int j2ri3;
    static double dtau, adothalf, dtri, adot0, dtri0, a;
    static int i;
    static double ahalf;
   
    static double thomc, a0, a2;
    static int j1, j2;
    static double x1, x2, thomc0, a02, fe, tb[NTHERMO+1], dtbdla, xe[NTHERMO+1];
    static int iv;
    static double tbhalf;
    static int ns;
    static double barssc, akthom;
   
    static double dtauri;
    static double cf1;
  
    static int ncount;
    static double a2t, vismax, sdotmu[NTHERMO+1], xe0, tg0;
    static double dtauri0, etc, tgh, vfi, tau, xod, tautemp;

    // used to be members of common thermod, however nothing in common, so local
    double tauminn, dlntau; 


    Spline &Dotmu = *Opac;  
    Spline &D2opac = *new Spline(splineCs2,"d2otmu",&SplineAnchor);
    Spline &Emmu = *Expmmu; // nicer to the eye to acess

    // Compute and save unperturbed baryon temperature and ionization fraction
    // as a function of time.  With nthermo=10000, xe(tau) has a relative 
    //
    // accuracy (numerical integration precision) better than 1.e-5. 
    //
    //

    // If not peebles, call recfast

    if (!peebles) recombination->compute();  

    ncount = 0;
 
    thomc0 = pow(T_cmb(), 4) * 5.0577e-8;
    akthom = (1 - Y_he()) * 2.3048e-5 * omega_b() * h2();
    tauminn = taumin * .05;
    dlntau = log(tau0 / tauminn) / (NTHERMO-1);

    // Initial conditions nothing assumed :-) just ask cosmo's tau2xx() 
    tau01 = tauminn;

    a0 = tau2a(tauminn); // COM a0 = adotrad * tauminn;
    adot0 = tau2adot(tauminn);        // COM    adot0 = adotrad;           
    a02 = a0 * a0;

    // Assume that any entropy generation occurs before tauminn. 
    // This gives wrong temperature before pair annihilation, but 
    // the error is harmless. 
    tb[0] = T_cmb() / a0;
    xe0 = 1.;
    x1 = 0.;
    x2 = 1.;
    xe[0] = xe0 + Y_he() * .25 / (1. - Y_he()) * (x1 + x2 * 2);
    barssc = (1. - Y_he() * .75 + (1. - Y_he()) * xe[0]) * 9.182e-14;

    splineCs2->set(tau01,  barssc * 1.3333333333333333 * tb[0]);

    sdotmu[0] = 0.;
    Dotmu.set( xe[0] * akthom / a02 );     
   
    for (i = 2; i <= NTHERMO; ++i) {
      tau = tauminn * exp((i - 1) * dlntau);
      dtau = tau - tau01;
      
      // Integrate Friedmann equation using inverse trapezoidal rule. 
      a = tau2a(tau);                    //a0 + adot0 * dtau;
      a2 = a * a;
      adot = tau2adot(tau);  
      // Baryon temperature evolution: adiabatic except for Thomson cooling. 
      // Use  quadrature solution. 
      tg0 = T_cmb() / a0;
      ahalf = (a0 + a) * .5;
      adothalf = (adot0 + adot) * .5;

      // fe=number of free electrons divided by total number of free baryon 
      // particles (e+p+H+He).  Evaluate at timstep i-1 for convenience;if 
      // more accuracy is required (unlikely) then this can be iterated with 
      // the solution of the ionization equation.
      
      fe = (1. - Y_he()) * xe[i - 2] / (1. - Y_he() * .75 + (1. - Y_he()) * xe[i - 2]);
      thomc = thomc0 * fe / adothalf / (ahalf*ahalf*ahalf );  
      etc = exp(-thomc * (a - a0));
      a2t = a0 * a0 * (tb[i - 2] - tg0) * etc - T_cmb() / thomc * (1. - etc);
      tb[i - 1] = T_cmb() / a + a2t / (a * a);

      // Integrate ionization equation. 
      tbhalf = (tb[i - 2] + tb[i - 1]) * .5;
      // If there is re-ionization, smoothly increase xe to the 
      // requested value. 
      if (reionizationZ() != (float)0. && tau > reionizationTau() * 9. / 10.) {
	if (ncount == 0) ncount = i - 1;
	// SMOOTH REIONIZATION 
	xod = (tau - reionizationTau()) * 150. / reionizationTau();
	if (xod > 100.) tgh = 1.; else tgh = (exp(xod) - exp(-xod)) / (exp(xod) + exp(-xod));
	xe[i - 1] = (reionizationFraction() - xe[ncount - 1]) * (tgh + 1.) / 2. + xe[ncount - 1];
      } else {
	if (peebles) {
	  ionize(&tbhalf, &ahalf, &adothalf, &dtau, &xe0);
	  ionhe(&tb[i - 1], &a, &xe0, &x1, &x2);
	  xe[i - 1] = xe0 + Y_he() * .25 / (1. - Y_he()) * (x1 + x2 * 2);
	} else {
	  //cout << "calling for: " << a << "  "  << a2z(a) << endl;
	  double z = -a2z(a);  // recombination spline has stored -z, such that x-data is in proper order
	  if (z < recombination->noteX->start()) xe[i-1] = recombination->noteX->front(); else
	    xe[i-1] = recombination->xe(z); 
	}
      }

      // Baryon sound speed squared (over c**2). 
      // CMBEASY:: we deviate here from the original version of finithermo for the following reason:
      // The last term in the equation below, (a * tb[i - 1] - T_cmb()) gives the temperature
      // difference between baryons and gammas. Unfortunately, this difference is very small.
      // in fact, it is so small that in double accuracy, one usually gets 0 back. Except for very
      // few times, when due to some roundoff, something at the edge of double accuracy is 
      // returned, i.e.  approx 10^{-16} !!!! 
      // While we believe that the correct i.e. full equation should be used,
      // it is in practice practically never in effect (due to the aforementioned canceling due to
      // truncation effects). 
      // Yet,. the few times at early times, when the bracket is non-zero, it is multiplied
      // by an extremely large pre-factor. Thus the rounding noise is amplified, leading
      // to 2 orders of magnitude larger values than usual, and even worse: randomly 
      // distributed at negative and positive values !!!!
      // Having at the very moment no better handle on this, I therefore will
      // multiply the second part of this equation by zero. 

      double ANIHILATE_EXPRESSION_TO_PREVENT_NOISE_AMP=0;

      dtbdla = -2.0*tb[i - 1] - ANIHILATE_EXPRESSION_TO_PREVENT_NOISE_AMP * thomc * adothalf / adot * (a * tb[i - 1] - T_cmb());

      barssc = (1. - Y_he() * .75 + (1. - Y_he()) * xe[i - 1]) * 9.182e-14;

      splineCs2->set(tau, barssc * tb[i - 1] * (1 - dtbdla / tb[i - 1] / 3.) );		
    
      // Calculation of the visibility function 

      Dotmu.set(xe[i - 1] * akthom / a2);

      if (tau < 0.001) sdotmu[i - 1] = 0.;
      else {
	sdotmu[i - 1] = sdotmu[i - 2] + dtau * 2. / (1. / Dotmu.y(i-1) + 1. / Dotmu.y(i - 2));
      }
      //      cout << 1/a - 1.0 << "  Dotmu: " << Dotmu.y(i-2) << "  kapppa: " << sdotmu[i-1] << endl;

      a0 = a;
      tau01 = tau;
      adot0 = adot;
    }
    
   
    // prepare the spline data for later interpolation
    splineCs2->arm();
    // if you like to se cs2, uncomment it here. 
    // to see the noise effect, set ANIHILATE.... =1 
    // splineCs2->dump("baryonsound",false);
    Dotmu.arm();

    if (xe[NTHERMO-1] < reionizationFraction() && reionizationZ() != 0.) {
      cout << "Warning: We use a smooth function to" << endl;
      cout << "approach your specified reionization" << endl;
      cout << "fraction. The redshift that is deduced from" << endl;
      cout << "youre input paprameters is so low that our" << endl;
      cout << "smooth function does not reach the required" << endl;
      cout << "value. You should go in to subroutine finithermo" << endl;
      cout << "and play with the shape of this smooth function." << endl;
      cout << "Search for SMOOTH REIONIZATION for the place where" << endl;
      cout << "the function is set." << endl;
    }
    for (j1 = 0; j1 < NTHERMO; ++j1) {
      Emmu.set(exp(sdotmu[j1] - sdotmu[NTHERMO-1]) + 1e-30);
      Visibility->set(Emmu[j1]*Dotmu[j1]);
    }

    iv = 0;
    vfi = 0.;

    // Getting the starting and finishing times for decoupling. 
    if (ncount == 0) {
	cf1 = 1.;
	ns = NTHERMO;
    } else {
	cf1 = exp(sdotmu[NTHERMO-1] - sdotmu[ncount - 1]);
	ns = ncount;
    }
    for (j1 = 1; j1 <= ns; ++j1) {
      vfi += Emmu.y(j1-1) * Dotmu.y(j1-1) * cf1; 
	if (iv == 0 && vfi > 1e-5) {
	    taurst = tauminn * .90000000000000002 * exp((j1 - 1) * dlntau);
	    iv = 1;
	}
	if (iv == 1 && vfi > 0.99) {
	    taurend1 = tauminn * exp((j1 - 1) * dlntau) * 1.5;
	    taurend1 = max(taurend1,  taurend1 * sqrt(2500. / (omega_c() + omega_b()) / ( 1e4* h2())));
	    iv = 2;
	}
    }
    if (iv != 2) taurend1 = tauminn * exp((ncount - 1) * dlntau) * 1.5;
    
    // Calculating the timesteps during recombination. 
    if (dtaurec != 0.0)  dtaurec = min( dtaurec , taurst / 40. ); else dtaurec = taurst / 40.;
    
    taurend2 = dtaurec / *dlntau0;
    *taurend = max(taurend1,taurend2);
    *taurend = min(*taurend,reionizationTau() * 9. / 10.);
    // In models where reionization starts very early so that 
    // it cut into what we call recombination we will 
    // make the timesteps there at least as small as the 
    // ones we were using during recombination. If not 
    // timesteps 1.5 times bigger suffice. 

    fact = 1.5;
    
    if (*taurend == reionizationTau() * 9. / 10.) fact = 1.;
    *n1 = (int) rint((*taurend - taurst) / dtaurec + 1);
    dtaurec = (*taurend - taurst) / (*n1 - 1);
    ++(*n1);
    // Calculating the timesteps after recombination (logarithmic 
    // outside re-ionization scattering surface). 
    *nstep =  (int) rint(log(taumax / *taurend) / *dlntau0);
    *dlntau0 = log(taumax / *taurend) / *nstep;
    *nstep += *n1;

    // Adjusting if there is reionization. 
    // There will be nri0 points to sample the quick rise in 
    // the free electron density. After that, timesteps of length 
    // dtauri until tauristp. 
 
    reionization2.nri0 = 50;
    if (reionizationZ() != 0.) {

      // modified for lensing 
      // taurist1=10.0d0/9.0d0*taurend 
      taurist1 = reionizationTau();
      reionization2.j2ri1 = (int) (log(taurist1 * 9. / *taurend / 10.) / *dlntau0 + *n1);
      tautemp = min( taurist1 * 21. / 20., reionizationTauStop());
      
      j2ri2 = (int) (log(tautemp / *taurend) / *dlntau0 + *n1);
      j2ri3 = (int) (log(reionizationTauStop() / *taurend) / *dlntau0 + *n1);
      
      dtri0 = *taurend * (exp(*dlntau0 * (j2ri2 - *n1)) - exp(*dlntau0 * (reionization2.j2ri1 - *n1)));
      dtri = *taurend * (exp(*dlntau0 * (j2ri3 - *n1)) - exp(*dlntau0 * (j2ri2 - *n1)));
      dtauri0 = dtri0 / (double) reionization2.nri0;
      dtauri = dtaurec * fact;
      dtauri = min( dtauri ,dtri / 10.);
      
      if (dtauri > 0.) {
	reionization2.nri = (int) (dtri / dtauri) + 1;
	dtauri = dtri / (double) reionization2.nri;
      } else {
	reionization2.nri = 0;
	dtauri = 0.;
	}
      *nstep = *nstep + reionization2.nri0 + reionization2.nri + reionization2.j2ri1 - j2ri3;     
    } else {
      reionization2.j2ri1 = 0;
      j2ri2 = 0;
    }
  
    timeStep.resize(*nstep);

    Dotmu.derive(*D2Kappa);
    D2Kappa->arm();
    D2Kappa->derive(D2opac);
    D2opac.arm();
    Emmu.arm();
    
   

    for (int i =0; i < Dotmu.size(); i ++) {
      double v2 = Dotmu[i]*Dotmu[i]; 
      DVisibility->set(Emmu[i]*(v2 + (*D2Kappa)[i]));
      DDVisibility->set(Emmu[i]*(Dotmu[i]* (v2 + 3* (*D2Kappa)[i]) + D2opac[i]));
    }
    DVisibility->arm();
    DDVisibility->arm();

    // Saving the tau  steps 
    vismax = 0.0;
    for (j2 = 2; j2 <= *nstep; ++j2) {
      if (j2 <= *n1) {
	tau = taurst + (double) (j2 - 2) * dtaurec;
      } else {
	if (reionizationZ() == 0. || j2 <= reionization2.j2ri1) {
	  tau = *taurend * exp(*dlntau0 * (double) (j2 - *n1)); 
	} else {
	  if (j2 < reionization2.j2ri1 + reionization2.nri + reionization2.nri0) {
	    if (j2 <= reionization2.j2ri1 + reionization2.nri0) {
	      tau = timeStep[reionization2.j2ri1 - 1] + dtauri0 * (double) (j2 - reionization2.j2ri1);
	    } else {
	      tau = timeStep[reionization2.j2ri1 + reionization2.nri0 - 1] + dtauri * (double) (j2 - reionization2.j2ri1 - reionization2.nri0);
	    }
	  } else {
	    tau = *taurend * exp(*dlntau0 * (double) (j2 - reionization2.j2ri1 - reionization2.nri0 - reionization2.nri + j2ri3 - *n1));
	  }
	} 
      }
      timeStep[j2 - 1] = tau;
    }
    timeStep[0] = 0.;
    timeStep[*nstep-1] = min(tau0 , timeStep[*nstep-1] );
 
    // getting the peak visibility and hence the time of last scattering    
    Visibility->arm();
    ValidThermo = true;
    z__ls = tau2z(Visibility->maximum());
    findSpecialMoments();    
    //Visibility->dump("visibility",false);  // if you like to see it, here it comes
    // DVisibility->dump("dvisibility",false);
    //    Dotmu.dump("dotmu",false);
}
*/

// CMBEASY Stuff, once again 

double Cosmos::powerSpectrum(const double k, double n, const double dnsdlnk) {     
  //return  5.0* pow(k*20.0, n-1 + 0.5*dnsdlnk*log(k*20.0)); 
  return   pow(k*20.0, n-1 + 0.5*dnsdlnk*log(k*20.0)); 
}

double Cosmos::scalarSpectrum(const double k, unsigned int q) {
  return powerSpectrum(k, InitialPower[q], InitialPower_dnsdlnk[q]);
}


double Cosmos::tensorSpectrum(const double k, int n, const ControlPanel& control) {
  return pow(k*20,InitialTensorPower[n]);
}


/*!
  Just a small loop for initializing InitialPower[] array with
  'steps' values from min to max.

  Also, it initializes the PowerNormalization[] factor for the cdm spectrum
  to unity.
*/
void Cosmos::setInitialPower(double min, double max,int steps) {
  InitialPower.clear();
  InitialPower_dnsdlnk.clear();
  PowerNormalization.clear();

  if (steps > 1) {
    for (int i=0; i < steps; i++) {
      InitialPower[i] = min + i*(max-min)/(steps-1);
      InitialPower_dnsdlnk[i] = 0.0;
      PowerNormalization[i] = 1.;
    }
  } else {
    InitialPower[0] = min;
    InitialPower_dnsdlnk[0] = 0.0;
    PowerNormalization[0] = 1.;
  }
}


void Cosmos::cobeNormalize(const ControlPanel& control,CL& cl, const vector<int>& lval) {
  cout << "cobeNormalize() is obsolete. If you would like to WMAP normalize" << endl;
  cout << "use the capabilities of AnalyzeThis." << endl;
  throw Bad_Error("Cosmos::cobeNormalize() support for function discontinued");

  /*
  static double delt, alnk, dlnk, xlnh, curv, d1ppr;
  static double r, s, x, xlog10, d2norm,  hc, ak;
  static double xl[200], sx, sy;
  
  static double  win, sxx, sxy;
  
  static double apowers, d1pr;
  
  double c10 = 1;

  sigma8.resize(InitialPower.size());
    
    // Function Body 
    xlog10 = log(10.);
    xlnh = log(h());

    //cout <<"output: " << TensorRatio[0] << endl;

    for (unsigned  int j = 1; j <= lval.size() ; ++j) xl[j]  = lval[j - 1];       // shortcut...
    
    // Curvature radius 
    if (abs(omega_k()) > .001) {
      hc = 2.998e1 / h();
      curv = -omega_k() / hc / hc;
      r = 1. / sqrt((abs(curv)));
    }
   
    // COBE normalization 
    // fit the spectrum to a quadratic around C_10 with equal weights n logl 
    for (unsigned int n = 0; n < InitialPower.size(); ++n) {

      // first, lets normalize to 1 (this has previously been done in cmbflat.f)      
      double cl2 = cl.ts[n]->y(0);  
      double normalize = cl2;  
      cout << "cl2: "<< cl2 << endl;
      if (normalize != 0.0) normalize = 1.0/normalize ;
      *cl.ts[n] *= normalize;      // use scalar product definition of spline with number
      *cl.es[n] *= normalize;     // from spline.h . This just multiplies all y-data with cl...
      *cl.cs[n] *= normalize;
      *cl.kk[n] *= normalize;
      *cl.tk[n] *= normalize;
      // and now the same for tensors...
      // if original cmbfast behaviour is requested (or a free chosen T/S ratio is requested),
      // or no scalars have been calculated
      // we divide by tensor quadrupole,
      // otherwise we divide by scalar quadrupole
      if ( control.originalTensorHandling || (! control.canonicalTensorRatio)  ||  (!control.scalar) ) normalize = cl.tt[n]->y(0);
      else normalize = cl2;

      if (normalize != 0.0) normalize = 1.0/normalize ;
      *cl.tt[n] *= normalize;
      *cl.et[n] *= normalize;
      *cl.bt[n] *= normalize;
      *cl.ct[n] *= normalize;

      // Ensuring scalar to tensor ratio (if cmbfast's original version is requested)
      if (control.scalar && control.tensor && control.originalTensorHandling ) {
      	(*cl.tt[n]) *= TensorRatio[n];   // skalar multiplication of spline y-data with number
	(*cl.et[n]) *= TensorRatio[n];
	(*cl.bt[n]) *= TensorRatio[n];
	(*cl.ct[n]) *= TensorRatio[n];
	(*cl.kk[n]) *= TensorRatio[n];
	(*cl.tk[n]) *= TensorRatio[n];
	cout << "I DID MULTIPLY WITH THE T/S : " << TensorRatio[n] << endl;

      }


      // we need clts's and cltt's interpolation capabilities, so we have to arm()
      // later, we will disarm() in order to change the y-data again (normalizing)
      cl.ts[n]->arm();       // we can ALWAYS arm, because driver has called setN...
      cl.tt[n]->arm();      


      c10 = (*cl.ts[n])(10) + (*cl.tt[n])(10);

      //cout << "clts[10] : " << (*cl.ts[n])(10) <<  "   " <<  (*cl.tt[n])(10) << endl; 
      
      double d1=((*cl.ts[n])(xl[2]) + (*cl.tt[n])(xl[2]))/c10-1.0;
      double d2=((*cl.ts[n])(xl[3]) + (*cl.tt[n])(xl[3]))/c10-1.0;
      double d3=((*cl.ts[n])(xl[5]) + (*cl.tt[n])(xl[5]))/c10-1.0;
      double d4=((*cl.ts[n])(xl[7]) + (*cl.tt[n])(xl[7]))/c10-1.0;
      double d5=((*cl.ts[n])(xl[10]) + (*cl.tt[n])(xl[10]))/c10-1.0;
      double d6=((*cl.ts[n])(xl[11]) + (*cl.tt[n])(xl[11]))/c10-1.0;
      double d7=((*cl.ts[n])(xl[12]) + (*cl.tt[n])(xl[12]))/c10-1.0;

      double x1 = log(xl[2]) / xlog10 - 1.;
      double x2 = log(xl[3]) / xlog10 - 1.;
      double x3 = log(xl[5]) / xlog10 - 1.;
      double x4 = log(xl[7]) / xlog10 - 1.;
      double x5 = log(xl[10]) / xlog10 - 1.;
      double x6 = log(xl[11]) / xlog10 - 1.;
      double  x7 = log(xl[12]) / xlog10 - 1.;

      sy = x1 * d1 + x2 * d2 + x3 * d3 + x4 * d4 + x5 * d5 + x6 * d6 + x7 * d7;
      s = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4 + x5 * x5 + x6 * x6 + x7 * x7;
      sx=x1*x1*x1 + x2*x2*x2 + x3*x3*x3 + x4*x4*x4 + x5*x5*x5 + x6*x6*x6 + x7*x7*x7;
      sxy=x1*x1*d1+x2*x2*d2+x3*x3*d3+x4*x4*d4+ x5*x5*d5+x6*x6*d6+x7*x7*d7;
      sxx=x1*x1*x1*x1 +x2*x2*x2*x2 + x3*x3*x3*x3 + x4*x4*x4*x4 + x5*x5*x5*x5 + x6*x6*x6*x6 + x7*x7*x7*x7;
      delt = s * sxx - sx * sx;
      d1pr = (sxx * sy - sx * sxy) / delt;
      d1ppr = (s * sxy - sx * sy) * 2. / delt;
      // Bunn and White fitting formula 
      c10 = (d1pr * .02282 + .64575 + d1pr * .01391 * d1pr - d1ppr * .01819 
	     - d1pr * .00646 * d1ppr + d1ppr * .00103 * d1ppr) / c10;
  

      // density power spectrum normalization; 
      // Look at delta.f for an example of how to use it. 
      if (control.power_cdm && control.cmb) { 
	// get the power spectrum (still unnormalized, however with the
	// intiial scalarSpectrum()  already multiplied in)
	Spline* cdm = createPower(n,"COBENORMALIZE_cdm", power_cdm() ); 
	d2norm = c10 * 1.1e-9 / cl2;  
	cout << "d2norm: " << d2norm << endl;

	// this is to be multiplied by (k/h)^{n+3} to get 4pik^3P(k)  
	
	double dsig8 = 0.;
	double dsig8o = 0.;
	double sig8 = 0.;
	double alnko = 0.;
	double ict0 = 0;
	
	for (int i = 0; i < cdm->size(); i++) {	  
	  ak = cdm->x(i)*h();  // power has been stored as k/h now we make k = k/h *h :-)
	  x = cdm->x(i) * 8.;
	  win = (sin(x) - x * cos(x)) * 3. / (x * x * x);
	  alnk = log(ak);
	  dlnk = alnk - alnko;
	  if (ict0 == 0) {
	    dlnk = 2. / (InitialPower[n] + 3.);
	    ict0 = 1;
	  }
	  if (fabs(omega_k()) < .001) {
	    apowers = powerSpectrum(ak, InitialPower[n]);
	  } else {
	    cout << "omega_k() :  " << omega_k() << endl;
	    throw Bad_Error("Cobenormalize() This is an open model, not supported :-)");
	    // powersopen_(&ak, &n, &apowers);
	    apowers /= ak;
	  }
	   
	  // tfc as defined by CMBFAST is deltac / k^2
	  //
	  // cdm holds (deltac / k^2)^2 * powerSpectrum(), hence this is tfc^2*apowers
	  
	  dsig8 = win * 12.566370600000001 * win * exp(alnk * 5.) * cdm->y(i)/(ak*ak);
	  sig8 += (dsig8 + dsig8o) * dlnk / 2.;
	  dsig8o = dsig8;
	  alnko = alnk;
	}
	
	sig8 = sqrt(sig8 * d2norm);
	cout << "sigma8 = " << sig8 << endl;
	sigma8[n] = sig8;
	//cout << "4pik^3P(k) norm: " << d2norm << endl;
	
	//       Change d2norm to be consistent with lensing, 
	//       where we use k/h 
	// 	     d2norm=d2norm*exp(xlnh*3.0d0) 

	// in principle, this should not happen, i.e. i will have
	// to move the line multiplying the cdm spectrum
	// above the 	d2norm *= exp(xlnh * 4.);     statement 
	// conventionally, we have h^3, not h^4


	d2norm *= exp(xlnh * 4.);       // this is just d2norm *= pow(h,4);

	// please do also note that this definition of P(k)
	// additionally differs from the the usual one 
	// by the fact that here, sigma8 = int 4*pi ...
	// and usually sigma8 = int 1/(2*pi*pi)
	// i.e. one has to multiply by 8*pi^3 to compare this 
	// We do this here:
	PowerNormalization[n]  = 8*pow(M_PI,3)*d2norm/h();	// this is the correct normalization
	// Please note that PowerNormalization is therefore proportional to h^3
	// For sigma8() determination in AnalyzeIt (and currently also somewhere in the gui, I believe)
	// this h^3 is divided out again :-). 
	delete cdm;

      }
      // C_l normalization; output l(l+1)C_l/twopi 
      c10 = c10 * 2.2e-9 / 12.566370600000001;

      cout << "c10 factor: " << c10 << " and with c2: " << c10/cl2 << endl;
      
      cl.ts[n]->disarm();
      cl.tt[n]->disarm();

      (*cl.ts[n]) *= c10;
      (*cl.es[n]) *= c10;
      (*cl.cs[n]) *= c10;
      (*cl.tt[n]) *= c10;
      (*cl.et[n]) *= c10;
      (*cl.bt[n]) *= c10;
      (*cl.ct[n]) *= c10;
      (*cl.kk[n]) *= c10;
      (*cl.tk[n]) *= c10;
 
    } 
*/
} // cobenormalize_ 


/*!
  Set a new point in (k/h, tau) space in each power_cdm() web known to cosmos.
  It uses the InitialPower[] slopes and you only have to provide k, tau and the
  'raw' density contrast (\delta_cdm) (i.e in cmbflat: y[2]). fillPower() will divide this
  by k^{-3} 

  Please note that from v4.0 on, this is k^{-3} as opposed to k^{-4}, because of
  the additional factor of k in  scalarPower() from v4.0
  
  In addition, we now store the logarithm and not the argument itself. Makes 
  interpolation much  more robust
*/
void Cosmos::fillPower(SplineWeb* power,double k, double tau, double raw) {
  //  power->set(k/h(), tau, raw*raw/(k*k*k*k)); 
  power->set(k/h(), tau, log(raw*raw/(k*k*k))); 
}

void Cosmos::fillPsiWeb(SplineWeb* power,double k, double tau, double raw) {
  power->set(k, tau, raw);
}

/*! \brief obsolete - moved to CmbCalc::fillPower()
  High-Level routine to fill power-web at given conformal time tau by
  asking Perturbation *pert;

  The points are stored in  (k/h, tau) web
*/
void Cosmos::fillPower(Perturbation *pert, double tau) {
  throw Bad_Error("Cosmos::fillPower  - this routine is obsolote; the functionality moved to cmbcalc.cc");
  double k = pert->k;
  fillPower(power_cdm(),k,tau,pert->delta_c());
  fillPower(power_baryon(),k,tau,pert->delta_b());
  fillPower(power_gamma(),k,tau,pert->delta_g());
  fillPower(power_nu(),k,tau,pert->delta_n());
  fillPsiWeb(power_psi(), k, tau, pert->psi());
  if (nuNR() != 0) {
    fillPower(power_nuNR(),k,tau,pert->delta_nr());
  }
}

#warning move LogStepper somewhere more appropriate / or use NewLogStepper
class LogStepper
{
  public:
    LogStepper( const double start, const double end): mLogStart(start), mLogEnd(log(end)), mLogX(log(start)) {}
    void setSteps(const int steps) { mLogStep = (mLogEnd-mLogStart)/(double)steps; }
    double next() {
      const double arg = mLogX;
      mLogX += mLogStep;
      if (mLogX > mLogEnd)
        mLogX = mLogEnd;
      return exp(arg);
    }
  private:
    double mLogStart, mLogEnd, mLogX, mLogStep;
};


Spline* Cosmos::growthFactor(const double k, std::string name, GrowthFactorType species, double* norm, Anchor* a)
{

  static const int points = 5555;
  const double last1pZ = 1205;
	
	// Get the time in PMmodels Units - H_0^1
    double H_t = 9.776/h();

    //cout << " H_t(): " << H_t << endl;
  Spline* retSplineA = new Spline(points, "growthFactorA", a);
  Spline* retSplineAZ = new Spline(points, "growthFactor_over_a_Z", a);
  Spline* retSplineZ = new Spline(points, "growthFactorZ", a);
  Spline* retSplineT = new Spline(points, "growthFactorT", a);

  Spline* retSplineTDot = new Spline(points, "growthFactorTDot",a);
  Spline* retSplineZDot = new Spline(points, "growthFactorZDot",a);
  Spline* retSplineADot = new Spline(points, "growthFactorADot",a);

  Spline* retSpline_Dot = new Spline(points, "growthFactor_Dot",a);

  LogStepper stepper(1, last1pZ);
  stepper.setSteps(points);
  double onepz = stepper.next();

  SplineWeb* powerSpline = 0;
  if (species == Cosmos::Cdm) {
    powerSpline = power_cdm();
  } else if (species == Cosmos::Baryons) {
    powerSpline = power_baryon();
  } else if (species == Cosmos::Weighted) {
    powerSpline = power_cdm();
  } else {
    throw Bad_Error("Cosmos::growthFactor() - unknown species.");
  }

    Spline* power_norm = createPower(0, "deltaDPower_norm", powerSpline, 0, z2tau(0));
	power_norm->arm();
 /*   
  *   double normD_0;
  *   normD_0 = sqrt(Pk2Delta2(power_norm, k));
  *   cout << " Power_norm D0: " << normD_0 << endl;
  */

  while (onepz < last1pZ) {

    Spline* power = createPower(0, "deltaDPower", powerSpline, 0,  z2tau(onepz-1.));
       	  if (power->size() == 0)
      throw Bad_Error("Cosmos::growthFactor():\n No  perturbations saved - forgot to compute them?\n");
    power->arm();

    double delta = sqrt(power->fastY(k)/power_norm->fastY(k));

    if (species == Cosmos::Weighted) {
      const double tau = z2tau(onepz-1.);
      Spline* power_b = createPower(0, "deltaDPower", power_baryon(), 0,  tau);
      power_b->arm();
      const double delta_b = sqrt(Pk2Delta2(power_b, k));
      delete power_b;

      const double omega_m = tau2rho_m(tau)/tau2rho(tau);
      const double omega_c = tau2rho_cdm(tau)/tau2rho(tau);
      const double omega_b = tau2rho_b(tau)/tau2rho(tau);
      delta *= omega_c/omega_m;
      delta += delta_b * omega_b/omega_m;
    }

    delete power;

    double time;
    const double a = z2a(onepz-1.);
  //  double gF = delta/a;
	double gF = delta;

	// Time is given in GYrs
    time = mpc2year()*tau2t(z2tau(onepz-1))*1.e-09;
	
    /* Fill up the three splines wrt a, z, t */
    retSplineA->set(a, gF);
    retSplineZ->set(onepz -1, gF);
    retSplineAZ->set(onepz -1, gF/a);
    retSplineT->set(time, gF);
    //cout << " z: " << a << " z: " << onepz-1 << " T: " << time << endl; 
    onepz = stepper.next();
  }
//	cout << "time(z=60): " << mpc2year()*tau2t(z2tau(60))*1.e-09 << " GYrs" << endl;
	
   	
		// Normalization is the same for all splines
	double normalize;  

	/* At zNorm normalize the GF*/
	double zNorm=20;

	retSplineZ->arm();
	normalize = (1./(zNorm + 1.))/retSplineZ->fastY(zNorm);

	/* Use this when using normalized splines*/
	//normalize = 1./retSplineZ->fastY(zNorm);
	//normalize = 0.788 ;
	// cout << " normalize^-1: " << retSpline->fastY(29.3030) << endl;
	
	*retSplineA *= normalize;
	*retSplineZ *= normalize;
	*retSplineT *= normalize;

	//cout << " ++++++++++++++ normalize: " << normalize << endl;

	retSplineT->flip();
	retSplineT->arm();
	retSplineT->dump("growthFactor_T");

	retSplineT->derive(*retSplineTDot);
	retSplineTDot->arm();
	retSplineTDot->dump("growthFactor_TDot");
       
	retSplineZ->dump("growthFactor_Z");
/*
	retSplineZ->derive(*retSplineZDot);
  	retSplineZDot->arm();
	retSplineZDot->dump("growthFactor_ZDot");
*/
	retSplineA->flip();
	retSplineA->arm();
	retSplineA->dump("growthFactor_A");
	retSplineA->derive(*retSplineADot);
	retSplineADot->arm();
	retSplineADot->dump("growthFactor_ADot");
	retSplineAZ->arm();
	retSplineAZ->dump("growthFactor_over_a_Z");

	/* get the time derivative using g(a) ---> (dg / da) * (da / dt)*/
  const double amax = 1;
  double aa=0.001; 
  //points = 5000;
  //double aa=1.;
  //double amax = 1000.;
	
	double step = (amax-aa)/points;
	
	while (aa<amax) {
	//double gdot = retSplineZDot->fastY(1./aa -1)*(Tau2ADot->fastY(aa))*(-1./(aa*aa));
	//double gdot = retSplineADot->fastY(aa)*(Tau2ADot->fastY(aa))*aa; //*(-1./(aa*aa));

	double zz = 1./aa - 1. ;
      double tt   = mpc2year()*tau2t(z2tau(zz))*1.e-09;
	     
      double gdot = retSplineTDot->fastY(tt);

//      cout << " tt: "  << tt << " gdot: " << gdot << endl; 

	/* Multiply gdot by H_t to get the correct ICgenerator Units*/
	retSpline_Dot->set(zz,gdot*H_t);
	step=((amax-step)/points);
	aa=exp(step+log(aa));
	
	//CnInvariant::printOSF5("LCDM_cosmicvector", tt, aa, retSplineT->fastY(tt), gdot, tau2Hubble(z2tau(zz))/H_0_cpm());	
	}


	retSpline_Dot->flip();
	retSpline_Dot->arm();
	retSpline_Dot->dump("growthFactor_dGdt_vs_Z");
  
	return retSplineT;
}


void Cosmos::dumpGrowthFactor(double k){
	Spline *s = new Spline(300, "gf", &SplineAnchor);
	Spline *sDot = new Spline(300, "gf_dot", &SplineAnchor); 

 // Spline* retSplineA = new Spline(points, "growthFactorA", a);
	std::string nam;
	nam="LCDM_growth_factor";
	Cosmos::GrowthFactorType f;
	f = Cosmos::Cdm;
	//f = Cosmos::Weighted;
	double norm;
	norm = 0.0915;
	s = growthFactor(k, nam, f, &norm, &SplineAnchor);
	//s->flip();
	//s->arm();
	//s->dump("growth_factor_t");
	//s->derive(*sDot);
	//sDot->arm();
	//sDot->dump("growth_factor_dot_t");
}





/*
Spline* Cosmos::growthFactor(const double k, std::string name, GrowthFactorType species, double* norm, Anchor* a)
{
  double normFactor = -1;

  if (norm) { // null pointer -> don't care about the used factor,
              //  value 0 -> put used factor in *norm, value != 0 ->use specific value
    if (*norm != 0) {
      normFactor = *norm;
    }
  }

  static const int points = 600;
  //const double last1pZ = 1205;
  const double last1pZ = 101;

  Spline* retSpline = new Spline(points, name, a);
  LogStepper stepper(1, last1pZ);
  stepper.setSteps(points);
  double onepz = stepper.next();

  SplineWeb* powerSpline = 0;
  if (species == Cosmos::Cdm) {
    powerSpline = power_cdm();
  } else if (species == Cosmos::Baryons) {
    powerSpline = power_baryon();
  } else if (species == Cosmos::Weighted) {
    powerSpline = power_cdm();
  } else {
    throw Bad_Error("Cosmos::growthFactor() - unknown species.");
  }
  while (onepz < last1pZ) {
    Spline* power = createPower(0, "deltaDPower", powerSpline, 0,  z2tau(onepz-1.));
    if (power->size() == 0)
      throw Bad_Error("Cosmos::growthFactor():\n No  perturbations saved - forgot to compute them?\n");
    power->arm();
    double delta = sqrt(Pk2Delta2(power, k));

    if (species == Cosmos::Weighted) {
      const double tau = z2tau(onepz-1.);
      Spline* power_b = createPower(0, "deltaDPower", power_baryon(), 0,  tau);
      power_b->arm();
      const double delta_b = sqrt(Pk2Delta2(power_b, k));
      delete power_b;

      const double omega_m = tau2rho_m(tau)/tau2rho(tau);
      const double omega_c = tau2rho_cdm(tau)/tau2rho(tau);
      const double omega_b = tau2rho_b(tau)/tau2rho(tau);
      delta *= omega_c/omega_m;
      delta += delta_b * omega_b/omega_m;
    }

    const double a = z2a(onepz-1.);
    //double gF = delta/a;
    double gF = delta;
    delete power;

    // at z=0
    if ( normFactor == -1 && a == 1.) {
      normFactor = gF;
      if (norm) {
        *norm = normFactor;
      }
    }
    if (normFactor == -1.)
      throw Bad_Error("Cosmos::growthFactor():\n\texpected the normalization to be set.");


    gF /= normFactor;
    retSpline->set(onepz, gF);

    onepz = stepper.next();
  }

  return retSpline;
}


void Cosmos::dumpGrowthFactor(double k){
	double norm=0.76; // LCDM value at z=0
	double a,t,z,tau,gf;
	double unit=9.776/h();
	int size=600;
	std::string name1="growthFactor2Z";
	std::string name2="growthFactorOverA2Z";
	std::string name3="growthFactor_T_2t";
	std::string name4="growthFactor_dt_2t";
	std::string name5="growthFactor_dt_2Z";

	Cosmos::GrowthFactorType f;
	f = Cosmos::Cdm;
	growth2Z      = new Spline(size, name1);
	growthOverA2Z = new Spline(size, name2);
	growthT2t  = new Spline(size, name3);
	growthDt2t = new Spline(size, name4);
	growthDt2Z = new Spline(size, name5);

	growth2Z= growthFactor(k, name1, f, &norm, &SplineAnchor);
	growth2Z->arm();
	growth2Z->dump();

	for(int k=0; k<size; k++) {
	gf=growth2Z->front(k);
	z=growth2Z->start(k)-1;
	a=1./(1.+z);
	tau=a2tau(a); 
	t=a2t(a)*mpc2year()*1.e-09;
	growthOverA2Z->set(z,gf/a);
	growthT2t->set(t,gf);
}
	growthT2t->flip();
	growthT2t->arm();
	growthT2t->derive(*growthDt2t);
	growthDt2t->arm();

	for(int k=0; k<size; k++) {
	z=growth2Z->start(k)-1;
	a=1./(1.+z);
	t=a2t(a)*mpc2year()*1.e-09;
	growthDt2Z->set(z,growthDt2t->fastY(t));
}

	growthDt2Z->arm();
	growthDt2Z->dump();
}
 */




/*!
  Saves power spectrum to disk. If SplineWeb *power is 0, then cdm_power
  is used as default
*/
void Cosmos::dumpPower(int n, string name, SplineWeb *power, double z) {
cout << "Cosmos::dumpPower() " << endl;
 if (! power) power=power_cdm(); 
double tau = z2tau(z); 
  //Spline *s = createPower(n,name,power,0,tau);
  Spline *s = allMatterPowerSpline(0);
  s->dump(name);
  delete s;
}

Spline* Cosmos::growthIndex(Spline *gf, Anchor *a, string name){
gf = growthFactor();
gf->arm();
Spline *gamma = new Spline(1000,name, a);
Spline *ln_gf = new Spline(1000,"ln_gf", a);
Spline *d_ln_gf = new Spline(1000,"d_ln_gf", a);
double a_min = 0.01; double a_max = 1.; 
double step = -log(a_min/a_max)/(1000); 
for(int i=0; i<1000; i++){
double a = a_min*exp(i*step);
double lngf = log(a*gf->fastY(1./a)/gf->fastY(100));
//cout << "a: " << a << " ln_gf: " << lngf << endl; 
ln_gf->set(log(a), lngf);
}
ln_gf->arm();
ln_gf->derive(*d_ln_gf);
d_ln_gf->arm();
for(int i=0; i<1000; i++){
double a = a_min*exp(i*step);
double a2t = a2tau(a);
double om = (tau2rho_cdm(a2t))/tau2rho(a2t);
double g = 1./log(om) * log(d_ln_gf->fastY(log(a)));
//cout << "a: " << a << " om: " << om << " d_ln_gf: " << d_ln_gf->fastY(log(a)) << endl; 
//cout << "a: " << a << " gf: " << a*gf->fastY(1./a)/gf->fastY(100) << endl;
gamma->set(a,g);
}
return gamma;
}


Spline *Cosmos::allMatterPowerSpline(double z){
double k, step_k, pk_tot, kmax, kmin; int size_k=1000;
double tau = z2tau(z); 
  Spline *s1 = createPower(0,"cdm",power_cdm(),0,tau);
  Spline *s2 = createPower(0,"baryons",power_baryon(),0,tau);
  Spline *s = new Spline(size_k, "all_power_spline", &SplineAnchor);
// Create spline containing the whole matter power spectrum
kmax=s1->stop(); kmin=0.001;
cout << " k_max: " << kmax << " k_min: " << kmin << endl;
step_k=log(kmax/kmin)/(size_k-1);
s1->arm(); s2->arm();
for(int i=0; i<size_k; i++){
k=kmin*exp(step_k*i);

if(pkWithBaryons()) {  
   pk_tot = s1->fastY(k) + s2->fastY(k);
} else {
pk_tot = s1->fastY(k);
}
 
//cout << " k: " << k << " p_k: " << pk_tot << endl;
	s->set(k,pk_tot);
}
return s;
}

MassFunction* Cosmos::dumpMassFunction(string name, double z, bool print, double norm) {
if(z==0) cout << "Cosmos::dumpMassFunction() " << endl;
// The mass function needs density at z=0 
	double tau0 = z2tau(0);

	stringstream red;	
	red << z;
	
        Spline *s = allMatterPowerSpline(z);
	MassFunction *mf;
	string tinker ("tinker");
	string tinker_z ("tinker_z");
	string sheth_tormen ("sheth_tormen");
	string sheth_tormen_z ("sheth_tormen_z");

	double alpha = -1.0;
	double Mmin = 1.e+9;
	double Mmax = 1.e+16*pow(1+z,alpha);
	//double Mmax = 1.e+16;
	double rho_m;

if(pkWithBaryons()){
	rho_m = (tau2rho_cdm(tau0) + tau2rho_b(tau0))*cpm2msun();
} else {	
	rho_m = (tau2rho_cdm(tau0))*cpm2msun();
}

if (name==tinker) {
	mf = new Tinker(s, Mmin, Mmax, rho_m);
} else if (name==tinker_z) {
	mf = new TinkerZ(s, Mmin, Mmax, rho_m);
} else if (name==sheth_tormen) {
	mf = new ShethTormen(s, Mmin, Mmax, rho_m);
} else {
 cout << "dumpMassFunction(), MassFunction:  " << name  << " not found." << endl;
}

	mf->set_z(z);
	mf->init(norm);
//name = name + red;
if(print==true){
	mf->print_differential(name);
	mf->print_integral(name);
}

  delete s;
  return mf;
}

double Cosmos::Sigma8_z(double z){
MassFunction *mf = dumpMassFunction("tinker",z, true, 1.);
return mf->sigma8();
}

Spline *Cosmos::dumpNumberDensity(string name, double Mass) {
cout << "Cosmos::dumpNumberDensity() " << endl;

// Calculation settings
double zp1, step_z, zmax, zmin; int size_z = 55;
double rho_m, time, time0;
zmax=3; zmin=0;

std::ostringstream oss;
oss << Mass;
std::string value = oss.str();

Spline *nd = new Spline(size_z, "number_density_M_z_"+name+"_M_"+value, &SplineAnchor);
step_z=(zmax-zmin)/(size_z-1);
zp1=zmin;

MassFunction *mf_0 = dumpMassFunction(name,0, true, 1.);

double norm = 0.8/mf_0->sigma8();
//double norm=1.;

cout << "Sigma8 file: " << mf_0->sigma8() << " normalization: " << norm << endl;

for (int j=0; j<size_z; j++) {
MassFunction *mf = dumpMassFunction(name, zp1, false, norm);
	nd->set(zp1, mf->integral_number_density(Mass));
zp1+=step_z;
//cout << j << " z: " << zp1 << " mf  " << mf->integral_number_density(Mass) << endl; //<< " pk(0.1): " << s->fastY(0.1) << " rho_m: " << rho_m << " tau: " << time << endl;
  delete mf;
} // end z-for

nd->dump();
nd->arm();
return nd;
}

double Cosmos::integrateOverComovingVolume(string name, double M, double z0, double z1, double the1, double the2, double phi1, double phi2){
double integ1, integ2; double angle;
double z, z_step; int size = 100;
Miscmath *mm;
angle = mm->solid_angle(the1,the2,phi1,phi2);
Spline *integral = new Spline(size, "integral_over_comvol", &SplineAnchor);
Spline *nd = dumpNumberDensity(name, M);

//z_step = (z1-z0)/(size-1);
//z=z0;
z_step = (z1-0)/(size-1);
z=0;

for(int i=0; i<size; i++){
double y = nd->fastY(z)*comovingVolume(0,z);
integral->set(z,y);
//cout << "z: " << z << " y: " << y << endl;
z+=z_step;
}

integral->arm();
integ1 = integral->integrate(0,z0);
integ2 = integral->integrate(0,z1);
delete integral;
return (integ2-integ1)*angle;
}

/*!
  As the power spectrum is encoded in a SplineWeb, one can at any
  time draw a section through the web. For convenience, this function
  will return a section along the tau-direction, i.e. the power spectrum
  at conformal time tau.
  If called with only two arguments
  \code
  createPower(0,"mypower");
  \endcode
  it returns a spline holding the power spectrum with first spectral index, i.e. spectral
  index number 0 and value InitialPower[0]

  Note: since v4.0, we store the logarithm in the powerwebs. Hence, to get a plot version,
  this routine exp()'s the argument. The reason: we first interpolate then exp(), cause Spline
  interpolation of the cdm-powerspectrum is numerically difficult, if sampling not dense 
  enough. However, interpolating the logarithmic powerspectrum is simple. 

  We therefore sample more points than previously. The underlying numerical integration
  is of course still the same at the same k values.

  The PowerNormalization is set in AnalyzeThis such that what you get from createPower()
  and the routines like dumpPower() returns a power spectrum spline with units
  P(k) Mpc^{-3} h^3 as a function of k/h in Mpc^{-1}
*/

Spline* Cosmos::createPower(int n, string name, SplineWeb* power, Anchor *a, double tau) {
  if (PowerNormalization[n] == 0.) {
    std::string errorString = "Cosmos::createPower() - PowerNormalization is zero";
    errorString += " AnalyzeThis::rescaleSpectra() should be called before calling createPower.";
    throw Bad_Error(errorString);
  }
  if (tau == 0) tau=tau_0();
  if (power==0) power=power_cdm();
  //  Spline *s = power->createAlongX(tau, name,a);
  Spline *s = power->createAlongX(tau, "createpower");
  s->arm();

//	cout << " start: " << s->start() << " stop: " << s->stop() << endl;

  int steps=400;
  Spline *output = new Spline(steps+1,name,a);
  double lnkh=log(s->start());
  do {
    if (lnkh > s->stop()) lnkh = s->stop();
    double kh = exp(lnkh);   // k/h
    double k = kh *h(); 
    double power = exp(s->fastY(kh)) * scalarSpectrum(k,n) * PowerNormalization[n];
    output->set(kh,power);
    lnkh += (log(s->stop())-log(s->start()) )/steps;
  } while (lnkh < log(s->stop()));
  delete(s);
  return output;
  /*
  for (int i = 0; i < s->size(); i++) {
    double k = s->x(i) *h();  // cause power stored  k/h, only needed for the multiplication
    // multiply with scalar spectrum and normalization. Result is in  (Mpc h^-1)^3
    // and as a function of (k/h)
    s->mul(i, scalarSpectrum(k,n) * PowerNormalization[n]);
  }
  return s;
  */
}

/*! 
  Given a power spectrum spline created from createPower(), return
  capital \Delta, i.e. the density contrast at a given scale k, where k
  is in Mpc^{-1} (not k/h!) and the power spectrum is a function of
  k/h as usual in this code. 
*/
double Cosmos::Pk2Delta2(Spline *P, double k) {
  return k*k*k*P->fastY( k/h() ) / (2*M_PI*M_PI * h() *h() *h() );
}

/*
Spline* Cosmos::createDelta2(string name, SplineWeb* power, double tau, Anchor* a) {
  if (tau==0) tau = tau_0();
  Spline *s = power->createAlongX(tau,name, a);
  double k,k2;
  for  (int i = 0; i < s->size(); i++) {
    k = s->x(i)*h(); // because power is stored for k/h
    k2 = k;
    s->mul(i, k2*k2);
  }
  return s;
}
*/

/*!
  Limber approximation for Lensing. Implementation by Georg Robbers along
  the lines of astro-ph/0502425 
*/
Spline * Cosmos::createLensingPowerLimber( int n, string name, int lstart, int lend, Anchor* a)
{
  //cout << "Using Limber approximation" << endl;
  const double chiRec = tau_0() - tau_ls();
  double tau0 = tau_0();
  Spline *Cpsi_l = new Spline( lend - lstart + 1, name, a );

  for (  int l = lstart; l <= lend; l += ( l<200 )?( l<60 )?( l<15 )?1:1:50:50 )
  {
    double chiIntegral = 0;
    //    double lastchi;
    double chiStepSize = (l<1000)?(l<200)?200:600:600;
    for ( double chi = 1, lastchi = 1; chi <= chiRec; chi += chiStepSize )
    {
      double k = l/chi;
      if ( k < power_psi()->x( 0 ) || k > power_psi()->x( power_psi()->igrid ) ) continue;
      Spline *tmpSpline = power_psi()->createAlongY( k, "powerSpline" );
      tmpSpline->arm();
      double psi2 = pow( tmpSpline->fastY( tau0-chi ),2 );
      double psiPower = PowerNormalization[ n ]*scalarSpectrum( k, n )*psi2/( 2 * M_PI * M_PI*h()*h()*h() );
      double value = 8*M_PI*M_PI*l*psiPower*(chiRec-chi)*(chiRec-chi)/chiRec/chiRec/chi;
      chiIntegral += (chi-lastchi)*value;
      lastchi = chi;
      delete tmpSpline;
    }
    Cpsi_l->set( l, chiIntegral );
    chiIntegral = 0;
  }
  return Cpsi_l;
}
