#include "coupledleaping.h"
#include "exponentialcoupling.h"
#include "spline.h"
#include "quintcosmos.h"
#include "coupledquintcosmos.h"
#include <list>

static bool verbose = false;

CoupledLeaping::CoupledLeaping(CoupledQuintCosmos &c) : Arthur(c)
{
  mQinitial = -1;
  mHasQInitial = false;
  mQinitialFactor = 1.;
  mHasQInitialFactor = false;

  reset();  // allocate splines 
  param.resize(3);
  param[0] = 0.1;
  param[1] = 278;
  param[2] = 1;

  start_invert = 0;
  stop_invert = 700;
  chiEnd = chiCross = phiCross = -100;

  AlphaOfPhi = 0;
  AlphaCanonical = 0;

  kinetic=0;
  CoupledCosmos=&c;
  tuned=false;
  early=false;

  coupling = 0;

  capitalK = new Spline(1010,"LeapCoupcapitalK",&anchor);
  capitalKInverse = new Spline(1010,"LeapCoupcaptialKInverse",&anchor);
}

void CoupledLeaping::reset()
{
  if (verbose)
    cout << "Coupled Leaping::reset()" << endl;
  anchor.kill();
  v = new Spline(1010,"LeapCoupV",&anchor);
  vp = new Spline(1010,"LeapCoupVprime",&anchor);
  vp2 = new Spline(v,"LeapCoupVprime2",&anchor);

  altV = new Spline(1010,"LeapCoupV",&anchor);
  altVp = new Spline(1010,"LeapCoupVprime",&anchor);
  capitalK = new Spline(1010,"LeapCoupcapitalK",&anchor);
  capitalKInverse = new Spline(1010,"LeapCoupcaptialKInverse",&anchor);
}

double CoupledLeaping::Vprime(const double q, const double a, const double) const
{
//X   const double exponential = (*vp)(q/M_p);
//X #warning should this really be handled this way?
//X   int leftIndex = vp->indexToLeft(q/M_p);
//X   int rightIndex = vp->indexToRight(q/M_p);
//X   static const double negInf = -1 * std::numeric_limits<double>::infinity();
//X   if ( isnan(exponential)
//X        && vp->y(leftIndex) < -140 && vp->y(rightIndex) < -140) {
//X     cout << "returning 0 as Vprime at q=" << (q/M_p) << endl;
//X     cout << "w/ left: " << leftIndex << " / "  << vp->y(leftIndex)
//X          << " and right: " << rightIndex << " / " << vp->y(rightIndex) << endl;
//X     exit(0);
//X     return 0.;
//X   }
//X   cout << "For " << (q/M_p) << "    I have " << (*v)(q/M_p) << " and at back: " <<  (*v).back() << endl;
//X   cout << "   and in vp " << (q/M_p) << "    I have " << (*vp)(q/M_p) << " and at back: " <<  (*vp).back() << endl;
//X   cout << "   and log: " << (log((*v)(q/M_p))) << endl;
//X   cout << "   returning: " << -Mp3*exp( (*vp)(q/M_p)) << endl;
//X   cout << "checks were: " << vp->y(leftIndex) << "  -  " <<  vp->y(rightIndex) << endl;
  const double Vp = (*vp)(q/M_p);
  return -Mp3*Vp;
}

void CoupledLeaping::prepare()
{
  if (verbose)
    cout << "CoupledLeaping::prepare()" << endl;
  //  double start = start_invert;
  double stop = stop_invert;
  double canonical_stop = 400;
  double step = 0;
  double y[2];
  double dy[2];
  double& phi = y[1];

  double chi = start_invert;
  phi = 0;
  //   cout << " chi -q_1() : " << chi-q_1() << endl;
  kinetic = new Spline(1000,"fieldspline",&anchor);


  do {
    if (chi + step > stop && phi > canonical_stop) break; // step = stop - chi;
    // in the next six lines or so, we find out how large the step_size
    // should be in order not to change the curve to much
    // this is important, cause we've got a very "jumpy" behaviour
    // in the kinetic term

    double step_change = 0.01;  // relative change of y (this is as we flip truely x)
    double PrimeOne;
    PrimeOne = kineticTerm(chi);

    // cout << "phi: " << phi << "  PrimeOne: " << PrimeOne << "  step: " << step << endl;

    step =  phi / PrimeOne * step_change; 
    if (step < 1e-5) step = 1e-5; 

    double MiniStep = step*0.5;
    int count = 0;
    while (kineticTerm(chi + step) / PrimeOne > 2.0 && count++ < 10) {
      step -= MiniStep;
      MiniStep *= 0.5;
      if (step < 1e-5) { step = 1e-5; break;} 
    }

    if (step > 0.1) step = 0.1;

    //    cout << "phi: " << phi << "   PrimeOne: " << PrimeOne <<   "   chi: " << chi << "  step: " << step << endl;

    Miscmath::rungeInt(y, 1, chi, chi+step, 1e-10, step*1e-2, 0, (moDerivs)&CoupledLeaping::leap, *this); 
    chi += step;

    leap(chi,y,dy);  // at start

    capitalK->set(chi, phi);
    capitalKInverse->set(phi, chi);

    v->set(phi,exp(-chi));
    double a_eq = cosmos.omega_relativistic()/cosmos.omega_m();
    double rho_m_eq = cosmos.initialRho_cdm(a_eq);
    if (phiCross < 0 && Mp4*exp(-chi) < rho_m_eq) {
      if (verbose)
        cout << "set phiCross to " <<  phi << "  and chiCross= " << chi << endl;
      phiCross = phi;
      chiCross = chi;
    }
//X     cout << phi << " " <<  chi << endl;
//X     cout << phi  << " " <<  exp(-chi) << " (exp(-chi) /dy[1])= " << (exp(-chi) /dy[1])
//X           << " and log " << log(exp(-chi) /dy[1])<< endl;
//X     cout << "set Vp = " << phi << " to " << log(exp(-chi) /dy[1]) << endl;
    vp->set(phi, exp(-chi) /dy[1]);
    altVp->set(phi, -exp(-chi) / dy[1]);
    kinetic->set(phi,kineticTerm(chi));

//X      cout << "chi: " << chi << "  step: " << step << endl;
//X      cout << "phi: " << phi<<  "  k(chi): " << kineticTerm(chi) <<endl;
//X      cout << "CoupledLeaping::prepare  kanonical field phi = " << phi
//X           << "  non-canonical: " << chi << "   ratio: " << phi / chi << "  leap: " << kineticTerm(chi) << endl;
//X      cout << chi << " kinetic: " << kineticTerm(chi) << endl;
  } while (chi < stop || phi < canonical_stop);

  chiEnd = stop = chi;

  capitalK->arm();
  capitalKInverse->arm();

  v->arm();
  vp->arm();

  altVp->arm();
  altVp->derive(*vp2);
  vp2->arm();

  kinetic->arm();



  //cout <<"CoupledLeaping::~prepare()" << endl;
  //printStatus();
}

void CoupledLeaping::dumpSplines(const std::string postfix) const
{
  capitalK->dump("K"+postfix);  // large K of Arthur, see paper "Natural Quintessence" by Hebecker & Wetterich
  kinetic->dump("kinetic"+postfix);
  v->dump("V"+postfix);
  vp->dump("Vp"+postfix,true);
  vp2->dump("Vp2"+postfix,true);
  altVp->dump("altVp"+postfix);
}

void CoupledLeaping::leap(const double x, const double*y, double* yp)
{
  yp[1] = kineticTerm(x);
}

double CoupledLeaping::canonicalToOriginalField(const double phi) const
{
  if (!capitalKInverse || capitalKInverse->size() < 2)
    throw Bad_Error("CoupledLeaping::canonicalToOriginalField() - call prepare() first.");
  return capitalKInverse->fastY(phi);
}

double CoupledLeaping::originalToCanonicalField(const double chi) const
{
  if (!capitalK || capitalK->size() < 2)
    throw Bad_Error("CoupledLeaping::originalToCanonicalField() - call prepare() first.");
  return capitalK->fastY(chi);
}

double CoupledLeaping::k_minRadAttractor(double omegaq_early)
{
  ExponentialCoupling* expCoupling = dynamic_cast<ExponentialCoupling*>(coupling);
  if (!expCoupling)
    throw Bad_Error("CoupledLeaping::k_minRadAttractor() - assumes constant coupling (i.e. ExponentialCoupling at the moment");
  const double beta = expCoupling->b_original(); // the coupling constant

  const double lambdaEarly = beta*(1.-1./(2.*omegaq_early))+sqrt(beta*beta/(4.*omegaq_early*omegaq_early)+4./omegaq_early);
  const double k_min = 1./lambdaEarly;
  if (verbose)
    cout << "k_min_analyt yields: k_Min = " << k_min << endl;
  return k_min;
}

//! Critical density on the radiation attractor for a true exponential potential
static double rhoCritRDE(const double rho_relativistic, const double lambda, const double beta)
{
  return rho_relativistic/(1.-(4.-beta*(lambda-beta))/(lambda-beta)/(lambda-beta));
}

double CoupledLeaping::initialQ(const double a) const {
  if (mHasQInitial)
    return mQinitial;

  if (!v)
    throw Bad_Error("CoupledLeaping::initialQ() - call prepare() first.");

   // like the uncoupled case:
  list <double> *nullen = v->getZeroList(suppression()*cosmos.initialRho_g(a)/(Mp4),1e-5);
  if (nullen->size() != 1)
    throw Bad_Error("CoupledLeaping::initialQ() I expect the potential to be single-valued");
  double r =  M_p*nullen->front();
  delete nullen; 

  const double phiAtSplineStart = v->start(); // first phi in Spline for potential, assumed to be deep in radiation domination
  const double lambdaEarly = -log(v->fastY(phiAtSplineStart))/phiAtSplineStart;
//X   cout << "lambda_rde recovered: " << lambdaEarly << endl;

  if (!coupling)
    throw Bad_Error("CoupledLeaping::initialQ() - call setCoupling() first.");

  ExponentialCoupling* expCoupling = dynamic_cast<ExponentialCoupling*>(coupling);
  if (!expCoupling)
    throw Bad_Error("CoupledLeaping::initialQ() - assumes constant coupling (i.e. ExponentialCoupling at the moment");
  const double beta = expCoupling->b_original(); // the coupling constant

  double phi_initial;

  if (cosmos.rho_nuNR0() != 0)
    throw Bad_Error("CoupledLeaping::initialQ() - take massive nu into account");
  const double rhoRDECrit = rhoCritRDE(cosmos.initialRho_nu(a)+cosmos.initialRho_g(a), lambdaEarly, beta);
  const double Omega_potRDE = (4.-3.*beta*(lambdaEarly-beta))/(3.*pow(lambdaEarly-beta, 2));
  phi_initial = -M_p/lambdaEarly*log(Omega_potRDE*rhoRDECrit/cosmos.M_p(4));

  if (mHasQInitialFactor)
    phi_initial *= mQinitialFactor;

  if (verbose) {
    cout << "LKT would've had initialQ of " << r << "  " << (r/M_p) << " M_p" << endl;
    cout << "CoupledLeaping returning initialQ of    " << phi_initial << "  " << (phi_initial/M_p) << " M_p" << endl;
    cout << "Omega_potRDE = " << Omega_potRDE << endl;
  }
  return phi_initial;
}

double CoupledLeaping::initialQDot(const double a)  const
{
  const double phiAtSplineStart = v->start(); // first phi in Spline for potential, assumed to be deep in radiation domination
  const double lambdaEarly = -log(v->fastY(phiAtSplineStart))/phiAtSplineStart;
  ExponentialCoupling* expCoupling = dynamic_cast<ExponentialCoupling*>(coupling);
  if (!expCoupling)
    throw Bad_Error("CoupledLeaping::initialQDot() - assumes constant coupling (i.e. ExponentialCoupling at the moment");
  const double beta = expCoupling->b_original(); // the coupling constant

  const double rhoRDECrit = rhoCritRDE(cosmos.initialRho_nu(a)+cosmos.initialRho_g(a), lambdaEarly, beta);
  const double Omega_potRDE = (4.-3.*beta*(lambdaEarly-beta))/(3.*pow(lambdaEarly-beta, 2));
  const double VEDE = rhoRDECrit*Omega_potRDE;
  const double qDot = 4.*a*sqrt(VEDE/fabs(4.-3.*beta*(lambdaEarly-beta)));
  if (verbose) {
    cout << (4.-3.*beta*(lambdaEarly-beta)) << " " << beta << " and lambdaEarly " << lambdaEarly << endl;
    cout << "Coupled Leaping is returning initialQDot of: " << qDot  << " , i.e. " << (qDot/M_p)
      << ", uncoupled LKT would return " << Arthur::initialQDot(a) << " i.e. " << (Arthur::initialQDot(a)/M_p) << endl;
    cout << "comparing Vs: " << (rhoRDECrit*Omega_potRDE) << " and V_true: " << V(initialQ(a),a) << endl;
  }
  return qDot;
}

void CoupledLeaping::printStatus() {
  cout << "This is a Coupled Leaping kinetic term quintessence" << endl; 
  if (tuned) { 
    cout << "Tuned to ";
    if(early) cout << "early beta with beta= " << tuneparameter << endl;
    else cout << "late beta with beta = " <<  tuneparameter << endl;
  }
  Quintessence::printStatus();
}


double CoupledLeaping::tuneHelper2(double k) {
  cout << "#"; cout.flush();
  // cout << "CHECKING 2: " << k<< endl;
  param[1] = k;
  cosmos.reset();
  reset();
  //CoupledCosmos->tuneCDM(); // Mona
  // cosmos.quickhistory(); // Mona
  cosmos.history(false); // Mona
  //cout << "phi_0 is: " << k << "  omega_q: " <<  cosmos.omega_q() << " requested: " << cosmos.omega_q(false) << endl;
  //  return cosmos.omega_q() - cosmos.omega_q(false); // Mona: tune rho_q instead
  cout << "rhoq_set: " << cosmos.rho_0()*cosmos.omega_q(false) << " rhoq_obtained: " << cosmos.tau2rho_q(cosmos.tau_0()) << endl;
  return cosmos.tau2rho_q(cosmos.tau_0()) - cosmos.rho_0()*cosmos.omega_q(false);

}


double CoupledLeaping::tuneHelper(double k) {
  cout << "%";
  cout.flush(); // so that % will be given out at once!
  //  cout << "CHECKING: " << k<< endl;
  param[0] = k;
  cosmos.reset();
  reset();
  //  cosmos.quickhistory(); // Mona - otherwise it doesn't use the coupling_correction!!
  cosmos.history(false); // Mona

  cout << "k_min is: " << k << " " << "cosmos.omebar() is: " << cosmos.omebar() << "  TuneHelpOmebar is: " << TuneHelpOmebar << endl;

 // Mona start -- give out some things
    
  cout << "tau_ls is: " << cosmos.tau_ls() << endl;
  /*
  ofstream check_omega;
  check_omega.open("check.dat");
  for(double zeit = 6.e-7; zeit < 4500; zeit += 10){
    check_omega << zeit << "  " << cosmos.tau2rho_q(zeit)/cosmos.tau2rho(zeit) << "  " << cosmos.tau2AvOmega_q(zeit) << endl; // Mona: xdriver: cosmos is an object of the type CoupledQuintCosmos! ??warum kann ich hier cosmos verwenden, muss aber in Zeile 313 CoupledCosmos->getInitialB() verwenden???
  }
  check_omega.close();
  
  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  */
 // Mona end


  cout.flush();
  return cosmos.omebar() - TuneHelpOmebar;
}


// tune to the right value; done iteratively, since with beta the expansion history changes which
// changes the tuned value of the Omegals and therefore k(phi), changing the true coupling. 
void CoupledLeaping::tuneQuintAndCoupling(double Omegabar, double coupling, bool bearly, double genauigkeit){
  double counter = 0; //Mona
  tuneQuintessence(Omegabar); // fills the Bspline, needed in the following
  if (bearly){
    double currentCoupling=CoupledCosmos->bspline->front();
    //  while(abs(currentCoupling-coupling)/abs(coupling)>1e-3){
    while(fabs(currentCoupling-coupling)/fabs(coupling)>genauigkeit){
      counter += 1; // Mona
      vector<double> vec(1);
      vector<double> temp(1);
      // vec[0]=param[0]*coupling;
      temp=CoupledCosmos->couple->getParameters(); // Mona: getParameters() is a fcn from kineticcoupling.cc. it sets vec[0]=_B and returns vec
      // cout << temp[0] << endl;
      double kearly=temp[0]/CoupledCosmos->bspline->front();
      vec[0]=kearly*coupling;
      CoupledCosmos->couple->setCouplingParameters(vec);
      // cout << param[0]*coupling << " " << param[0] << endl;
      tuneQuintessence(Omegabar);
      currentCoupling=CoupledCosmos->bspline->front();
      //      cout << "currentCoupling is: " << currentCoupling << " coupling is: " << coupling << endl;

      if(counter > 20) throw Bad_Error("too many iterations in tuneQuintAndCoupling!");
    }
    early=true;
  }
  else{ 
    double currentCoupling=CoupledCosmos->bspline->back();
    while(fabs(currentCoupling-coupling)/fabs(coupling)>1e-3){
      vector<double> vec(1);
      vector<double> temp(1);
      // vec[0]=param[0]*coupling;
      temp=CoupledCosmos->couple->getParameters(); 
      // cout << temp[0] << endl;
      double klate=temp[0]/CoupledCosmos->bspline->back();
      vec[0]=klate*coupling;
      CoupledCosmos->couple->setCouplingParameters(vec);
      // cout << param[0]*coupling << " " << param[0] << endl;
      tuneQuintessence(Omegabar);
      currentCoupling=CoupledCosmos->bspline->back();
      // cout << "Current" << currentCoupling << " target: " << coupling << endl;
    }
    early=false;
  }
  tuned=true;
  tuneparameter=coupling;
}

void CoupledLeaping::tuneQuintessence(double Omebar) {


 
  //  double k_middle = sqrt(Omebar/3.0);  // Mona: without coupling

  // Mona start
  
 // get right k_middle from omega_q in Christof's attractor, which takes the coupling into account
  //first get the b_frueh from cqc: this is really b, not btilde!!


  double b_frueh = CoupledCosmos->getInitialB(); // Mona
  // cout << "b_frueh: " << b_frueh << endl; // Mona

  // cout << "binitial_set is: " << b_frueh << ", Omebar_set is: " << Omebar << endl;
  double k_middle = (-b_frueh + 2*b_frueh*Omebar - sqrt(b_frueh*b_frueh+12*Omebar))/(2*(-3-b_frueh*b_frueh+b_frueh*b_frueh*Omebar));  // including coupling 
  
  // Mona end

  // cout << "k_middle is: " << k_middle << endl;  // Mona
 TuneHelpOmebar = Omebar;

 double tol = 1e-4;
 // double tol = 1e-6; // Mona: have increased the interval for zbrent => have effectively decreased the precision!!
 double dk = 1.9/6.0 * Omebar / k_middle;

 // double dk = sqrt(1.0/(12*Omebar)) * Omebar*0.1;

 // double a = k_middle - dk; // Mona -- better to fix the lower bound a
 // double b = k_middle + dk; // Mona
 // double a = 1.e-5; // Mona -- this one works!!!!!!!!!!!!!!!!!!!!
 double a = 1.e-3; // Mona -- works quicker!!!!!!!!!!!!!!!!!
 double b = k_middle + dk; // Mona 
 //cout << "dk: " << dk << "  a: " << a << "  b: " << b << endl;


 /*
 ofstream of;  // Mona
 of.open("tunehelper.dat");  // Mona
 for(double testx = a; testx <= b; testx += 1e-3)  // Mona
   of << testx << " " << tuneHelper(testx) << endl;  // Mona
 of.close();  // Mona
 */

 double k = Miscmath::zbrent( (moSingle) &Arthur::tuneHelper, *this, a,b, fabs(a-b)*tol);
 param[0] = k;
 
 //cout << "a: " << a << "b: " << b << endl;
 // cout << "k is: " << k << endl;
 
 
 // +++

 a = 270;
 b = 290;
 
 
 double fi = Miscmath::zbrent( (moSingle) &Arthur::tuneHelper2, *this, a,b, fabs(a-b)*tol);
 cout << endl;
     // cout << "phi_0: " << fi << endl;
 param[1] = fi;
 cosmos.reset();
 reset();
 
 //cosmos.history(true);
 //cout << "DONE HISTORY l2 " << endl;
 return;
} 

// VM start 


void CoupledLeaping::tuneOmegaQ(double accuracy){

  //  cout << "this is tuneOmegaQ()" << endl;

  cosmos.reset();
  reset();
  cosmos.history(false);

  int Qcount = 0;
  double phinot = param[1];
  double phinot_middle = param[1];
  double phinot_later;
  double phinot_earlier;

  double rhoQset = cosmos.rho_0()*cosmos.omega_q(false);
  double rhoQobtained = cosmos.tau2rho_q(cosmos.tau_0());
  double rhoQobtained_middle = cosmos.tau2rho_q(cosmos.tau_0());
  double rhoQobtained_upper;
  double rhoQobtained_lower;

  if(fabs(rhoQset-rhoQobtained)/fabs(rhoQset) < 0.05){
    // phinot_later = phinot + 100; // use for chains
    phinot_later = phinot + 100;
    //    phinot_earlier = phinot - 2;
    phinot_earlier = phinot - 10; // Mona: for constant b = 0.5
  }
    else{
    phinot_later = phinot + 10;
    //    phinot_earlier = phinot - 1;
    phinot_earlier = phinot - 10; // Mona: for constant b = 0.5
  }
   
   
  while(fabs(rhoQobtained-rhoQset)/fabs(rhoQset) > accuracy){
    
    cout.flush();
    cout << "-";
    

    //cout << "phinot_later is: " << phinot_later << endl;
    //cout << "phinot_earlier is: " << phinot_earlier << endl;
    //cout << "rhoq_set is : " << rhoQset << " rhoq_obtained is: " << rhoQobtained_middle << endl;
    // cout << "phi_0 is: " << param[1] << endl; 

    Qcount += 1;

    // check rho_q today for three different value of phi_0:

    param[1] = phinot_middle;
    cosmos.reset();
    reset();
    cosmos.history(false);
    rhoQobtained_middle = cosmos.tau2rho_q(cosmos.tau_0());
    // cout << "rhoQobtained_middle is: " << rhoQobtained_middle << endl;

    param[1] = phinot_earlier;
    cosmos.reset();
    reset();
    cosmos.history(false);
    rhoQobtained_upper = cosmos.tau2rho_q(cosmos.tau_0());
    //cout << "rhoQobtained_upper is: " << rhoQobtained_upper << endl;

    param[1] = phinot_later;
    cosmos.reset();
    reset();
    cosmos.history(false);
    rhoQobtained_lower = cosmos.tau2rho_q(cosmos.tau_0());
    //cout << "rhoQobtained_lower is: " << rhoQobtained_lower << endl;

    // comment for chains!!!:
    /*
    ofstream check_omega;
    check_omega.open("check.dat");
    for(double zeit = 6.e-7; zeit < 11700; zeit += 10){
      check_omega << zeit << "  " << cosmos.tau2rho_q(zeit)/cosmos.tau2rho(zeit) << endl; 
    }
    check_omega.close();
    */
    //
    


    // choose the interval in which rhoQset lies:

    if(rhoQset > rhoQobtained_upper) throw Bad_Error("rhoq_today out of range -- decrease phinot_earlier in tuneOmegaQ()");
    else if(rhoQset < rhoQobtained_lower) throw Bad_Error("rhoq_today out of range -- increase phinot_later in tuneOmegaQ()");
    else if(rhoQobtained_lower <= rhoQset && rhoQset <= rhoQobtained_middle){
      phinot_earlier = phinot_middle;
    }
    else if(rhoQobtained_middle < rhoQset && rhoQset <= rhoQobtained_upper){
      phinot_later = phinot_middle;
    }

    //  cout << "phinot_later is: " << phinot_later <<  " phinot_earlier is: " << phinot_earlier << endl;

    phinot_middle = fabs(phinot_later + phinot_earlier)/2.;
    //  cout << "phinot_middle  is: " << phinot_middle << endl;

    // set phinot and rhoQobtained to the values in the middle of the new interval:
    param[1] = phinot_middle;
    cosmos.reset();
    cosmos.history(false);
    rhoQobtained = cosmos.tau2rho_q(cosmos.tau_0());

    if(Qcount > 20) throw Bad_Error("too many iterations in tuneOmegaQ()");  

  }

  cout << "tuneOmegaQ set phi_initial to: " << param[1];
  cout << endl;
}


void CoupledLeaping::tuneCosmos(){

  //CoupledCosmos->tuneCDM(0.001);
  tuneOmegaQ(0.001);

  //CoupledCosmos->tuneCDM(0.001);
  tuneOmegaQ(0.001);

  //CoupledCosmos->tuneCDM(0.001);
  tuneOmegaQ(0.001);

  //  CoupledCosmos->tuneCDM(0.001);
  //  tuneOmegaQ(0.001);
}

// VM end
