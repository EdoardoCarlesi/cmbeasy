#include "arthur.h"
#include "spline.h"
#include "quintcosmos.h"
#include <list>
Arthur::Arthur(QuintCosmos &c) : Quintessence(c) {
  reset();  // allocate splines 
  param.resize(3);
  param[0] = 0.1;
  param[1] = 278;
  param[2] = 1;

  Mp2 = M_p*M_p;
  Mp3 = M_p*Mp2;
  Mp4 = Mp2*Mp2;

  start_invert = 0;
  stop_invert = 400;

  AlphaOfPhi = 0;
  AlphaCanonical = 0;  
}

void Arthur::reset() {
  anchor.kill();
  v = new Spline(1010,"ArthurV",&anchor);
  vp = new Spline(1010,"ArthurVprime",&anchor);
  vp2 = new Spline(v,"ArthurVprime2",&anchor);

  altV = new Spline(1010,"ArthurV",&anchor);
  altVp = new Spline(1010,"ArthurVprime",&anchor);
}


double Arthur::kineticTerm(const double q) const { 
  return (k_min() + 1.0) + tanh(beta()*(q - q_1()));   
}

void Arthur::prepare() {
  //  double start = start_invert;
  double stop = stop_invert;
  double step = 0;
  double y[2];

  // cout <<"Arthur::prepare()" << endl;


  //Spline *forFun = new Spline(2000,"forFun::arthur");
 
  double x = start_invert;
  y[1] = 0;
  double dy[2]; // just to sniff around 
  // v->set(x,1); PARDON ?
  //   cout << " x -q_1() : " << x-q_1() << endl;
  do {
    if (x + step > stop) break; // step = stop - x;
    // in the next six lines or so, we find out how large the step_size
    // should be in order not to change the curve to much
    // this is important, cause we've got a very "jumpy" behaviour
    // in the kinetic term
  
    double step_change = 0.01;  // relative change of y (this is as we flip truely x)
    double PrimeOne;
    PrimeOne = kineticTerm(x);
    step =  y[1] / PrimeOne * step_change; 
    if (step < 1e-5) step = 1e-5;
    
    double MiniStep = step*0.5;
    int count = 0;
    while (kineticTerm(x + step) / PrimeOne > 2.0 && count++ < 10) {
      step -= MiniStep;
      MiniStep *= 0.5;
      if (step < 1e-5) { step = 1e-5; break;}
    } 

    if (step > 0.1) step = 0.1;

    //    cout << "y[1]: " << y[1] << "   PrimeONe: " << PrimeOne <<   "   x: " << x << "  step: " << step << endl;


    Miscmath::rungeInt(y, 1, x, x+step, 1e-10, step*1e-2, 0, (moDerivs)&Arthur::leap, *this);
    x += step;
    
    leap(x,y,dy);  // at start
 
    //forFun->set(x,y[1]);
    v->set(y[1],exp(-x));
    vp->set(y[1], log(exp(-x) /dy[1]));
    altVp->set(y[1], -exp(-x) / dy[1]);
    

    
    // cout << "Arthur::prepare  kanonical field phi = " << y[1] << "  non-canonical: " << x << "   ratio: " << y[1] / x << "  leap: " << kineticTerm(x) << endl;
    
    
    // cout << x << " kinetic: " << kineticTerm(x) << endl;
  } while (x < stop);
  
 
  //cout << "Letztes Lambda von Arthur: " << -log(v->back())/v->stop() << endl;

  //forFun->arm();
  //forFun->dump("arthurK");  // large K of Arthur, see paper "Natural Quintessence" by Hebecker & Wetterich
  
  v->arm();
  vp->arm();

  //vp->printStatus();
  altVp->arm();
  altVp->derive(*vp2);
  vp2->arm();
  
  /*
  v->dump("arthurV");
  vp->dump("arthurVp",true);
  vp2->dump("arthurVp2",true);
  */
  //altVp->dump("altVp");
  
 
  
 
  //  altVp->dump("altarthur",false);

  //cout <<"Arthur::~prepare()" << endl;
  //    throw Bad_Error("adsD");
  //delete forFun;
  //printStatus();
}

void Arthur::leap(const double x, const double*y, double* yp) { yp[1] = kineticTerm(x); }
 
double Arthur::suppression() const {
  double suppress = 1;
  if (k_min() < 0.1) {
    suppress = 4*k_min()*k_min(); 
  }
  return suppress;
} 

#define DOF 1e3
double Arthur::initialQ(const double a) const {
  list <double> *nullen = v->getZeroList(suppression()*cosmos.initialRho_g(a)/(Mp4*DOF),1e-5);
  if (nullen->size() != 1) throw Bad_Error("Arhur::initialQ() I expect the potential to be single-valued");
  double r =  M_p*nullen->front();
  delete nullen; 
  return r;
}

double Arthur::initialQDot(const double a)  const {
  return sqrt(suppression()*cosmos.initialRho_g(a)/DOF ) *a;
}

void Arthur::printStatus() {
  cout << "Leaping kinetic term" << endl;
  cout <<endl;
  Quintessence::printStatus();
}


double Arthur::tuneHelper2(double k) {
  cout << "#"; cout.flush();
  // cout << "CHECKING 2: " << k<< endl;
  param[1] = k;
  cosmos.reset();
  reset();
  cosmos.history(false);
  return cosmos.omega_q() - cosmos.omega_q(false);
}

double Arthur::tuneHelper(double k) {
  cout << "%"; cout.flush();
  //cout << "CHECKING: " << k<< endl;
  param[0] = k;
  cosmos.reset();
  reset();
  cosmos.history(false);
  return cosmos.omebar() - TuneHelpOmebar;
}

void Arthur::tuneQuintessence(double Omebar) {
 
 double k_middle = sqrt(Omebar/3.0);
 TuneHelpOmebar = Omebar;

 double tol = 1e-4;
 double dk = 0.3/6.0 * Omebar / k_middle; 
 double a = k_middle  -dk;
 double b = k_middle + dk;

 double k = Miscmath::zbrent( (moSingle) &Arthur::tuneHelper, *this, a,b, fabs(a-b)*tol);

 param[0] = k;

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


vector<QPName> Arthur::parameterNames() const {
  vector<QPName> s(3);
  s[0].name = "k_min";
  s[0].tooltip = "Initial value of kinetic term";
  s[0].determined = true;
  s[1].name = "phi_0";
  s[1].tooltip = "Field value at which transition takes place";
  s[1].determined=true;
  s[2].name = "alpha";
  s[2].tooltip = "Controls steepness of transition. Try 0.1 (slow), 10 (steep)";
  s[2].determined=false;
  return s;
}
