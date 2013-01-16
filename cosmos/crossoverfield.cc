#include "crossoverfield.h"
#include "spline.h"
#include "quintcosmos.h"
#include "minmax.h"
#include <list>

CrossoverField::CrossoverField(QuintCosmos &c) : Arthur(c)  {
  reset();
  param.resize(5);
  param[0] = 10; // E
  param[1] = 20; // J
  param[2] = 0.1; // C
  param[3] = 0.1; // D
  param[4] = 278; // phi_crit
  Tuning = false;
  callFromHistory=false;
}

void CrossoverField::reset() {
  Arthur::reset();
  delta = new Spline(1000,"CrossoverField::delta",&anchor);
}

double CrossoverField::kineticTerm(const double q) const{ 
  //  cout << "kinetic term for q " << q << "  " <<  sqrt(delta->fastY(q))/2 << endl;
  //cout << "the right one" << endl;
  return sqrt(delta->fastY(q))/2;
}

 
void CrossoverField::prepare() {
  // if the flow equation for delta would be  d delta / d ln(chi) = E*delta^2 then 
  // the solution would be delta = 2/(E ( phi_critical - phi)) 
  // where phi = 2 ln (chi/M_p)
  // hence, given phi_critical, we can obtain delta at any other value of phi 
  // analytically (well, for phi << phi_crit , this should be ok).

  cout << "CROSSOVER: prepare()" << endl;

  double E = param[0];
  //double J = param[1];
  // double C = param[2];
  // double D = param[3];
  double phi_critical = param[4];
  double phi = 1;
  double del = 2 / ( E*(phi_critical - phi));


  // now we integrate the flow equation 
  double realstep = 1e-1;
  double y[2] ,dy[2];
  double last=del;
  y[1] = del;
  delta->set(phi, del); 
  double large_delta = 300; // field value at which delta crosses 1 
  do {
    Miscmath::rungeInt(y, 1, phi, phi+realstep, 1e-10, realstep*1e-2, 0, (moDerivs)&CrossoverField::flowWrapper, *this);
    phi += realstep;
    delta->set(phi, y[1]);
    if (y[1] > 1 && large_delta == 300) large_delta = phi;

    double change = fabs((y[1] - last) / y[1]);
    last = y[1];
    if (change > 0.01) realstep *= 0.5; else realstep *= 1.8;
    realstep = min(realstep,1.0);

    //cout << change <<" realstep: " << realstep << " phi: " << phi << "  delta:  " << y[1] << "  k: " << sqrt(y[1])*0.5 << endl;
    // cout << "   2/(E/2*(phi - phi_crit)): " << 1.0 / (0.5*E*(phi_critical - phi) );
    flow(phi, y, dy);
    //cout << "     d delta / dphi = " << dy[1] << " diff in flow: " << dy[1] - 0.5*y[1]*y[1]*E << endl;
    //    cout << "phi: "<< phi << "  delta: "<< y[1] << endl;
  } while (y[1] < 1e10 && phi < 300);
  cout << "done: phi is: " << phi << "  phi_critical: " << phi_critical << " realstep : " << realstep << " exp(-phi) = " << exp(-phi) << "  large_delta: " << large_delta <<  endl;
  delta->arm();
  //delta->dump("crossdelta",false);


  start_invert = 1;
  stop_invert = phi -5*realstep;
  stop_invert = min(300.0, stop_invert);


  double totalRho = 3*Mp2*pow(cosmos.H_0_cpm(),2);
  double expo = cosmos.omega_q()*totalRho / Mp4 *0.1;  // safety
  cout << "TOTALRHO : " << totalRho << "  expo: " << expo << "   field-value: " << -log(expo) << endl;
  if (large_delta < 270){
    if (!callFromHistory) throw PhiCritTooLow();
  }// Bad_Error("Crossoverfield::prepare() phi for this model remains to low to account for the small late time energy density"); 
  cout << "STARTINVERT: " << start_invert << "  STOPINVERT: " << stop_invert << endl;
  
  if (!Tuning) getAlpha();

  cout << "Calling Arthur's prepare" << endl;
  Arthur::prepare();
  
  //throw Bad_Error("stopli");
  
}

/*!
  The flow equation for delta (in terms of phi = 2*ln(chi)) :
  
  d delta / d phi = 1/2 * d delta / d ln(chi) = E/2 * delta^2 / ( 1 + J*delta)
*/
void CrossoverField::flow(const double x, const double* y, double *dy) {
  double E = param[0];
  double J = param[1];
  double C = param[2];
  double D = param[3];
  double del = y[1];
  dy[1] = 0.5*E*del*del*(1.0 + C*del) / (1  + J*del*(1+D*del));
  //  dy[1]= 0.5*E*pow(del,gamma) / (1 + J*del);
  if (dy[1] > 1e6) dy[1] = 1e6; // 0.5*E*pow(1e3,gamma) / (1 + J*del);
}


double CrossoverField::estimatedRhoQ(const double a) const {
 double k = kineticTerm(1); 
 double omega = 3*k*k; // very early field -> almost konstant k
 double totalRho = cosmos.initialRho_g(a); // not all, but still....
 return totalRho*omega;
}



double CrossoverField::initialQ(const double a) const {
  double k = kineticTerm(1); 
  double omega = 3*k*k; // very early field -> almost konstant k
  double frac_potential = 1.0/3.0; // to account for the fact that only 1/3 at that time in potential energy 
  
  //  v->dump("asseenbycross");
  
  list <double> *nullen = v->getZeroList( estimatedRhoQ(a)*frac_potential/Mp4 ,1e-5);
  cout << "Crossing should be: " << estimatedRhoQ(a)*frac_potential/Mp4 << endl;
  cout << "size of nullen;" << nullen -> size() << endl;
  if (nullen->size() != 1) throw Bad_Error("Crossoverfield::initialQ() I expect the potential to be single-valued");
 
  double q = nullen->front();
  delete nullen;


  cout << "888888888888888888888888888888888888888888888888" << endl;
  cout << "CROSSOVER :: initialQ: omega: " << omega << endl;
  cout << "total rho: " << cosmos.initialRho_g(a) << endl;
  cout << "estimatedRho() : " << estimatedRhoQ(a) << endl;
  cout << "kinetic term: " << k << endl;

  cout << "q: " << q << "     v(q) : " << Mp4*v->fastY(q) << endl;
  return M_p*q;
}

double CrossoverField::initialQDot(const double a)  const {
  double frac_kinetic = 2.0/3.0;
  cout << "INITIALQDOT: " <<  a*sqrt(2*frac_kinetic*estimatedRhoQ(a)) << endl;
  return a*sqrt(2*frac_kinetic*estimatedRhoQ(a));
} 



void CrossoverField::tuneQuintessence(double Omebar) {
  Tuning = true;
  if (AlphaOfPhi) { delete AlphaOfPhi; AlphaOfPhi = 0;}
  cout << "CROSSOVER TUNE" << endl;
  double RequiredPrecision = 1e-3;
  cosmos.reset();
  reset();
  //double s[5];
  int last=0,now=0;
  vector <double> l(5); 
  bool virgin;
  
  l = parameters();
  l[4] = 100;
  setParameters(l);
  

  do {
    try {
      cout << "TRYING WITH PHI_CRIT = " << l[4] << endl;
      virgin = false;
      l = parameters();
      prepare();
      reset();
    } catch (PhiCritTooLow) {
      cout <<"TOOOOOOOOO LOWWWWWWWWWWWWWWWW" << endl;
      l[4] += 10;
      setParameters(l);
      virgin = true;
      reset();
    }
  } while (virgin);
  

  
  double step_phicrit = 2;

  //  if (Omebar < 5e-3) Omebar = 5e-3;  // at least a little bit, please
 
  
  while (step_phicrit > 0.001) {
    try {
      cout << endl << endl << "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
      cosmos.history(false);
    
      if (fabs(cosmos.omega_q() - cosmos.omega_q(false))/cosmos.omega_q(false) < RequiredPrecision) break;
      
      l = parameters();
      cout << endl << endl;
      cout << "GEFORDERT: " << cosmos.omega_q(false) << "    wirklich: " << cosmos.omega_q() << "   phicrit: " << l[4] << "  step: "<< step_phicrit << endl;
      if (cosmos.omega_q() < cosmos.omega_q(false) )  { l[4] -= step_phicrit; now = -1;} else { l[4] += step_phicrit; now = 1;}
   
      cosmos.reset();
      reset();
      setParameters(l);
      if (last != now) step_phicrit *= 0.5; else step_phicrit  *= 1.3;
      step_phicrit = min(step_phicrit,5.0);
      last = now;
    } catch (PhiCritTooLow) {
      cout << "CAUGHT PhiCritTOOLow in tune !!!!!!!!!!!!!!!!!!!!!!" << endl;
      l[4] += step_phicrit;  // going back to higher value
      step_phicrit *= 0.5; // much smaller step 
      last = 1; 
      cosmos.reset();
      reset();
      setParameters(l);
    }
  }
  cosmos.reset();
  reset();
  Tuning = false;
  cout << "END OF CROSSFIELD-TUNING  " << endl;
} 


void CrossoverField::getAlpha() {

  if (AlphaOfPhi) { delete AlphaOfPhi; AlphaOfPhi =0;}
  AlphaOfPhi = new Spline(2000,"alphaofphi",&anchor);
  
  double alpha = 100;
  double y[2];
  y[1] = alpha;
  double hnext = delta->x(1) - delta->x(0); 
  for (int i = 0; i < delta->size() -1; i++) {
    double phi_start = delta->x(i);
    double phi_end = delta->x(i+1);
    hnext = Miscmath::rungeInt(y, 1, phi_start, phi_end, 1e-8,hnext, 0, (moDerivs)&CrossoverField::betaFunction, *this);
    alpha = y[1];
    //    cout << "phi: " << phi_end << " alpha: " << alpha << endl;
    //note.set(phi_end,alpha);

    AlphaOfPhi->set(phi_end,alpha);

  }
  AlphaOfPhi->arm();
  //note.arm();
  //note.dump("alpha",false);
  


  
}

void CrossoverField::betaFunction(const double phi, const double *y, double *dy) {
  
  double alpha = y[1];
  double del = delta->fastY(phi);
  double J =  param[1];

  double b2 = 0.2;
  double b4 = 40*b2;
  
  double b6 = 0.15;

  dy[1] =  b2*alpha - b4*alpha*alpha;



  double S = 6.0 - del / (1 + J*del); 
  dy[1] -= b6 * alpha*alpha*alpha*S;
  
  
  //  cout << "phi: " << phi << " alpha: " << alpha << " del: " << del << "  S: " << S << endl;
  

}

