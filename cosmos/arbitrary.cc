#include "arbitrary.h"
#include <iostream>
#include "quintcosmos.h"
#include "miscmath.h"


Arbitrary::Arbitrary(QuintCosmos &c) : Quintessence(c) , ValidPhi(false) {
  param.resize(1);  // only two parameters for this
  param[0] = 1;  // power 4...  
  //param[1] = 280.6406251;  // q_1 = 278  sensible values

  Mpinv2 = 1.0/(M_p*M_p); 
}

void Arbitrary::prepare() {
  ValidPhi = false;
  anchor.kill();
  Rho = new Spline(1000,"Arbitrary_rho",&anchor);
  PhiDot = new Spline(Rho,"Arbitrary_phidot",&anchor);

  
  //cout << "Arbitrary::prepare" << endl;
  //cout << "RHO_M 0: " << cosmos.rho_b0() + cosmos.rho_cdm0() << endl;
  // double totalRho = cosmos.rho_b0() + cosmos.rho_cdm0() + cosmos.rho_nu0() + cosmos.rho_g0() + cosmos.tau2rho_v();
   //  totalRho *= cosmos.omega_q()/(1 - cosmos.omega_q());
  double totalRho = cosmos.rho_0();
  totalRho *= cosmos.omega_q();

  double energy[2];
  energy[1]=0;  
  double step=0.001;

  //  ofstream monitor("monitorarbitray.dat");

  // first backward from today into the past
  for (double loga = 0.0; loga > -16; loga -= step) {
    double a = pow(10,loga);
    double b = pow(10,loga - step);
    double en = totalRho*exp(energy[1]);  
    double  pd = a*sqrt( fabs( (1+w(a))*en ) ); // absolute value ! For w < -1 phantom menace 
 
    Rho->set(a,en);
    PhiDot->set(pd);
    Miscmath::rungeInt(energy,1,b, a,1e-8,(a-b)*0.01,0,   (moDerivs) &Arbitrary::integrateWa, *this);
  }
  Rho->flip(); 
  PhiDot->flip(); 
  
  // now forward from today into the future (just a little bit)
  energy[1]=0;  // reset
  double hnext=-1;
  for (double loga = 0; loga < 0.1; loga += step) {
    double a = pow(10,loga);
    double b = pow(10,loga + step);
    if (hnext == -1) hnext = (b-a)*0.01;
    hnext = Miscmath::rungeInt(energy,1,b, a,1e-8,hnext,0, (moDerivs) &Arbitrary::integrateWa, *this);
    double en = totalRho*exp(energy[1]);  
    double  pd = a*sqrt( fabs( (1+w(b))*en ) ); // absolute value ! For w < -1 phantom menace 
    Rho->set(b,en);
    PhiDot->set(pd);
  }
  Rho->arm(Spline::thoseReady);

//X   Rho->dump("arbitrary_rho");
//X   PhiDot->dump("arbitrary_phidot");

}
 
double Arbitrary::q(double a) const {
  if (ValidPhi) return Phi->fastY(a);
  return 0; // It is impossible to know this at the time that quintcosmos asks for it...
}

void Arbitrary::reconstructPhi() {
//X   cout << "reconstruct phi to a = " << cosmos.a_max() << endl;
  ValidPhi = true;
  Phi =  new Spline(1000,"Arbitrary_phi",&anchor);
  Spline PhiDotTau(1000,"phidottau"); 
  for (int k=0; k < PhiDot->size()-1; k++) {
    double a = PhiDot->x(k);
    if (a >= cosmos.a_min() && a < cosmos.a_max() ) {
      PhiDotTau.set(cosmos.a2tau(a), PhiDot->y(k)); // translate to tau
    }
  }
  PhiDotTau.arm();
//X   PhiDotTau.dump("phidottau");
  double phi=0,x1=0,x2=0;
  const int kmax = PhiDotTau.size()-1;
  for (int k=0; k < kmax ; k++) {
    x1=PhiDotTau.x(k);
    x2=PhiDotTau.x(k+1);
    Phi->set(cosmos.tau2a(x1),phi);
    phi += PhiDotTau.integrate(x1,x2);
  }
  Phi->set(cosmos.tau2a(x2),phi);
  Phi->arm();
  double Phi0 = Phi->fastY(1.0);
//X   cout << "PHI0: " << Phi0 << endl;

  Phi->disarm();
  *Phi += -Phi0;
  Phi->arm();

  /*
  ofstream demo("demo.dat");
  int c = 0;
  for (int i =0; i < Phi->size(); i++) {
    double g1 =   -10*    Phi->y(i)/cosmos.M_p()  +240.616;
    double g2 = -11.54* Phi->y(i)/cosmos.M_p() +    236.169;
    double a = Phi->x(i);
    double V = 0.5*(1 - w(a))*rho(a);
    c++;
    if (c == 20) {
      demo << 1/a - 1 << " " << cosmos.a2tau(a) << " " << Phi->y(i)/cosmos.M_p() << " " << V << " " << log(V) << " "<< g1 << " " << " " << g2 <<  " " << w(a) << endl;
      c = 0;
    }
  }
  cout << "Phi(a_equ): " << Phi->fastY(cosmos.a_equ())/cosmos.M_p() << endl;
  */

//X   Phi->dump("arbitrary_phi");
}



double Arbitrary::Vprime(const double q,const  double a,const double adotoa) const  {
  double W = w(a);
  double W1 = wdot(a,adotoa);

  return -( 1.5*(1-W)*adotoa + 0.5*W1/(1+W) )*PhiDot->fastY(a)/(a*a);
}

double Arbitrary::Vprime2(const double q, const double a, const double adotoa, const double tau) const {
  // this is d^2 a / dtau^2 over a from Einstein's equation 
  // double d2a2a = -0.5*( cosmos.tau2rho(tau)/3.0 + cosmos.tau2p(tau) ) *Mpinv2*a*a + adotoa*adotoa; 
  double second = Mpinv2*a*a*( cosmos.tau2rho(tau) - 3* cosmos.tau2p(tau) )/6.0;
  double W = w(a);
  double W1 = wdot(a,adotoa);
  double W2 = wdot2(a,adotoa,second);

  double v = -1.5*(1-W)*(second-adotoa*adotoa*(3.5+1.5*W) );
  v += ( W1*W1/(4*(1+W)) - 0.5*W2 + W1*adotoa*(3*W +2) )  / (1+W);  
 
  return v/(a*a);
}

