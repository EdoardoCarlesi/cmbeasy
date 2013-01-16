#include "ratra.h"
#include "quintcosmos.h"

#include <iomanip>

Ratra::Ratra(QuintCosmos &c) 
  : Quintessence(c) , A_A(1), w0_wish(100)
{
  param.resize(1);  // only two parameters for this
  param[0] = 1;  // power 4...  
  //param[1] = 280.6406251;  // q_1 = 278  sensible values

  Mp2 = M_p*M_p;
  Mp3 = M_p*Mp2;
  Mp4 = Mp2*Mp2;

}

double Ratra::wq_rd() const {
    double beta = 2 / alpha();
    return (0.333 - beta) / (1+beta);
}

/*!
  Returns an estimate of rho_q at scale factor
  a, uses formulae for w during radiation and
  matter domination and scales back today's rho_q
  using this
*/
double Ratra::estimatedRhoQ(const double a) const  {
  double rho_0 = 3*Mp2 *pow(cosmos.H_0_cpm(), 2);  // today
  double beta = 2 / alpha();  // shortcut
  double wq_md = - beta / (1 + beta);  // w_q during matter dom.
  double a_eq = cosmos.omega_relativistic() / cosmos.omega_m(); // equality a 
  double rho_q = rho_0 * pow(a_eq, -3*(1+wq_md));  // scale until equality
  rho_q *= pow(a / a_eq, -3*(1+wq_rd()));  // scale after equality
  return rho_q*0.8;  // return sort of lower value
}

/*!
  initial field value at scale factor a. Will be close to tracker
  value. Actually, we just take q = qdot * tau(a)
*/
double Ratra::initialQ(const double a) const {
  //double rho_q = estimatedRhoQ(a);
  // double V = (1 - wq_rd()) / 2 * rho_q; 
  //double result = M_p*exp(-log(V / A() ) / alpha());
  // cout << "INTIALQ: first estimate: " << result << endl;

  double rho_0 = 3*Mp2 *pow(cosmos.H_0_cpm(), 2);  // today
  double a_eq = cosmos.omega_relativistic() / cosmos.omega_m(); // equality a 
  rho_0 *= pow(a_eq,-3);  // md
  rho_0 *= pow(a/a_eq,-4); // rd
  
  double Hubble = sqrt(rho_0/3.0) / M_p;
  //  cout << "estimated conformal time at a = " << a << "  : " << 1/(a*Hubble) << endl;
  double conformaltime =  1/(a*Hubble);
  double result = conformaltime*initialQDot(a);
  //  cout << "INITIALQ: second estimate " << result << endl;
  return result;
}

double Ratra::initialQDot(const double a)  const {
  double rho_q = estimatedRhoQ(a); 
  double T = (1 + wq_rd())/2 * rho_q; // kinetic energy
  return sqrt(2*T) *a; 
}

void Ratra::setParameters(const vector<double> &l) {
  Quintessence::setParameters(l);
  Alpha = param[0];
}

void Ratra::printStatus() {
  cout << "Ratra Peebles inverse Power Law" << endl;
  cout << "alpha() " << alpha() << endl;
  cout << "A() " << A()<< ", i.e. in M_p^2-units: " << (A()/pow(baseCosmos::M_p(), 2)) << endl;
  cout <<endl;
  Quintessence::printStatus();
}

double Ratra::tuneHelper(double a) {
//X   cout << "CHECKING: " << scientific  << setprecision(14) << a << endl;
  A_A = a;
  cosmos.reset();
  reset();
  cosmos.history(false);
  return cosmos.omega_q() - cosmos.omega_q(false);
}

void Ratra::tuneQuintessence(double Omebar) {
  A_A =  Mp4*pow(alpha(),alpha())*exp(-277.0);
  bool up = false;
  double direction = 1.0;
  cosmos.history(false);
  if (cosmos.omega_q() < cosmos.omega_q(false) ) up =true;
  double a,b;
  double tmpA = A_A;
  int watchdog = 0;
  for (;;) {
    cosmos.reset();
    reset();
    if (up) direction *= 3; else direction *=0.3;
    a = min(tmpA,tmpA*direction);
    b = max(tmpA,tmpA*direction);
//X     cout << "WATCHDOG: " << watchdog << " direc: " << direction << " a: " << a << "  b: " << b << endl;
    if (up) A_A = b; else A_A=a;
    cosmos.history(false);
//X     cout << "               " << cosmos.omega_q() << "  --  " << cosmos.omega_q(false) << endl;
    if (up && cosmos.omega_q() > cosmos.omega_q(false) ) break;
    if (!up && cosmos.omega_q() < cosmos.omega_q(false) ) break;
    watchdog++;
    if (watchdog > 50) throw Bad_Error("Ratra::tuneQuintessence() impossible to bracket");
    cout << '^';
  }
 cosmos.reset();
 reset();

 //  cout << " direc: " << direction << " a: " << a << "  b: " << b << endl;
  double tol=1e-4;
  double zerois =  Miscmath::zbrent( (moSingle) &Ratra::tuneHelper, *this, a,b, fabs(a-b)*tol);
  A_A = zerois;
  //  cout << "zerois: " << zerois << endl;
  cosmos.reset();
  reset();
  cosmos.history(false);
  double difference = cosmos.omega_q() - cosmos.omega_q(false);
  if (fabs(difference)/cosmos.omega_q(false) > 1e-2) throw Bad_Error("Ratra::tuneQuintessence: Desired precision not reached");
  cout << "Relative Tuning  Precision: " << difference/cosmos.omega_q(false) << "  for omega_q: " << cosmos.omega_q(false) << endl;
  cosmos.reset();
  reset();
  return;
}

double Ratra::tuneHelperForw0(double a)
{
  cout << "CHECKING: " << std::scientific  << std::setprecision(14) << a << endl;
  Alpha = a;
  cosmos.reset();
  reset();
  cosmos.history(false);
  cout << "getting: " << cosmos.tau2w_q(tau_0()) << " while wish was: " << w0_wish << endl;
  return cosmos.tau2w_q(tau_0())-w0_wish;
}

double Ratra::tau_0() const
{
  return cosmos.tau_0();
}


void Ratra::tuneAlphaForw0(double w_0)
{
  w0_wish = w_0;
  //A_A = (1.0/2.0)*cosmos.rho_0()*cosmos.omega_q(false)*(1.0 - w0_wish);
  Alpha = 1;
  bool up = false;
  double direction = 1.0;
  cosmos.history(false);
  if (cosmos.tau2w_q(tau_0()) < w0_wish ) up =true;
  double a,b;
  double tmpAlpha = Alpha;
  int watchdog = 0;
  for (;;) {
    cosmos.reset();
    reset();
    if (up) direction *= 1.1; else direction *=0.9;
    a = min(tmpAlpha,tmpAlpha*direction);
    b = max(tmpAlpha,tmpAlpha*direction);
    cout << "test no: " << watchdog << " direc: " << direction << " a: " << a << "  b: " << b << endl;
    if (up) Alpha = b; else Alpha=a;
    cosmos.history(false);
    cout << "               " << cosmos.tau2w_q(tau_0()) << "  --  " << w0_wish << endl;
    if (up && cosmos.tau2w_q(tau_0()) > w0_wish ) break;
    if (!up && cosmos.tau2w_q(tau_0()) < w0_wish) break;
    watchdog++;
    if (watchdog > 50) throw Bad_Error("Ratra::tuneQuintessenceForw0() impossible to bracket");
    cout << '^';
  }
  cosmos.reset();
  reset();

  cout << " direc: " << direction << " a: " << a << "  b: " << b << endl;
  double tol=1e-4;
  double zerois =  Miscmath::zbrent( (moSingle) &Ratra::tuneHelperForw0, *this, a,b, fabs(a-b)*tol);
  Alpha = zerois;
  cout << "Alpha is: " << zerois << endl;
  cosmos.reset();
  reset();
  cosmos.history(false);
  double difference = cosmos.tau2w_q(tau_0()) - w0_wish;
  if (fabs(difference)/w0_wish > 1e-2) throw Bad_Error("Ratra::tuneQuintessence: Desired precision not reached");
  cosmos.reset();
  reset();
  return;
}

vector<QPName> Ratra::parameterNames() const {
  vector<QPName> s(1);
  s[0].name = "alpha";
  s[0].tooltip="The (absolute value) of the power of the inverse power law";
  s[0].determined = false;
  return s;
}
