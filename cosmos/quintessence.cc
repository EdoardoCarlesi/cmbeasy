#include "quintessence.h"
#include "global.h"
#include "quintcosmos.h"
#include <iostream>

Quintessence::Quintessence(QuintCosmos &c) :  Q(0), Qdot(0), SpeedOfSound2(1) , cosmos(c),param(vector<double>(1))
{
  Gpi8 = cosmos.Gpi8();
  M_p = cosmos.M_p();
  param[0] = 0.0;
}



double Quintessence::rho(const double a) const {  
	//cout << "Quint-rho: " << rho_kin(a) << endl;
	//cout << "Quint_V: " << V(q(a),a)<< endl;
	return rho_kin(a) + V(q(a),a); }


/*!
  As phiDot is derivative wrt. conformal time, we need scale factor a 
  to get the kinetic energy right. Please do provide this
*/
double Quintessence::rho_kin(const double a) const {return 0.5* qDot(a) * qDot(a) / (a*a);} 

double Quintessence::p(const double a) const {  return rho_kin(a) - V(q(a),a);}

double Quintessence::w(const double a) const { return p(a)/rho(a); }

void Quintessence::setInitial(const double a) {
  setQ(initialQ(a));
  setQDot(initialQDot(a));
}

void Quintessence::touch(const double tau) {
  if (cosmos.validHistory()) {
    setQ(cosmos.tau2phi(tau));
    setQDot(cosmos.tau2phidot(tau));
  } else throw Bad_Error("Quintessence::touch() cosmos does not have valid history. Call history() first");
}


void Quintessence::setParameters(const double p1, const double p2,const double p3,const double p4,const double p5) {
  vector<double> v(5);
  for(;;) { 
    v[0] = p1;
    if (p2==Q_MV) {v.resize(1); break;}
    v[1] = p2;
    if (p3==Q_MV) {v.resize(2); break;}
    v[2] = p3;
    if (p4==Q_MV) {v.resize(3); break;}
    v[3] = p4;
    if (p5==Q_MV) {v.resize(4); break;}
    v[4] = p5;
    break;
  }
  setParameters(v);
}

void Quintessence::setParameters(const vector<double> & l) {
  if (l.size() < param.size()) { 
    parameterResize(); // virtually return here, if it is ok for Q, default will throw an Bad_Error
  }

  if (l.size() > param.size()) {
    printStatus();
    throw Bad_Error("Quintessence::setParameters() too many parameters specified");
  }
  for (unsigned int i =0; i < l.size();i++) param[i] = l[i];
}

void Quintessence::parameterResize() {
  throw Bad_Error("Quintessence::parameterResize()\nLess Parameters specified than needed");
}

void Quintessence::printStatus() {
  cout << "===== Quintessence status =====\n";
  for (unsigned int i = 0; i < param.size(); i++)
    cout << "param[" << i << "]: " << param[i] << endl; 
  //cout << "q:    " << q() << "  in M_p: " << q()/M_p << endl;
  // cout << "qdot: " << qDot() << endl;

  //cout << "V() : " << V() << endl;
  //cout << "Vprime: " << Vprime() << endl;
  //cout << "Vprime2: " << Vprime2() << endl;
  cout << "\n=========================\n";
}


