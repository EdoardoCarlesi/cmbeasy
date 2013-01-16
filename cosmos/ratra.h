#ifndef RATRA_H
#define RATRA_H

#include "quintessence.h"
#include "mathobject.h"

#include <cmath>

/*!
  Ratra-Peebles inverse Power law potential V = A*phi^alpha. Make sure to call tuneQuintessence() to adjust the constant A such that today the correct amount of quintessence is present.

  \param alpha The inverse power
*/
class Ratra : public Quintessence,Mathobject {
  double A_A;
  double Alpha;
  double Mp4;
  double Mp3;
  double Mp2;

  double w0_wish;

  double wq_rd() const;

 protected:
  double tau_0() const;

 public:
  Ratra(QuintCosmos&);

  virtual double V(const double q,const double a =1) const {
    return A()*pow(q/M_p,-alpha());
  } //!< Potential
  virtual double Vprime(const double q,const double a  =1,const double adotoa=1) const {
    return -alpha()*A()/M_p*pow(q/M_p,-alpha()-1.0);
  } //!< dV/dphi
  virtual double Vprime2(const double q, const double a =1,
                         const double adotoa=1,const double tau=1)  const {
    return alpha()*(alpha()+1)*A() / Mp2* pow(q/M_p,-alpha()-2.0);
  } //!< d^2 V/ dphi^2

  virtual double getKinetic(double q) const { return 1.; }

  double A() const { return A_A; }
  void setA(const double a) { A_A = a; }
  double alpha() const { return Alpha; }
  void setalpha(const double a) { param[0] = Alpha = a; }
  double estimatedRhoQ(double a) const;
  virtual Type type() { return ipl;}
  virtual void setParameters(const vector<double> &);
  double initialQ(const double a) const;
  virtual double initialQDot(const double a)  const;
  //double suppression() const; //!< return suppression wrt. gammas for initial conditions of Q-field
  double tuneHelper(double);
  void tuneAlphaForw0(const double w0);
  double tuneHelperForw0(double);
  virtual void tuneQuintessence(double omebar=0.0);
  virtual void printStatus();
  virtual vector<QPName> parameterNames() const;
};

#endif
