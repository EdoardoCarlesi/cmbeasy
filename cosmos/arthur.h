#ifndef ARTHUR_H
#define ARTHUR_H

#include "quintessence.h"
#include "mathobject.h"
#include "spline.h"


class Spline;

/*! 
  Quintessence subclass for generic 
  exponential potentials and a none-canonical
  kinetic term. 

  Of course this term can easyly be transformed
  to a exponential potential with different
  slopes etc. 

  The kinetical term in this realization is
  
  (kmin+1) + tanh(phi - phi1) 
  
  In practice, the kinetical term is evaluated
  when prepare() is called and the potential
  for a canonical kinetic term is constructed.

  Parameters are kmin and phi1 with 
  setParameters() setting them in this order.

  Interesting parameter sets for h = 0.65, Omega_0=1

  0.1   277      -> Omega_quintessence0 = 61.001%
  0.1   277.03 -> Omega_quintessence0 = 59.943%
  0.1   277.028 -> Omega_quintessence0 = 60.009%

  leading to Omega_q(equ) approx 3e-2

  0.2  277.1025  -> Omega_quintessence0 = 60.004%

  leading to Omega_q(ls) 0.124, w_0 = -0.77

  0.3  277.246 -> Omega_quintessence0 = 59.9987%

  leading to Omega_q(ls) 0.279, w_0 = -0.71


 

*/
class Arthur : public Quintessence,  public Mathobject {

 protected:
  Spline* v, *vp, *vp2;
  Spline* altV,*altVp;
  Spline *AlphaOfPhi, *AlphaCanonical;
  Anchor anchor;
  double Mp4;
  double Mp3;
  double Mp2;

  double TuneHelpOmebar;

 
 public:  

  double start_invert, stop_invert;
  
  Arthur(QuintCosmos&);
  virtual void reset();
  virtual Type type() { return leaping; }
  Spline* vSpline() const { return v; }
  double V(const double q,const double a =1) const { return Mp4*(*v)(q/M_p); }                 //!< Potential
  double Vprime(const double q,const double a  =1,const double=1) const {return -Mp3*exp( (*vp)(q/M_p)); }         //!< dV/dphi
  double Vprime2(const double q,const double a =1,const double adotoa=1,const double tau=1)  const {return Mp2*(*vp2)(q/M_p); }       //!< d^2 V/ dphi^2
  
  virtual double initialQ(const double) const;       //!< phi(a_initial )
  virtual double initialQDot(const double) const;  //!< d phi/ dtau (a_initial)
  
  double suppression() const;

  /*!
    The two parameters of the model: k_min and q_1
    are set via setLambda(k_min, q_1)
  */

  virtual void prepare();
  void leap(const double, const double *, double*);
  virtual double kineticTerm(const double q) const;
  virtual void printStatus();
  virtual double k_min() const { return param[0];}
  virtual double q_1() const { return param[1];}
  virtual double beta() const { return param[2];} 
  double tuneHelper(double);
  double tuneHelper2(double);
  virtual void tuneQuintessence(double omebar=0.0);
  virtual vector<QPName> parameterNames() const;
  virtual string name() const { return "LKT";}
};


#endif 
