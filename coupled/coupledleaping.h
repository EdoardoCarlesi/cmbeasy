#ifndef COUPLEDLEAPING_H
#define COUPLEDLEAPING_H

// #include "quintessence.h"
#include "arthur.h"
#include "mathobject.h"
#include "spline.h"
#include "coupledquintcosmos.h"

class Spline;


class CoupledLeaping : public Arthur {

  double TuneHelpOmebar;
  double mQinitial, mQinitialFactor;
  bool mHasQInitial; //!< true if QInitial has been explicitely set to some value
  bool mHasQInitialFactor; //!< true if QInitialFactor has been explicitely set to some value

  Spline *capitalK, *capitalKInverse; //!< K(phi) and K^-1 from hep-ph/0008205

 public:

  Coupling* coupling;
  void setCoupling(Coupling* c) {coupling = c;}
  double start_invert, stop_invert;
  double phiCross, chiCross, chiEnd;

  double canonicalToOriginalField(const double phi) const;
  double originalToCanonicalField(const double chi) const;

  CoupledLeaping(CoupledQuintCosmos&);

  virtual void reset();

  virtual Type type() { return coupledleaping; }

  double V(const double q,const double a =1) const { return Mp4*(*v)(q/M_p); }                 //!< Potential
  double Vprime(const double q,const double a  =1,const double=1) const; //!< dV/dphi
  double Vprime2(const double q,const double a =1,const double adotoa=1,const double tau=1)  const {return Mp2*(*vp2)(q/M_p); }       //!< d^2 V/ dphi^2

  virtual void   setInitialQ(const double qi) { mHasQInitial = true;  mQinitial = qi; }
  virtual void   setInitialQFactor(const double f) { mQinitialFactor = f; mHasQInitialFactor = (f!=1.); }
  virtual double initialQ(const double) const;       //!< phi(a_initial )
  virtual double initialQDot(const double) const;  //!< d phi/ dtau (a_initial)
  double getKinetic(const double q) const {
    if (!kinetic)
      throw Bad_Error("CoupledLeaping::getKinetic : need to call prepare() first");
    return kinetic->fastY(q/M_p);
  }

  double phiStart() const { return v->start(); }
  double phiStop() const { return v->stop(); }

  virtual void prepare();
  void leap(const double, const double *, double*);

  virtual void printStatus();
  void dumpSplines(const std::string postfix = "-coupledlkdt") const;
  virtual double k_min() const { return param[0]; }
  virtual void setk_min(const double kmin) { param[0] = kmin; }

  virtual double k_minRadAttractor(double omegaq_early);

  virtual double q_1() const { return param[1];}
  virtual double beta() const { return param[2];} 
  double tuneHelper(double);
  double tuneHelper2(double);

  /*!This function can be used to overcome a problem with a non-canonical kinetic term+coupling. It
   is mainly a problem in interpretation, as the coupling changes over time and the physical interpretation
of the raw parameter beta is difficult since the real coupling is determined to be beta/k(phi) with k(phi)
evaluated at a certain epoch. In order to overcome this problem, it is possible to tune the coupling such that the 
coupling is fixed at early or late times, e.g. when the kinetic term does not change. 
 \param Omegabar the amount of early quintessence 
 \param coupling the desired 'true' coupling
 \param bearly if true, tune coupling deep in radiation domination, otherwise tune to today */
  void tuneQuintAndCoupling(double Omegabar, double coupling, bool bearly, double genauigkeit);

  virtual void tuneQuintessence(double omebar=0.0);
  virtual string name() const { return "CoupledLKT";}
  Spline *kinetic;
  bool tuned;
  bool early;
  CoupledQuintCosmos *CoupledCosmos;
  double tuneparameter;


  void tuneOmegaQ(double accuracy); // VM
  void tuneCosmos(); // VM 

};


#endif
