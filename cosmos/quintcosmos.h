#ifndef quintcosmos_H
#define quintcosmos_H 

#include "global.h"
#include "spline.h"
#include <iostream>
#include "quintessence.h"
#include "anchor.h"
#include "safevector.h"
#include "controlpanel.h"
#include "cl.h"

#include "cosmos.h"
#include "quintessence.h"
//class Quintessence;




/*!
  QuintCosmos is a special kind of Cosmos. 
  In order to accomodate universes including 
  scalar quintessence, the background evolution
  history() has been re-implemented

  In addition, convenient functions to 
  access the background behaviour of 
  quintessence has been added. Also, cluster
  abundance constraints can be calculated. 

  This class is an example of implementing
  a different cosmology derived from the
  standard cosmos class. Most of 
  the history() function is 
  cut-copy-paste from the cosmos version.
  
  In order to support all kinds of different
  Quintessence models, a quintessence
  base class as well as some specific realizations
  are provided.


  It is quite easy to use a cosmos (remember,
  most things in CMBEASY don't really care
  what kind of cosmos, so e.g. handing over
  a Cosmos or QuintCosmos to CmbCalc 
  is the same).
 
  \code
  QuintCosmos cosmos;

  cosmos.setOmega_b(0.05);
  cosmos.setOmega_quintessence(0.5);  // 50% quintessence
  cosmos.setOmega_cdm(0.45); 

  // now we specify the type:
  // we take an Inverse Power Law
  cosmos.setQuintessence(QuintCosmos::Ratra);
  cosmos.setQParameters(2); // inverse power = 2
  

  // now we ask the quintessece to set
  // its remaining parameters properly
  cosmos.tuneQuintessence(); 
  
  // lets do the background
  cosmos.history();
\endcode
If above, a Leaping Kinetic Term would have
been chosen
\code
cosmos.setQuintessence(QuintCosmos::Arthur);
\endcode
we could have used the tuneQuintessence() call
to fine-tune the amount of quintessence at
last scattering (actually the tau-average)
\code
cosmos.tuneQuintessence(0.1); // 10% at last scattering
\endcode
*/
class QuintCosmos : public Cosmos {

  double Omega_quintessence;
#define CACHESPLINE(a) Spline* a; double a##X, a##Y      /* a define for the spline variables for convenience: one refernce and two caching double variables */
 public:  
  CACHESPLINE(Tau2Rho_qpot);
  CACHESPLINE(Tau2Rho_qkin);
  CACHESPLINE(Tau2Rho_q);
  
  CACHESPLINE(Tau2P_q);
  CACHESPLINE(Tau2Vprime);
  CACHESPLINE(Tau2DotVprime);
  
  CACHESPLINE(Tau2VprimeOverRho);
  CACHESPLINE(Tau2DotVprimeOverRho);
  
  CACHESPLINE(Phi2Tau);
  CACHESPLINE(Tau2Phi);
  CACHESPLINE(Tau2PhiDot);

  CACHESPLINE(Tau2Rho_qDot);
  CACHESPLINE(Tau2P_qDot);

  CACHESPLINE(Tau2W_q);
  CACHESPLINE(Tau2WDot_q);
  
#undef CACHESPLINE
  
  // here come some not so important splines
  // just for visualization and some other stuff
  // keep adding splines for the quantities you
  // want to integrate etc. here :-)

  Spline* Tau2Omega_qw; // <! internal w_q * O_q
  Spline* Z2W;  // internal for determining Steinhard't slow or fast condition dw/dlnz ^2 bar 
  Spline* Z2Omega_qkin;  
  Spline* Z2dWdZ; 
  Spline* Z2OmegaWZ;
  Spline* A2Omega_qw;
  Spline* A2Hubble;
  Spline* Z2Hubble_factor;
  Spline* Z2Hubble_ratio;
  Spline* A2Mass_factor;
  Spline* Tau2Mass_factor;
  Spline* A2Hubble_factor;

  Spline* Tau2Omega_q; //<! this may however be upgraded to a cachespline easily :-)
  Spline* A2Omega_q;
  Spline* A2Omega_qkin;  
  Spline* LogA2Omega_q, *LogA2Omega_qw; //!< for yet another weff

  //enum QuintType { none, exponential, arthur, ratra, rohm, arbitrary,celestine, crossover}; //!< quintessence type
  enum WeffType { tauWeff, logAWeff, todayWeff}; //!< effective equation of state type, i.e. averaging over what ...
 
  Quintessence *quint;   //!< pointer to quintessence object
  Quintessence::Type userRequest;  //!< setQuintessence() stores what your wish has been. If omega_q == 0, this wish is not fullfilled :-) see also setOmega_q() which will create your wish quintessence if needed
 
  double cross_W25;
 public: 
 
 
  QuintCosmos(Quintessence::Type = Quintessence::none);
  ~QuintCosmos();
  void initSplines(); // re-implemented
	void dumpExtraSplines();
  double tau2rho_qpot(const double); //!< quintessence potential energy
  double tau2rho_qkin(const double); //!< quintessence kinetic energy 
  double tau2rho_q(const double);

  double tau2vprime(const double); //!< derivative of potential w.r.t field phi
  double tau2dotVprime(const double); //!< d^2/(dtau dphi) V
 
  double tau2vprimeOverRho(const double); 
  double tau2dotVprimeOverRho(const double);
 
  double tau2AvOmega_q(const double); //!< return tau-averaged omega_q until time tau

  double tau2p_q(const double); //!< pressure of quintessence field
  double tau2w_q(const double);  //!< pressure/rho of quintessence
  double tau2p_qdot(const double); //!< derivative of quintessence field pressure w.r.t. tau
  double tau2rho_qdot(const double); //!< derivative of quintessence field density w.r.t. tau
  double tau2wdot_q(const double); 
  double phi2tau(const double);
  double tau2phi(const double);  //!< phi, including the factor of M_p, so you have to divide by M_p to
                                 //!< get  phi in units of M_p
  double tau2phidot(const double);  //!< derivative of phi w.r.t. to conformal time tau

  double omega_q(bool reality=true);  //!< omega quintessence

// Omega CDM today (i.e. not the one set in the initial conditions, which can be larger due to the coupling or varying mass)
  double omega_cdm_final; 

  void setOmega_quintessence(double x);
  void setOmegaH2_quintessence(double x) { setOmega_quintessence(x/h2());}   //!< Omega * h^2 of quintessence
  virtual void setOmega_quintessence_flat() {
     double omegaQuint = 1.-omega_cdm()-omega_b(false)-omega_g();
     if (nuR()>0)
        omegaQuint-=omega_nu();
     if (nuNR()>0)
       omegaQuint-=omega_nuNR(false);
     if (omegaQuint<0.)
       throw Bad_Error("Cosmos::setOmega_quint_flat() - would need negative Omega_cdm");
     setOmega_vacuum(0.);
     setOmega_quintessence(omegaQuint);
  } //!< set omega_quintessence such that a flat universe comes out

  double weff(WeffType t) { 
    if (omega_q(false) == 0.0) return -1.0;
    switch(t) {
    case tauWeff:
      return Tau2Omega_qw->average(tau_0()) / Tau2Omega_q->average(tau_0());
    case logAWeff:
      return 1.0 /   ( LogA2Omega_qw->average(cross_W25, 0) / LogA2Omega_q->average(cross_W25,0));
    case todayWeff:
      return tau2w_q(tau_0());
    default:
      return 1.0;
    }
  }
  double omebar() { return tau2AvOmega_q(tau_ls()); } //!< Averaged Omega_q until last scattering
  virtual double omesf(double a_tr=1./3.); //!< log a averaged omega_dark. This is relevant for structure formation
  double omega_sf(double a_tr=1./3.) { return omesf(a_tr); } //!< same as omesf(a_tr)
  double omega_0() { return Cosmos::omega_0() + omega_q(); }
  //! set omega_cdm such that a flat universe comes out
  virtual void setOmega_cdm_flat() {
    setOmega_cdm(1.-omega_q()-omega_v()-omega_b()-omega_nuNR()-omega_nu()-omega_g());
  }
  void reset(); //! reset cosmos & quintessence

   /*! set the quintessence to q; bypasses the
    *  userRequest mechanism, i.e. this Quintessence will be used whether or not omega_q > 0
    */
  void setQuintessence(Quintessence* q) { delete quint; quint=q; userRequest=quint->type(); }
  void setQuintessence(Quintessence::Type);
  Quintessence::Type quintessenceType() { return userRequest; } //!< return the type of quintessence the user requested  
  Quintessence* quintessence() { return quint; }
  void setQParameters(const vector<double> &l) { if (quint) quint->setParameters(l);}
  void  setQParameters(const double p1, const double p2=Q_MV,const double p3=Q_MV,const double p4=Q_MV,const double p5=Q_MV) { 
    if (quint) quint->setParameters(p1,p2,p3,p4,p5);}
  void setEnableQuintessencePerturbations(bool b); //!< hmmm. Not really used. Obsolete ?

  void tuneQuintessence(double omebar=0, double weff=-0.7); //!< calls quintessence->tuneQuintessence()
  virtual void propagateHistoryInTau(const double , const double*,double *); //!< Return array ziel[] with the tau derivatives of a, rho_r etc for the background equations
  void history(bool inform=false); //!< re-implemented from Cosmos
  virtual void initializeDensities(double* y, double* ydot); //!< zero both integration arrays, and fill y[] with initial values
  virtual Background fillHistorySplines(double tau, const double *y, const double *ydot, bool writeSplines=true); //!< re-implemented from Cosmos
  virtual void finalizeHistorySplines(bool inform); //!< arm the history splines, recalculate phi etc. with knowledge of history if necessary
  double sigma8Omega(WeffType t=tauWeff, int pn=0); //!< return sigma8 * Omega_m ^ gamma from Steinhardt's Cluster Constraints. Takes intialpower and sigma8 from index number n
  double sigma8Omega(double n, double sig8, WeffType t=tauWeff);  //!< same as above, just takes the explizit values of spectral index and sigma8
  double estimateSigma8Q(double sig8, double tau0);
  virtual void printStatus(ostream& o=cout);
};

//#define CACHE(myspline)  if (myspline##X == tau) return myspline##Y; myspline##X = tau; return myspline##Y = (*myspline)(tau)
#define CACHE(myspline)  return (*myspline)(tau)


inline double QuintCosmos::tau2rho_qpot(const double tau)             { CACHE(Tau2Rho_qpot); }
inline double QuintCosmos::tau2rho_qkin(const double tau)             { CACHE(Tau2Rho_qkin); }
inline double QuintCosmos::tau2rho_q(const double tau)             { CACHE(Tau2Rho_q); }
inline double QuintCosmos::tau2p_q(const double tau)  { CACHE(Tau2P_q); }

inline double QuintCosmos::tau2vprime(const double tau)  { CACHE(Tau2Vprime); }
inline double QuintCosmos::tau2dotVprime(const double tau)  { CACHE(Tau2DotVprime); }

inline double QuintCosmos::tau2vprimeOverRho(const double tau)  { CACHE(Tau2VprimeOverRho); }
inline double QuintCosmos::tau2dotVprimeOverRho(const double tau)  { CACHE(Tau2DotVprimeOverRho); }

inline double QuintCosmos::tau2rho_qdot(const double tau)  { CACHE(Tau2Rho_qDot); }
inline double QuintCosmos::tau2p_qdot(const double tau)  { CACHE(Tau2P_qDot); }

inline double QuintCosmos::phi2tau(const double tau)  { CACHE(Phi2Tau); }
inline double QuintCosmos::tau2phi(const double tau)  { CACHE(Tau2Phi); }
inline double QuintCosmos::tau2phidot(const double tau)  { CACHE(Tau2PhiDot); }
inline double QuintCosmos::tau2w_q(const double tau)  { CACHE(Tau2W_q); }
inline double QuintCosmos::tau2wdot_q(const double tau)  { CACHE(Tau2WDot_q); }

#undef CACHE

#endif

