#ifndef vdecosmos_H
#define vdecosmos_H 

#include "global.h"
#include "spline.h"
#include <iostream>
#include "quintessence.h"
#include "anchor.h"
#include "safevector.h"
#include "controlpanel.h"
#include "cl.h"

#include "cosmos.h"
#include "vectorde.h"
class Vde;

class VdeCosmos : public Cosmos {

  double Omega_vde;
 
  // here come some not so important splines
  // just for visualization and some other stuff
  // keep adding splines for the quantities you
  // want to integrate etc. here :-)

  Spline* Tau2Rho_v;
  Spline* Tau2P_v;
  Spline* Tau2Omega_v;
  Spline* Tau2Omega_k;
  Spline* Tau2W_v;
  Spline* Tau2Omega_vw; // <! internal w_v * O_v
  Spline* a2W;  
  Spline* a2A_0; 
  Spline* a2dWdZ; 
  Spline* a2OmegaWZ;
  Spline* Z2Gamma;  // internal for determining Steinhard't slow or fast condition dw/dlnz ^2 bar 
  Spline* A2Omega_vw;

  Spline* a2Omega_v;
  Spline* a2Omega_k;
  Spline* a2Omega_l;
  Spline* a2Omega_c;
  Spline* a2Omega_b;
  Spline* a2Omega_g;
 
  Spline* a2Rho_v;
  Spline* a2Rho_k;
  Spline* a2Rho_l;
  Spline* a2Rho_c;
  Spline* a2Rho_b;
  Spline* a2Rho_g;

  Spline* A2Omega_v;
  Spline* LogA2Omega_v, *LogA2Omega_vw; //!< for yet another weff
  //enum QuintType { none, exponential, arthur, ratra, rohm, arbitrary,celestine, crossover}; //!< vde type
  enum WeffType {tauWeff, logAWeff, todayWeff}; //!< effective equation of state type, i.e. averaging over what ...
 
 
  double cross_W25;
 public: 
 
  VdeCosmos();
  ~VdeCosmos();
  void initSplines(); // re-implemented
  double tau2rho_v(const double);
 
  Vde *vde;   //!< pointer to vector dark energy object
  double tau2AvOmega_v(const double); //!< return tau-averaged omega_v until time tau

  double tau2p_v(const double); //!< pressure of vde field
  double tau2w_v(const double);  //!< pressure/rho of vde
  double tau2omega_v(const double); 
  
  double tau2rho_vdot(const double); //!< derivative of vde field density w.r.t. tau
  double tau2wdot_v(const double); 
  double v2tau(const double);
  double tau2v(const double);  //!< phi, including the factor of M_p, so you have to divide by M_p to
                                 //!< get  phi in units of M_p
  double tau2vdot(const double);  //!< derivative of phi w.r.t. to conformal time tau

  double omega_vde();  //!< omega qu
  double omega_k(){1-omega_0()-omega_vde();}
  void initVde();	

  void dumpSplines();

  void setOmega_vde(double x);
  void setOmegaH2_vde(double x) { setOmega_vde(x/h2());}   //!< Omega * h^2 of vde
  void setOmega_vde_flat();

  double weff(VdeCosmos::WeffType t) { 
    //if (omega_vde == 0.0) return -1.0;
    switch(t) {
    case tauWeff:
      return Tau2Omega_vw->average(tau_0()) / Tau2Omega_v->average(tau_0());
    case logAWeff:
      return 1.0 /   ( LogA2Omega_vw->average(cross_W25, 0) / LogA2Omega_v->average(cross_W25,0));
    case todayWeff:
      return Tau2W_v->fastY(tau_0());
    default:
      return 1.0;
    }
  }


  //double omebar() { return tau2AvOmega_v(tau_ls()); } //!< Averaged Omega_v until last scattering
  virtual double omesf(double a_tr=1./3.); //!< log a averaged omega_dark. This is relevant for structure formation
  double omega_sf(double a_tr=1./3.) { return omesf(a_tr); } //!< same as omesf(a_tr)
  double omega_0() { return Cosmos::omega_0() + omega_v(); }
  //! set omega_cdm such that a flat universe comes out
  virtual void setOmega_cdm_flat() {
    setOmega_cdm(1.-omega_v()-omega_vde()-omega_b()-omega_nuNR()-omega_nu()-omega_g());
  }
  //void reset(); //! reset cosmos & vde

  void setEnableVdePerturbations(bool b); //!< hmmm. Not really used. Obsolete ?
  double initialRhoVde(double);
  double initialRho_k(double){;}
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

#endif

