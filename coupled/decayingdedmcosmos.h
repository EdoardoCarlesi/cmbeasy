#ifndef decayingdedmcosmos_H
#define decayingdedmcosmos_H 

#include "global.h"
#include "spline.h"
#include <iostream>
#include "quintessence.h"
#include "anchor.h"
#include "safevector.h"
#include "controlpanel.h"
#include "cl.h"
#include "coupling.h"
#include "cosmos.h"
#include "quintcosmos.h"
#include "quintessence.h"
#include "exponentialcoupling.h"

struct BackgroundDEDM {
  double A, rho_g, rho_nu, rho_b, rho_c;
  double rho_c1, rho_c2, t, rs;
  double rho_nuNR, p_nuNR;
  double quint;
};

/*! Class for quintessence coupled non-minimally to a DM specie resulting from 
 * the decay of a primary DM component
*/
class DecayingDEDMCosmos : public QuintCosmos {


 public:
  DecayingDEDMCosmos(Quintessence::Type = Quintessence::none);
  ~DecayingDEDMCosmos();

  Spline* A2Rho_cdm_1;
  Spline* A2Rho_cdm_2;
  Spline* A2Rho_cdm_tot;
  Spline* A2Rho_b;
  Spline* A2Rho_g;
  Spline* A2Rho_q;

  /*!This is where the background cosmology will be computed */
  virtual void propagateHistoryInTau(const double , const double*,double *);

  /*! This defines the coupling. For a different coupling, just create your own Coupling class.
   *  One example for use here is the ExponentialCoupling class
   */
  virtual void setCoupling(Coupling* c){couple=c;}

  /*! Returns a pointer to the coupling class. */
  Coupling* coupling() const { return couple; }

  int counter;

 BackgroundDEDM fillHistorySplinesDEDM(double, const double*, const double*, bool);

  double omebar() { return tau2AvOmega_q(tau_ls()); } 

  double tau2AvOmega_q(const double tau) {  return Tau2Omega_q->integrate(Tau2W_q->start(),tau)/(tau - Tau2W_q->start()); }

  virtual void history(bool inform=false);

  virtual void findSpecialMoments(); //!<  more robust than the one in cosmos

  virtual void seth(double h);

  virtual void printStatus(ostream& o=cout);

  //double omega_cdm_1(){return b.rho_c1/rho_crit();}
  //double omega_cdm_2(){return b.rho_c2/rho_crit();}
 
  void reset();

  void dumpAll();

  /*! This sets the constant coupling that are needed for some quintessence classes. Be careful
   since this is not necessarily the coupling used in your Coupling class unless you have set it
  to this value (aplogies for any inconvenience)*/
  virtual void setInitialB(double b) {_B=b;}
  
  /* This gives the constant value of the coupling. However, this may not be the same as what is contained
   in your instance of Coupling. */
  virtual double getInitialB() {return _B;}
  
  double _B;
  Spline* bspline;
  Spline* bprimespline;
  Spline* bprimeprimespline;
  Coupling *couple;

  void initSplines();

	Anchor anchor;
	double expDecay(double, double);

	double tau_dm; // Decay constant for the first dm specie into the second
	// Check out the units of measure
	void setTauDM(double tau) {
	tau_dm = tau;
	//cout << " set TauDM : " << tau_dm << endl; 
	}


};

#endif
