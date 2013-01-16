#ifndef PERTURBATION_H
#define PERTURBATION_H 

class ControlPanel;

#include "cosmos.h"
#include "massiveneutrinos.h"

#define LMAX0 30 // 8  /* Number of Photon Momenta */
#define LMX0 60  //!< maximum multipole moments in photon etc distribution expansion
#define NQMAX0 15 //! From CMBFAST
#define LMAXT0 10  //!< maximum multipole moment of tensor distribution·
#define LMAXNR0 50 //25
#define LMAXNU0 50 //25
#define NQMAX0 15

/*! 
  Perturbation parent class.  
  
  This class provides the frame-work for the perturbation evolution.
  The synchronous gauge and gauge invariant classes are derived from
  it.  
  Maybe most notable is the propagateScalar() function that will call
  virtual functions defined in the sub-classes to actually propagate
  the perturbations through conformal time. 
  PropagateScalar() detects if there is a switch from tightcoupling
  to non-tightcoupling within the interval of propagation. The 
  integration is the split in two parts. In this way, all derivatives 
  are always smooth.

  It also defines convenience functions to access the Cosmos
  object holding the background evolution. This is done to provide
  nicer looking formulas in the perturbation evolution. 

  You can, of course use anything Cosmos provides explicitly. In this
  case, you need to do something like
  \code
  cosmos->Omega_b() 
  \endcode
*/
class Perturbation  : public Mathobject {
  private: 
  double GPI8;
  bool SeedFlag;
 protected:
  double hnext; // the next of for odeint()
  double hnext_tensor; // the same for the tensors
  bool ForceTightCouplingOverrulePhoton;  //!< sometimes, odeint would go a tiny tick further than it should and a switch in the regimes occurs within one integration step. This causes numerical difficulties, because of none-analytic behaviour. This flag forces isTightCoupling() to return MustCoupleTightly 
  bool ForceTightCouplingOverruleBaryon;  //!< sometimes, odeint would go a tiny tick further than it should and a switch in the regimes occurs within one integration step. This causes numerical difficulties, because of none-analytic behaviour. This flag forces isTightCoupling() to return MustCoupleTightly 
  bool MustCoupleTightlyPhoton; //!< see  above ForceTightCoupling.  Don't mess with it: leave it to propagateScalar()
  bool MustCoupleTightlyBaryon; //!< see  above ForceTightCoupling.  Don't mess with it: leave it to propagateScalar()
 protected:
  bool inform; // for testing
  Cosmos* cosmos;
  double *y, *yt;
  double *yprime, *ytprime;
  int nvar,nvart; // number of scalar + tensor equations
  double denl[LMX0+1];
  double dlfdlq[NQMAX0];  //!< dlnf_0 / dlnq array up to NQMAX0-1
  int nqmax;
  double k2;  //!< needs to be set. This is done in initialScalarPerturbations and propagteScalar() / propagateTensor()
  int lmaxg; //!< maximum number of photon momentum l running from 0 ... lmaxg
  int lmaxt; //!< maximum multipole l for tensors
  int lmaxnr; //<! maximum l for relativstic neutrinos
  int lmaxnu; //<! maximum l for relativstic neutrinos


  /*! keep track of whether we already printed out a warning for the case that the source terms
   * apparently have to be calculated at very early times. This bool is set to true after the
   * first such warning, in order to avoid excessive output.
   */
  bool warnOnEarlySource;

 public:
  enum Species { photon, baryon, octopole, multipole, NoSpecies};
  double k;
  double h() const {
    if (cosmos) return cosmos->h(); throw Bad_Error("Perturbation::h() - no cosmos set.");
  }
  double epsw;  // maybe this is just synchronous stuff, maybe should be private
  Perturbation(Cosmos* c);
  virtual ~Perturbation();

  virtual void initialScalarPerturbations(const ControlPanel& control, const double tau)=0;
  virtual void initialTensorPerturbations(); 
  virtual void fderivs(const double tau, const double *y, double *yprime)=0;
  virtual void fderivsTensor(const double tau,  double *y, double *yprime);

  virtual void scalarSources(double tau, double *d, double *dp, double *dk)=0;
  virtual void tensorSources(double tau, double *dt, double *dte, double *dtb);

  void fderivsWrapper(const double tau, const double *f, double *fprime) { fderivs(tau,f,fprime); }
  void fderivsTensorWrapper(const double tau,  double *f, double *fprime) { fderivsTensor(tau,f,fprime);}

  virtual void propagateScalar(double *tau, const double tauend, const double precision);
  void propagateTensor(double *tau, const double tauend, const double precision);

  virtual void setWarnOnEarlySources(bool b) { warnOnEarlySource = b; } //!< if true, warn when cmb sources are calculated earlier than tau = 130
  virtual bool isTightCoupling(const double tau, Species); //!< return true, if at tau, we have tight coupling regime

  pair<double,double> tauTightCoupling(double min, double max, Species); //!< return a tiny intervall in which tightCoupling() switches from true to false

  virtual double delta_c()=0;
  virtual double delta_b()=0;
  virtual double delta_g()=0;
  virtual double delta_n()=0;
  virtual double delta_nr()=0;
  virtual double psi()=0;

  virtual double delta_c_longit() const { throw Bad_Error( "This perturbation subclass does not implement delta_c_longit()"); }
  virtual double delta_b_longit() const { throw Bad_Error( "This perturbation subclass does not implement delta_b_longit()"); }
  virtual double delta_g_longit() const { throw Bad_Error( "This perturbation subclass does not implement delta_g_longit()"); }
  virtual double delta_n_longit() const { throw Bad_Error( "This perturbation subclass does not implement delta_n_longit()"); }
  virtual double delta_nr_longit(double tau) const { throw Bad_Error( "This perturbation subclass does not implement delta_nr_longit()"); }
  virtual double v_c_longit() const { throw Bad_Error( "This perturbation subclass does not implement v_c_longit()"); }
  virtual double v_b_longit() const { throw Bad_Error( "This perturbation subclass does not implement v_b_longit()"); }


  virtual void getReady(const ControlPanel&);

  void setSeedFlag(bool x) { SeedFlag = x; }
  bool seedFlag() { return SeedFlag; } //!< needed by fderivs()
  
  virtual void calcPerturbationNr(const ControlPanel &)=0; //! calculate the number of scalar and tensor equations 
  int scalarPerturbationCount() { return nvar;}
  int tensorPerturbationCount() { return nvart;}
  //int startof_additional() { return scalarPerturbationNr() - additional() +1; }  //!< the first index of y[] that corresponds to an additional perturbation

  
  void initDlnf();  // initializaes dlfdlq, i.e. dlnf_0 / dlnq array up to NQMAX0-1
 
  double tau2z(const double tau) const { return cosmos->tau2z(tau); }
  double tau2a(const double tau) { return cosmos->tau2a(tau); }
 
  double tau2adot(const double tau) { return cosmos->tau2adot(tau); }

  double tau2adotoa(const double tau) { return cosmos->tau2adotoa(tau); }
   
  double soundSpeed(const double tau) { return cosmos->splineCs2->fastY(tau);}

  double tau2p(const double tau) { return cosmos->tau2p(tau); }

  double tau2p_nu(const double tau) { return cosmos->tau2p_nu(tau); }

  double tau2p_g(const double tau) { return cosmos->tau2p_g(tau); }
  
  double tau2p_nuNR(const double tau) { return cosmos->tau2p_nuNR(tau); }
  
  double tau2rho(const double tau) { return cosmos->tau2rho(tau); }

  double tau2rho_cdm(const double tau) { return cosmos->tau2rho_cdm(tau); }
 
  double tau2rho_m(const double tau) { return cosmos->tau2rho_cdm(tau); }

  double tau2rho_b(const double tau) { return cosmos->tau2rho_b(tau); }
  
  double tau2rho_g(const double tau) { return cosmos->tau2rho_g(tau); }
  
  double tau2rho_nu(const double tau) { return cosmos->tau2rho_nu(tau); }
  
  double tau2rho_nuNR(const double tau) { return cosmos->tau2rho_nuNR(tau); }

  double tau2w_nu(const double tau) { return cosmos->tau2w_nu(tau); }
  double tau2w_nuNR(const double tau) { return cosmos->tau2w_nuNR(tau); }

  double tau2grhodot(const double tau) { return cosmos->tau2grhodot(tau); }
  double tau2rhodot_g(const double tau) { return cosmos->tau2rhodot_g(tau); }
  double tau2rhodot_nu(const double tau) { return cosmos->tau2rhodot_nu(tau); }
  
  

  double tau_0() { return cosmos->tau_0();}
  double tau_ls() { return cosmos->tau_0();}
 
  double omega_b() { return cosmos->omega_b();}

  double omega_c() { return cosmos->omega_c();}

  double visibility(const double tau) { return cosmos->visibility(tau); }
  double dvisibility(const double tau) { return cosmos->dvisibility(tau); }
  double ddvisibility(const double tau) { return cosmos->ddvisibility(tau); }
  virtual double opac(const double tau) { return cosmos->opac(tau); }
  double dopac(const double tau) { return cosmos->d2kappa(tau); }

  double expmmu(const double tau) { return cosmos->expmmu(tau); }
  

  double rho_nu0() { return cosmos->rho_nu0(); }
  double rho_0() { return cosmos->rho_0();}
  double nuNR() { return cosmos->nuNR(); }
  double nuR() { return cosmos->nuR(); }
 
  double Gpi8() { return GPI8; }
  static int WatchCount;

  friend class PerturbationTracker;
  friend class TrackerHelper;
};

#endif 
