#ifndef SPEEDYINVARIANT_H
#define SPEEDYINVARIANT_H

#include "perturbation.h"

#include <fstream>
#include <vector>

/*!
  
  Gauge invariant perturbation class implementing the speed ups of
  astro-ph/0503277

  The propagateScalar() routine is a complete rewrite and re-implentation
  of the one from the Perturbation parent class.

*/

class SpeedyInvariant : public Perturbation {
 protected:
  double Phi,Psi,PhiDot,PsiDot;
  ofstream *ofs,*ofs2, *sourcesFile;
  double denom1[LMX0+1], denom2[LMX0+1];
  double denomE1[LMX0+1], denomE2[LMX0+1];
  double lasttau;
  /*! Photon, Polarization and massless/massive neutrinos
    multipoles. These are just synonyms within
    the y and yprime arrays. Makes code much
    more readable. Set in calcPerturbationNr()
  */
  const double *M, *E, *N, *Nm;
  const double *NR_0, *NR_1, *NR;
  double *Mprime, *Eprime, *Nprime, *Nmprime;
  double *NRprime_0, *NRprime_1, *NRprime;
  double delta_nu_nr_longit, onepwNRVn_nr_longit;
  double rho_nuNR, P_nu_NR;
  double Pi_nu_NR, sigma_nu_nr, del_rho_nu_NR, f_nu_NR, del_P_nu_NR;

 protected:
  double sklm(double s, double l, double m);
  double  s5,s6;  // sqrt(5) and sqrt(6)

  // We do use a different mechanism than the usual perturbation classes 
  // for bookkeeping of 
  // tight coupling switches. 
  // Here are the relevant variables
  map<double, pair<double, Species> >  Thresholds;
  bool CoupledPhoton, CoupledBaryon; //!< Photons and Baryons tightly coupled ?
  bool CoupledOctopole, CoupledMultipole; //!< CoupledOctopole=true: use analytic formula, CoupledMultipole=true: truncate l > 3


  double TauDilute; //!< Tau at which rho_gamma is negligble for certain things
  double Combi; //!< Recursion relation factor: 
  pair<double,double> TauJumpPhoton, TauJumpBaryon;

 public:
  int called,called2,called3;
  SpeedyInvariant(Cosmos* c);
  virtual ~SpeedyInvariant();

  virtual void fderivs(const double tau, const double *y, double *yprime);
  //virtual void fderivsTensor(const double tau,  double *y, double *yprime);
  virtual void initialScalarPerturbations(const ControlPanel &control, const double tau);
  //void initialTensorPerturbations();

  virtual void fillPhiDot(double tau); //!< calculate d/dtau of phi and psi
  virtual void scalarSources(double tau, double *d, double *dp, double *dk);
  //virtual void tensorSources(double tau, double *dt, double *dte, double *dtb);

  virtual void calcPerturbationNr(const ControlPanel &); //! calculate the number of scalar and tensor equations 
  //void foutputt_(const ControlPanel& control,int *n, double *y, double *ypr,  double tau, double *dt, double *dte, double *dtb);
  virtual void setPointers(const double *y, double *yprime); //!< set E,N,M,Eprime ...
  virtual void propagateScalar(double *tau, const double tauend, const double precision);

  virtual void getReady(const ControlPanel&);
  bool isTightCoupling(const double tau, Species); //!< return true, if at tau, we have tight coupling regime
  void insertThreshold(Species);
  double fadeOut(double tau);

 virtual double delta_c() { return y[3];}
 virtual double delta_b() { return y[1];}
 virtual double delta_g() { return 4*y[5];}
 virtual double delta_n() { return y[2*lmaxg + 5];}
 virtual double delta_nr() { return delta_nu_nr_longit; }
 virtual double psi() { return Psi; }
};

#endif 
