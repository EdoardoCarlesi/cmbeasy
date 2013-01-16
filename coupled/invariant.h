#ifndef INVARIANT_H
#define INVARIANT_H 

#include "perturbation.h"
#include <fstream>

/*!
  Invariant parent class. Uses Variables 
  M:   Temperature Anisotropy
  E,B:  Polarization Anisotropy

  For their definition, see either

  [] M. Zaldarriaga, U. Seljak et. al, various papers
  [] Ruth Durrer:  astro-ph/0109522,
  and references therein, especially for E and B
  [] my thesis, in which all perturbation equations are
  summarized in an appendix
*/

class Invariant : public Perturbation {
 protected:
  double Phi,Psi,PhiDot,PsiDot;
  ofstream* ofs;
  double denom1[LMX0+1], denom2[LMX0+1];
  double denomE1[LMX0+1], denomE2[LMX0+1];
  /*! Photon, Polarization and massless neutrinos
    multipoles. These are just synonyms within 
    the y and yprime arrays. Makes code much
    more readable. Set in calcPerturbationNr()
  */
  const double *M,*E,*N;
  double *Mprime, *Eprime, *Nprime;
  double epsw; //!< fderivs needs this from CMBFAST
  int  lmaxnr, lmaxnu, iq0, iq1, iq2,lmaxt;

  double qdn[NQMAX0]; //!< For CMBFAST Neutrinoarray only member of common Fermi    
  double hdot, dgshear, rhonu, shearnu;

  double sklm(double s, double l, double m);
  double  s5,s6;  // sqrt(5) and sqrt(6)

 public:
  Invariant(Cosmos* c);
  virtual ~Invariant() {};

  virtual void fderivs(const double tau, const double *y, double *yprime);
  virtual void fderivsTensor(const double tau,  double *y, double *yprime);
  virtual void initialScalarPerturbations(const ControlPanel &control, const double tau);
  void initialTensorPerturbations();

  virtual void fillPhiDot(double tau); //!< calculate d/dtau of phi and psi
  virtual void scalarSources(double tau, double *d, double *dp, double *dk);
  virtual void tensorSources(double tau, double *dt, double *dte, double *dtb);
 
  virtual void calcPerturbationNr(const ControlPanel &); //! calculate the number of scalar and tensor equations 
  //void foutputt_(const ControlPanel& control,int *n, double *y, double *ypr,  double tau, double *dt, double *dte, double *dtb);
  virtual void setPointers(const double *y, double *yprime); //!< set E,N,M,Eprime ...

  
  void nu2(const double a, double *drhonu, double *fnu, double *dpnu, double *shearnu, const double *psi0, const double * psi1, const double *psi2);
  void nuder(const double, const double, const double, double*, double*, const double*, const double*);
  virtual void getReady(const ControlPanel&);


 virtual double delta_c() { return y[3];}
 virtual double delta_b() { return y[1];}
 virtual double delta_g() { return 4*y[5];}
 virtual double delta_n() { return y[2*lmaxg + 5];}
 virtual double delta_nr(double tau);
 virtual double psi() { return Psi; }
};

#endif 
