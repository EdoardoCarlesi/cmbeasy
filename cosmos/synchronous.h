#ifndef SYNCHRONOUS_H
#define SYNCHRONOUS_H

#include "perturbation.h"

/*!
  Perturbation solver for synchronous gauge variables.
*/

class Synchronous : public Perturbation {
 protected:
  double epsw; //!< fderivs needs this from CMBFAST
  int iq0, iq1, iq2;
  double denl[LMX0+1];
  double qdn[NQMAX0]; //!< For CMBFAST Neutrinoarray only
  double adotoa, hdot, dgshear, rhonu, shearnu, Psi;
 public:
  int called;
  ofstream *ofs;
  Synchronous(Cosmos* c);
  virtual ~Synchronous() { }

  virtual void fderivs(const double tau, const double *y, double *yprime);
  //  virtual void fderivsTensor(const double tau,  double *y, double *yprime);
  virtual void initialScalarPerturbations(const ControlPanel &control, const double tau);
  //void initialTensorPerturbations();

  virtual void scalarSources(double tau, double *d, double *dp, double *dk);
  //virtual void tensorSources(double tau, double *dt, double *dte, double *dtb);

  virtual void calcPerturbationNr(const ControlPanel &); //! calculate the number of scalar and tensor equations 

  virtual void getReady(const ControlPanel&);

 virtual double delta_c() { return y[2];}
 virtual double delta_b() { return y[4];}
 virtual double delta_g() { return y[6];}
 virtual double delta_n() { return y[2*lmaxg + 8];}
 virtual double delta_nr();
 virtual double psi() { return Psi; }
};

#endif
