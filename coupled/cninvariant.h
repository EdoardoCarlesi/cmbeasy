#ifndef CnINVARIANT_H
#define CnINVARIANT_H

#include "speedydeinvariant.h"
#include "massiveneutrinosnew.h"
#include "perturbationtracker.h"
#include <iostream>
/*!
  
  Gauge invariant perturbation class implementing the speed ups of
  astro-ph/0503277, including massive neutrinos coupled to dark energy

  The dark energy is modelled as a scalar fluid with rest-frame
  speed of sound mu^2 which you can change in fderivs().
*/

class CnCosmos;

class CnInvariant : public SpeedyDEInvariant
{
 public:
  CnInvariant(CnCosmos* c);
  virtual ~CnInvariant() { if (false) cout << "called: " << called << " thereof called2: " << called2 << " called3: " << called3 << endl; }

  virtual void fderivs(const double tau, const double *y, double *yprime);
  virtual void initialScalarPerturbations(const ControlPanel &control, const double tau);

  static void printOSF(string n, double x, double y);
  static void printOSF5(string n, double x, double y, double z, double w, double k);

  virtual void fillPhiDot(double tau); //!< calculate d/dtau of phi and psi


  virtual void calcPerturbationNr(const ControlPanel &); //! calculate the number of scalar and tensor equations 
  //void foutputt_(const ControlPanel& control,int *n, double *y, double *ypr,  double tau, double *dt, double *dte, double *dtb);

  virtual void getReady(const ControlPanel&);

  virtual double psi() const { return Psi; }
  virtual double psiDot() const { return PsiDot; }
  virtual double phi() const { return Phi; }
  virtual double phiDot() const { return PhiDot; }

  virtual double delta_c_longit() const { return mDelta_c_longit; }
  virtual double delta_b_longit() const { return mDelta_b_longit; }
  virtual double v_c_longit() const { return mV_c_longit; }
  virtual double v_b_longit() const { return mV_b_longit; }
  virtual double delta_q_longit() const { return mDelta_q_longit; }
  virtual double delta_nuNr_longit() const { return mDelta_nu_nr_longit; }
  virtual double v_q_longit() const { return mV_b_longit; }
  virtual double v_nu_nr_longit() const { return mV_nu_nr_longit; }

  void stopAllPerturbations(int dim, double *yprime){  // set derivatives to zero i.e. stop growth of the perturbations for deltas beyond a given limit
        for(int j=0; j<dim; j++) { yprime[j]=0; }
  }

  void stopPerturbation(double yprime){ yprime=0.; }

 protected:
  CnCosmos *quintcosmos;
  int qidx; //!< the index position of the delta phi and delta phidot variable in the y[] array
  int xidx; //!< the index of the field perturbation derivatives
  int nuNRidx; //!< the index position of the fluid-perturbation variables for massive neutrinos
  bool IsPhantomCrossing;

 private:
  double mDelta_c_longit, mDelta_b_longit, mV_b_longit, mV_c_longit;
  double mDelta_q_longit, mDelta_nu_nr_longit, mV_q_longit, mV_nu_nr_longit;

 public:
  static double maxGravPotential;
  static double earliestStopZ;
  static double maxEta;
  static double maxDelta; // delta overdensity cutofff
	static double cutoffIndex; // power law index dependence on k scale for the power-law 
 //static int printVariable;

	};

#endif
