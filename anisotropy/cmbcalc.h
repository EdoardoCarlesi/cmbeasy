#ifndef CMBCALC_H
#define CMBCALC_H

#include "cl.h"
#include "safevector.h"
#include "cleanvector.h"
#include "gauge.h"
#include "anchor.h"


#include <string>
#include <map>

class Integrator;
class Cosmos;
class Spline;
class ControlPanel;
class SplineWeb;

/*!
  Central class for all calculations. 
  It starts off the background evolution, overviews the
  perturbation evolution and schedules the integrators to
  perform the final C_l integration 

  Each mode and each multipole are calculated seperately. 
  However, for convenience, a wrapper behaving much 
  the same as CMBFAST's cmbflat() function is also
  defined.

  As there are many Splines and SplineWebs involved, 
  it uses an Anchor to prevent memory leaks.
*/
class CmbCalc  {
 private:
  ofstream *debug;
  int kCountCmb, kCountTransfer; //! number of k - modes to be calculated
  vector<double> mExtraKs;
  SafeVector<double> cmbK, transferK;
  double akmin,akmax; //! minimum k-value needed

  vector<double> mExtraTimeSteps;
  vector<double> mTimeSteps;
  double arel;
  double stpt; // stop tensor: for k*tau > stpt, tensor sources are 0

  map<int, Spline*> bessel;   //!<stores the bessel functions as bessel[multipole l]
  SafeVector<int> j2l; //!< maps sequential numbers 0,1,2,.. to the multipoles given by the bessel file 
  double maxBesselX; //!< the maximum argument for which one of the bessel splines provides a function value 

  int nstep; 
  int n1;
  double ak10, t10;
  Gauge::gauge  gauge;  //! The gauge for the perturbations

  Anchor anchor;   //!< general anchor for things like the webs below, at the end of one calculation, we call kill()
  Anchor BesselAnchor; //!< anchor for the bessel functions. only killed by initjl
  SplineWeb *dtWeb,*dteWeb,*dtbWeb; //!< SplineWebs for temperature, polarization etc. sources
  SplineWeb *dWeb, *dpWeb;

   Integrator* ScalarInt;  //! an integrator class for the final cl - integration. initialized in prepare() and called in oneCl()
   Integrator* TensorInt; //!< an integrator class for the final cl integration, this time for tensors

#ifdef _OPENMP
   omp_lock_t mPowerSplinesLock, mClLock, mTensorSourcesLock, mScalarSourcesLock;
#endif

 public:
   CmbCalc();
   ~CmbCalc();
   //vector<double> cmbflat_(Cosmos* cosmos,string name,ControlPanel&, CL&,const vector<int>& lval);
   void setGauge(Gauge::gauge g) { gauge = g;}  //!< sets the gauge, see also class Gauge
   void cmbflat(Cosmos* cosmos,string name,ControlPanel& control,CL& cl);   //!< convenience wrapper, behaves like cmbfast's  cmblflat
   int go(Cosmos* cosmos,string name,ControlPanel&,bool interactive=false); //!< Start the whole calculation, preparations
   int goTransfer(Cosmos* cosmos,string name,const ControlPanel&,bool interactive=false); //!< prepare for transfer k-modes calculations
   int prepareCl(Cosmos* cosmos,const ControlPanel& control,CL&,bool interactive=false); //!< prepare for Cl - integration
   int prepareScalarCl(Cosmos* cosmos,const ControlPanel& control,CL&,bool interactive=false); //!< prepare scalar for Cl - integration
   int prepareTensorCl(Cosmos* cosmos,const ControlPanel& control,CL&,bool interactive=false); //!< prepare scalar for Cl - integration
   void ridOfScalarIntegrator(); //!< delete the scalar integrator (and with it the vast amount of interpolation splines
   void ridOfTensorIntegrator(); //!< delete the scalar integrator (and with it the vast amount of interpolation splines
   void finish(Cosmos*, const ControlPanel&,CL&); //!< finish up all calculations
   //! compute perturbations for the extraKs, in addition to the standard k values
   void addExtraKValues(const vector<double>& extraKs) { mExtraKs = extraKs; }
   //! compute additional timesteps
   void addExtraTimeSteps(const vector<double>& extraTaus) { mExtraTimeSteps = extraTaus; }
   void oneK(int kNo,Cosmos*,const ControlPanel&); //!< perform the perturbation evolution for one single k-mode
#ifndef PRERELEASE
   //! evolve only this one specific scalar k-mode (for debugging)
   void runOneScalarK(double k,Cosmos* cosmos,const ControlPanel& control, bool checkInitialConditions=false,
                                              std::vector<double>* steps=0, double prec=0);
#endif
   void oneKTransfer(int kNo, Cosmos*,const ControlPanel&); //!< same as oneK() but for higher k-modes
   void oneCl(Cosmos*,int j,CL&,const ControlPanel&); //!< perfrom C_l integration for one multipole l
   void oneTensorCl(Cosmos*,int j,CL&,const ControlPanel&); //!< perfrom C_l for one tensor multipole l

   vector<int> lval(); // !< return vector of all multipole values...

   void initjl(string filename,int MaxMultipole,bool recursive=false);  //!< initialize (i.e. read from disk) the spherical bessel functions 
   static list<int> generateMultipoleValues(int MaxMultipole); //!< generate a list of multipole values at which bessel functions should be pre-calculated
   static void precalculateBesselFunctions( string filename, int MaxMultipole, int keta ); //!< precalculate bessel functions and store in file "jlgen.dat"

   void dumpTransfer(Cosmos *cosmos,ControlPanel &control,map<int,string>&); //!< write transfer spectra to disk
   void dumpCl(const vector<double>& InitialPower, const CL &cl, const ControlPanel &control,
               const string fname, string tname, bool dumpbs = false) const; //!< write cl's to disc
};

#endif 
