#include "global.h"

#include "spline.h"
#include "cosmos.h"
#include "cleanvector.h"
#include "tensortrapezoidintegrator.h"
#include "scalartrapezoidintegrator.h"
#include "controlpanel.h"
#include "cmbcalc.h"
#include "perturbationfactory.h"
#include "perturbation.h"
#include "model.h"
#include "bessel.h"
#include "minmax.h"

#include <string>
#include <cmath>
#include <fstream>
#include <iterator>
#include <algorithm>

//#define OMP_DEBUG

/*!
  Easy wrapper for all the parts of the cmb calculation.
  This subsequently calls all functions needed to calculate the
  Cl's and the transfer functions. 
*/
void CmbCalc::cmbflat(Cosmos* cosmos,string name,ControlPanel& control,CL& cl) {
  // brief plausibility check
  if (control.cmb && ! (control.scalar || control.tensor)) throw Bad_Error("CmbCalc::cmbflat() cmb but neither scalar or tensor wanted");
  if ((control.power_cdm) && !control.scalar) throw Bad_Error("CmbCalc::cmbflat() cdm power but no scalars wanted");
  // if no bessel functions are in, we throw bad error
  if (control.cmb && bessel.empty()) throw Bad_Error("CmbCalc::cmbflat() No Bessel functions read in.\nCall initjl() first");   
  
  go(cosmos,name,control);
 
  if (control.power_cdm) goTransfer(cosmos,name,control);
  if (! (control.cmb || control.power_cdm)) return;  // If go has been used to set up cosmos only, go back 
  if (control.minK != 0) return; // Can't do cmb with too high of a kmin
  prepareScalarCl(cosmos,control,cl);
  prepareTensorCl(cosmos,control,cl);
  finish(cosmos,control,cl);
}


#ifdef _OPENMP
struct OnScopeExit
{
  omp_lock_t* mLock;
  OnScopeExit(omp_lock_t* l) : mLock(l) {}
  ~OnScopeExit() { omp_unset_lock(mLock); }
};
#endif

typedef double (Perturbation::*ValueOfPerturbation)();
#define CALL_MEMBER_METHOD(object, ptrToMember)  ((object).*(ptrToMember))

struct PowerSplines
{
  typedef map<SplineWeb*, ValueOfPerturbation> WebToPertMap;
  WebToPertMap mSpectra;
#ifdef _OPENMP
  omp_lock_t* mPowerSplinesLock;
  void setLock(omp_lock_t* l) { mPowerSplinesLock=l; }
  omp_lock_t* lock() { omp_set_lock(mPowerSplinesLock); return mPowerSplinesLock; }
#endif

  PowerSplines(Cosmos& c) {
   addSpectrum( c.power_cdm(),    &Perturbation::delta_c);
   addSpectrum( c.power_baryon(), &Perturbation::delta_b);
   addSpectrum( c.power_gamma(),  &Perturbation::delta_g);
   addSpectrum( c.power_nu(),     &Perturbation::delta_n);
   addSpectrum( c.power_psi(),    &Perturbation::psi);
   addSpectrum( c.power_nuNR(),   &Perturbation::delta_nr);
  }

  void addSpectrum(SplineWeb* w, ValueOfPerturbation v) {
    mSpectra.insert(make_pair(w, v));
  }

  typedef map<SplineWeb* /*PowerWeb*/, map<double /*tau*/, double /*pertval*/> > PowerMaps;
  PowerMaps maps(Cosmos& c) const {
    PowerMaps retVal;
    int idx = 0;
    WebToPertMap::const_iterator it, end=mSpectra.end();
    for (it = mSpectra.begin(); it!=end; ++it) {
     SplineWeb *power = it->first;
     retVal.insert(make_pair(power, map<double, double>()));
    }
    return retVal;
  }

  void fillPowerMaps(PowerMaps& m, Perturbation& p, double tau) {
    PowerMaps::const_iterator it, end=m.end();
     for (it = m.begin(); it!=end; ++it) {
        double pertVal = CALL_MEMBER_METHOD(p, mSpectra[it->first])();
        m.at(it->first)[tau]=pertVal;
     }
  }


  void powerMapsToSplineWebs(Cosmos& c, PowerMaps& m, double k) {
#ifdef _OPENMP
    lock();
    OnScopeExit releaseLock(mPowerSplinesLock);
#endif
    PowerMaps::const_iterator it, end=m.end();
    map<double, double>::const_iterator vit, vend;
    for (it = m.begin(); it!=end; ++it) {
      SplineWeb *power = it->first;
      const PowerMaps::mapped_type& pertVals = it->second;
      vend = pertVals.end();
      if (it->first == c.power_psi()) {
        for(vit=pertVals.begin(); vit != vend; ++vit)
          c.fillPsiWeb(power, k, vit->first /*tau*/, vit->second /*pertVal*/);
      } else if ((it->first == c.power_nuNR()) && (c.nuNR() != 0)) {
        for(vit=pertVals.begin(); vit != vend; ++vit)
          c.fillPower(power, k, vit->first /*tau*/, vit->second /*pertVal*/);
      } else {
        for(vit=pertVals.begin(); vit != vend; ++vit)
          c.fillPower(power, k, vit->first /*tau*/, vit->second /*pertVal*/);
      }
    }
  }
};

/*!
  Constructor for cmbcalc. Please note that maxBesselX is initialized such that 
  cmbcalc can be "misused" to get the background evolution and thermal history only.
  Usually, when pertubations are requested maxBesselX will be initialized properly
  in initjl()
*/
CmbCalc::CmbCalc()
      : arel(1), maxBesselX(3e3), nstep(0),
        gauge(Gauge::synchronous), ScalarInt(0), TensorInt(0)
{
  // if you like to have some output stream for debugging,
  // uncomment the line below.
  // debug = new ofstream("debug.dat");  

#ifdef _OPENMP
   omp_init_lock(&mPowerSplinesLock);
   omp_init_lock(&mTensorSourcesLock);
   omp_init_lock(&mScalarSourcesLock);
   omp_init_lock(&mClLock);
#endif
}

CmbCalc::~CmbCalc()
{
#ifdef _OPENMP
   omp_destroy_lock(&mPowerSplinesLock);
   omp_destroy_lock(&mTensorSourcesLock);
   omp_destroy_lock(&mScalarSourcesLock);
   omp_destroy_lock(&mClLock);
#endif
}

/*! 
  The very first step in calculating the power spectrum and CMB anisotropies.
  As always, all units are in Mpc. Prominent players:
  k:  comoving Wave vector this k, *not* k/h. 
  (The cold dark matter powerspectrum is written to the Cosmos powerspectrum 
  as a function of  k/h).
  
  If interactive is true, it will not call oneK() itself, but return immediately. It is up
  to the GUI to call oneK(), then.
*/

int CmbCalc::go(Cosmos* cosmos,string name,ControlPanel& control, bool interactive) {
  double dlnk;
  double zst;  // zst is the value of z where the program stops the calculation.
  if ( fabs(cosmos->omega_m() - 1) < .001 && cosmos->reionizationZ() == 0. && cosmos->optDistanceLss()
       == 0. && (!control.tensor) && cosmos->h() > 0.40)  zst = 10.; else  zst = 0.;

  if (control.power_cdm) zst = min(zst,  control.transferZ.back() );

#ifdef OMP_DEBUG
  double time_before = omp_get_wtime();
#endif
  cosmos->getReady();
#ifdef OMP_DEBUG
  double time_after = omp_get_wtime();
  cout << "get ready: " << (time_after-time_before) << endl;
  time_before = omp_get_wtime();
#endif
  cosmos->history(control.verbose);   // Calculate history of the cosmos, true means: print information, dump several splines to files
#ifdef OMP_DEBUG
  time_after = omp_get_wtime();
  cout << "history: " << (time_after-time_before) << endl;
#endif

  // Initialize massive neutrinos.
  if (cosmos->nuNR() != 0. || cosmos->omega_nuNR() != 0.0) {
    arel = .001 / cosmos->amnu;
    MassiveNeutrinos::initnu1(cosmos->amnu);
  }
  // Calculating when reinization starts.
  if (cosmos->optDistanceLss() > 0) cosmos->reiopar();

  //    Maximum and minimum k-values. dlnk is the logarithmic spacing 
  //     used for low k. 
  akmax = maxBesselX / cosmos->tau_0();  // as x = k*tau -> k_max = x_max / tau_0
  akmin = .15 / cosmos->tau_0();
  if (control.minK !=0 ) {
    akmin = control.minK;
  }
  if (! control.tensor) dlnk = .1; else dlnk = .05;


  // Thermal history and timesteps etc. are now calculated. nstep is set
  double MaxK = akmax;
  if (control.power_cdm) MaxK = max(MaxK, control.transferMaxK);
#ifdef OMP_DEBUG
  time_before = omp_get_wtime();
#endif
  cosmos->thermo(MaxK*2,control);
#ifdef OMP_DEBUG
  time_after = omp_get_wtime();
  cout << "thermo: " << (time_after-time_before) << endl;
#endif
  //  cout << "Got timestpes " << endl;
  mTimeSteps.resize(0);
  using namespace std;
  merge(cosmos->TimeStep.begin(), cosmos->TimeStep.end(), mExtraTimeSteps.begin(), mExtraTimeSteps.end(),
        back_inserter(mTimeSteps));
  unique(mTimeSteps.begin(), mTimeSteps.end());

  nstep = mTimeSteps.size();

  //cosmos->printStatus(("cosmos->status"));  // if you want the status of cosmos in a file, uncomment this

  // stop doing anything if nothing is wanted
  // this is useful if you just want to use finithermo conveniently
  // but not calculate perturbations
  if (! (control.cmb || control.power_cdm)) return 0;

  //    Calculation of the CMB sources.
  //     set k values for which the sources for the anisotropy and
  //     polarization will be calculated. For low values of k we
  //     use a logarithmic spacing.

  double dlnk0;
  if (cosmos->reionizationZ() )  dlnk0 = dlnk * 2.; else dlnk0 = dlnk * 5.;

  int nk, nk2;
  if (control.cmb) {
    double dkn1 = .8 / cosmos->BeginRecombination;
    double dkn2 = 1.5 / cosmos->BeginRecombination;
    double nk1 = (int) (log(dkn1 / akmin / dlnk0) / dlnk0) + 1;
    if (akmax > 1500. / cosmos->tau_0()) {
      nk2 = (int) ((1500. / cosmos->tau_0() - akmin * exp((nk1 - 1) * dlnk0)) / dkn1) + (int)nk1 + 1;
      nk = (int) ((akmax - 1500. / cosmos->tau_0()) / dkn2) + nk2 + 1;
    } else {
      nk = (int) ((akmax - akmin * exp((nk1 - 1) * dlnk0)) / dkn1)  + (int)nk1 + 1;
      nk2 = nk + 1;
    }
    cmbK.resize(nk);
    for (int k = 1; k <= nk; ++k) {
      if (k <= nk1) {
        cmbK[k-1] = akmin * exp((k - 1) * dlnk0);
      } else {
        if (k > nk2) {
          cmbK[k-1] = cmbK[nk2-1] + (k - nk2) * dkn2;
        } else {
          cmbK[k-1] = cmbK[(int)nk1-1] + (k - nk1) * dkn1;
        }
      }
    }
  }

  if (!mExtraKs.empty()) {
    double lastCmbK = cmbK.empty()?0:cmbK.back();
    sort(mExtraKs.begin(), mExtraKs.end());
    vector<double> tempCmbK;
    merge(cmbK.begin(), cmbK.end(), mExtraKs.begin(), mExtraKs.end(), back_inserter(tempCmbK));
    vector<double>::iterator upper = upper_bound(tempCmbK.begin(), tempCmbK.end(), lastCmbK);
    cmbK = SafeVector<double>(tempCmbK.begin(), upper);
    nk = cmbK.size();
  }

  /*
  ofstream ts("timesteps2.dat");
  for (int i =0; i < nstep; i++) ts << cosmos->TimeStep[i] << "  " << cosmos->Visibility->fastY(cosmos->TimeStep[i]) << endl;
  */


  if (control.cmb) {
    int sourcesize = nstep; // counting is a hard task :-) 
    // now we initialize the SplineWeb's for the sources
    if (control.tensor) {
      dtWeb = new SplineWeb("dtWeb", &anchor, nk, sourcesize); 
      dteWeb = new SplineWeb("dteWeb", &anchor, nk, sourcesize);
      dtbWeb = new SplineWeb("dtbWeb", &anchor, nk, sourcesize);
    }
    if (control.scalar) {
      dWeb = new SplineWeb("dWeb", &anchor, nk, sourcesize);
      dpWeb = new SplineWeb("dpWeb", &anchor, nk, sourcesize);
    }
  }
  // if no cmb is calculated, there is only the power spectrum left
  // we leave the whole k-integration to goTransfer, setting our
  // part in it to 0 ...
  if (! control.cmb) nk =0;
  if (interactive) return nk;

  // Loop over wavenumbers, calculating perturbations for each k-mode
#ifdef OMP_DEBUG
  time_before = omp_get_wtime();
#endif
#ifdef WITH_OMP
#pragma omp parallel for
#else
  for (int kNo = 0; kNo < nk; ++kNo)  {
    try {
      oneK(kNo, cosmos, control);
    } catch (Bad_Error e) {
#ifdef _OPENMP
      cout << "caught Bad_Error from CmbCalc::oneK()\n" << e.s << endl;
#endif
      throw;
    }
  }
#ifdef OMP_DEBUG
  cout << "wtime onek: " << (omp_get_wtime()-time_before) << endl;
#endif
  std::vector<SplineWeb*> allPowerWebs = cosmos->powerSplineWebs();
  int size = allPowerWebs.size();
#ifdef WITH_OMP
#pragma omp parallel for
#endif
  for (int i = 0; i<size; ++i) {
    allPowerWebs[i]->makeProper();
  }
#ifdef OMP_DEBUG
  cout << "wtime onK+fixing: " << (omp_get_wtime()-time_before) << endl;
#endif
  return 0;
}



/*!
  Essentially the same as go() but for the k-modes that are
  needed in addition for the power-spectrum and transfer functions (and
  sigma8 
*/
int CmbCalc::goTransfer(Cosmos* cosmos,string name,const ControlPanel& control,bool interactive) {
  // If transfer functions are requested, set the remaining k-values.
  int nkt, nk=cmbK.size();
  double akdone=0;
  if (control.power_cdm) {
    if (control.scalar && control.cmb) { // if there has already been some work done in oneK() on scalars
      if (control.transferMaxK > akmax) {
        nkt = nk + (int) ((log(control.transferMaxK) - log(akmax)) * control.transferPerLog) + 1;
        akdone = cmbK.back();
      } else nkt = nk;
    } else {  // so it is our task to calculate everything
      nkt = (int) ((log(control.transferMaxK) - log(akmin)) * control.transferPerLog)  + 1;
      akdone = akmin;
    }

    transferK.resize(0);
    for (int k = nk + 1; k <= nkt; ++k) {
      double kVal = (akdone * exp((double) (k - nk) / control.transferPerLog));
      transferK.push_back(kVal);
    }

    if (!mExtraKs.empty()) {
      double firstTransferK = transferK.front();
      sort(mExtraKs.begin(), mExtraKs.end());
      vector<double> tempTransferK;
      merge(transferK.begin(), transferK.end(), mExtraKs.begin(), mExtraKs.end(), back_inserter(tempTransferK));
      vector<double>::iterator lower = lower_bound(tempTransferK.begin(), tempTransferK.end(), firstTransferK);
      transferK = SafeVector<double>(lower, tempTransferK.end());
    }
    nkt = transferK.size()+nk;

    if (interactive) return nkt;
    //   Loop over wavenumbers.
#ifdef WITH_OMP
#pragma omp parallel for
#endif
    for (unsigned int kNo = 0; kNo < transferK.size(); ++kNo) {
      try {
        oneKTransfer(kNo, cosmos,control);
      } catch (Bad_Error e) {
#ifdef _OPENMP
        cout << "caught Bad_Error from CmbCalc::oneKTransfer\n" << e.s << endl;
#endif
        throw;
      }
    }
    std::vector<SplineWeb*> allPowerWebs = cosmos->powerSplineWebs();
    for_each(allPowerWebs.begin(), allPowerWebs.end(),
             mem_fun(&SplineWeb::makeProper));
  }
  return 0;
}

/*!
  This prepares for scalar cl-integration, i.e. it allocates the cl.ts etc splines
  and creates a ScalarIntegrator. Then it integrates, if it isn't in interactive
  mode. CAUTION: The creation of 
  the ScalarIntegrator consumes quite some memory. 
  CAUTION 2: Call this *** even if you don't want the scalar cl's ***. Or even
  better, if you can afford the expense of having tensor and scalar memory
  allocated at once: call prepareCl(). The reason is, because in prepareScalarCl(),
  the cl.ts spline family is generated. The scalar integrator on the other hand
  will only be created if scalar cl's are requested. In that way, prepareTensorCl() 
  will have the cl.ts mother spline up and running, by the time it is called.

  Even in the case that you don't want any cmb spectrum at all, go ahead
  and call this, it will ask the cdm power web to shrinkToFit(). Very important!
  
*/
int CmbCalc::prepareScalarCl(Cosmos* cosmos, const ControlPanel& control,CL &cl,bool interactive) {
  cl.clear();
  if (control.cmb) {
    Spline::generateFamily(cl.ts,cosmos->InitialPower.size(), bessel.size(), "clts");
    Spline::generateChildren(cl.es,cl.ts[0], cosmos->InitialPower.size(), "cles");
    Spline::generateChildren(cl.cs,cl.ts[0], cosmos->InitialPower.size(), "clcs");
    Spline::generateChildren(cl.bs,cl.ts[0], cosmos->InitialPower.size(), "clbs");
    Spline::generateChildren(cl.kk,cl.ts[0], cosmos->InitialPower.size(), "clkk");
    Spline::generateChildren(cl.tk,cl.ts[0], cosmos->InitialPower.size(), "cltk");
  }

  if (control.power_cdm)
    cosmos->power_cdm()->shrinkToFit();

  if (control.scalar && control.cmb) {
    //cout << endl;
    dWeb->makeProper();
    dpWeb->makeProper();
    //    ScalarInt = new ScalarIntegrator(cosmos, akmin, akmax, dWeb, dpWeb);
     ScalarInt = new ScalarTrapezoidIntegrator(cosmos, akmin, akmax, dWeb, dpWeb);
      if (interactive) return bessel.size();  // if interactive, we don't call oneCl etc ourselves
#ifdef OMP_DEBUG
   double time_before = omp_get_wtime();
#endif
#ifdef WITH_OMP
#pragma omp parallel for
#endif
      for (unsigned int j = 0; j < bessel.size() ; ++j)
          oneCl(cosmos, j, cl, control);  // calculate scalar cl's
      cl.ts[0]->makeProper();

#ifdef OMP_DEBUG
    cout << "wtime cl: " << (omp_get_wtime()-time_before) << endl;
#endif
    ridOfScalarIntegrator();
  }
  return -1;
}

/*!
  Same as prepareScalarCl(), but for tensors.
*/
int CmbCalc::prepareTensorCl(Cosmos* cosmos, const ControlPanel& control,CL &cl,bool interactive) {
  if (control.cmb) {
    Spline::generateChildren(cl.tt,cl.ts[0], cosmos->InitialPower.size(), "cltt");
    Spline::generateChildren(cl.et,cl.ts[0], cosmos->InitialPower.size(), "clet");
    Spline::generateChildren(cl.bt,cl.ts[0], cosmos->InitialPower.size(), "clbt");
    Spline::generateChildren(cl.ct,cl.ts[0], cosmos->InitialPower.size(), "clct");
  }

  if (control.tensor && control.cmb) {
    cout << endl;
    //    TensorInt = new TensorIntegrator(cosmos, akmin,akmax,stpt,dtWeb,dteWeb,dtbWeb);
    dtWeb->makeProper();
    dteWeb->makeProper();
    dtbWeb->makeProper();
    TensorInt = new TensorTrapezoidIntegrator(cosmos, akmin,akmax,stpt,dtWeb,dteWeb,dtbWeb);
    if (interactive) return bessel.size();  // if interactive, we don't call oneCl etc ourselves
    for (unsigned int j = 0; j < bessel.size() ; ++j) oneTensorCl(cosmos,j,cl,control); // calculate all cl's 
    ridOfTensorIntegrator();
  }
  return -1;
}

/*!
  All-In-One Wrapper that prepares (and as both of them integrate also,
  integrates to produce the C_l's)
*/
int  CmbCalc::prepareCl(Cosmos* cosmos, const ControlPanel& control,CL &cl,bool interactive) {
  if (control.cmb) {
    prepareScalarCl(cosmos,control,cl,interactive);
    prepareTensorCl(cosmos,control,cl,interactive);
    if (interactive) return bessel.size();
  }  // if this function has no sense anyhow 
  return -1;// meaning :: do not calculate any cl for a interactive program ...
}

/*!
  Clean up the tensor integrator
*/
void CmbCalc::ridOfTensorIntegrator() {
  if (TensorInt) {
    delete TensorInt; TensorInt=0;
  }
}
/*!
  Clean up the scalar integrator
*/
void CmbCalc::ridOfScalarIntegrator() {
  if (ScalarInt) {
    delete ScalarInt; ScalarInt=0;
  }
}

/*!
  The last step in the cmb calculation, frees memory by asking anchor to kill
  For safety purposes (not that you want to forget to do it before :-) also, preparexxxxCl()
  does it for you, it calls ridOfScalarIntegrator() etc.
*/
void CmbCalc::finish(Cosmos* cosmos,const ControlPanel& control,CL& cl) {
  ridOfScalarIntegrator();
  ridOfTensorIntegrator();
  if (control.cmb) {
    //     Final calculations for CMB output.
    //     Scalar case
    if (! control.scalar) {
      // if the scalar spectrum has not been calculated, we have a problem,
      // cause cl.ts[0] holds the x - values for all return splines and it has never
      // been set. So, if control.scalar is not set, we fill the cl.ts[0] spline with (l , 0)
      // values      
      for (map<int,Spline*>::iterator i=bessel.begin(); i != bessel.end(); i++) 
	cl.ts[0]->setForce(i->first,0); 
    }
  }
  anchor.kill();  // this kills all splinewebs and frees the memory
} 



// ############ ALTERNATIVE ONE K
void CmbCalc::oneK(int k,Cosmos* cosmos,const ControlPanel& control) {
  Perturbation *pert = PerturbationFactory::perturbation(gauge,cosmos);
  pert->k=cmbK[k];

  if (k>1) {
    pert->setWarnOnEarlySources(false);
  }

  double taustart0;
  double precision = 1e-4; 
  if (cosmos->omegaH2_b() < 0.01) precision = 1e-5;

  pert->getReady(control);


   double taustart = cosmos->tauStartforKMode(pert->k);
   taustart0 = taustart;
   //     Start when massive neutrinos are strongly relativistic. 
   if (cosmos->haveMassiveNeutrinos()) {
     arel = .001 / cosmos->amnu;

     taustart = min(taustart0, arel / cosmos->tau2adot(taustart));
     double taustartLate = max(taustart, cosmos->Tau2A->start());
     if (taustartLate>taustart) {
//X        cout << "CmbCalc::oneK() - using late taustart" << endl;
       taustart=taustartLate;
     }
     //cout << "taustart: " << taustart<< " z=" << cosmos->tau2z(taustart) << endl;

   }
   //cout << ".";   // indicate progress
   //cout.flush();


   // if you are interested in some special k-mode
   // then maybe for convenience, you would like to interrupt
   // the whole thing at some k-value. 
   // if so, give ControlPanel.StopPert some k-value at which to stop

   //   taustart *= 10;

   //   cout << "taustart: " << taustart << " xstart: " << taustart*pert->k << "opac: " << cosmos->opac(taustart) <<  endl;

   if (control.scalar || control.power_cdm) {
     cout.flush();
     pert->initialScalarPerturbations(control,taustart);
     double tau = taustart;

     //  Begin timestep loop. 
     //
     //   d contains the sources for the anisotropy and dp 
     //   for the polarization. t means tensor. 
     // Unfortunately, taustart may not be the same for all modes
     // so we can't keep the very first (deep RD) initial perturbation here.

     int TimeSteps = mTimeSteps.size();
     int n = 0;
     while (mTimeSteps[n]<tau) ++n;
     int firstn = n;
     PowerSplines powerSplineSaver(*cosmos);
#ifdef _OPENMP
     powerSplineSaver.setLock(&mPowerSplinesLock);
#endif
     PowerSplines::PowerMaps pMaps = powerSplineSaver.maps(*cosmos);
     map<double, double> dValues, dpValues;
     for (; n < TimeSteps; n ++) {
       double tauend = mTimeSteps[n];
       pert->propagateScalar(&tau,tauend,precision);
       if (control.scalar && control.cmb) {
         if (n == TimeSteps-1) { // so for all but the last
           dValues[tauend]=0.;
           dpValues[tauend]=0.;
         } else {
           double t,e,te;
           pert->scalarSources(tau, &t, &e, &te);   // get the scalar sources
           dValues[tau]=t;
           dpValues[tau]=e;
         }
       }
       if (control.power_cdm){
         powerSplineSaver.fillPowerMaps(pMaps, *pert, tauend);
       }
     }
#ifdef _OPENMP
     {
       omp_set_lock(&mScalarSourcesLock);
       OnScopeExit releaseLock(&mScalarSourcesLock);
#endif
       for (unsigned int n=firstn; n < TimeSteps; n ++) {
         double tau = mTimeSteps[n];
         dWeb->set(pert->k, tau, dValues[tau]);
         dpWeb->set(pert->k,tau, dpValues[tau]);
       }
#ifdef _OPENMP
       omp_unset_lock(&mScalarSourcesLock);
     }
#endif
     if (control.power_cdm)
       powerSplineSaver.powerMapsToSplineWebs(*cosmos, pMaps, pert->k);
   }

   if (pert->k > control.StopPert) throw Bad_Error("k-breakpoint reached");

   //   time loop for tensors if requested. 	    
   if (control.tensor) {
     pert->initialTensorPerturbations();
     double tau = taustart0;

     //   Begin timestep loop. 
     //   dt contains the sources for the anisotropy and dtp
     //   for the polarization. t means tensor. 
     stpt = 50;
     int TimeSteps = mTimeSteps.size();
     int j = 0;
     while (mTimeSteps[j]<tau) ++j;
     int firstj = j;
     map<double, double> dtValues, dteValues, dtbValues;
     for ( ; j < TimeSteps; ++j) {
       double tauend = mTimeSteps[j];
       //   Tensors will only be calculated if k*tau< stsp. 
       if (pert->k * tauend > stpt || j == TimeSteps-1) {
         dtValues[tauend] = 0;
         dteValues[tauend] = 0;
         dtbValues[tauend] = 0;
       } else {
         pert->propagateTensor(&tau,  tauend, precision);
         double dt,dte,dtb;
         pert->tensorSources(tau, &dt, &dte,&dtb);
         dtValues[tauend] = dt;
         dteValues[tauend] = dte;
         dtbValues[tauend] = dtb;
       }
     }
#ifdef WITH_OMP
#pragma omp critical (TensorSources)
#endif
     for (j=firstj; j < TimeSteps; ++j) {
       double tauend = mTimeSteps[j];
       dtWeb->set(pert->k, tauend, dtValues[tauend]);
       dteWeb->set(pert->k, tauend, dteValues[tauend]);
       dtbWeb->set(pert->k, tauend, dtbValues[tauend]);
     }
   }
   delete pert; pert=0;
}

#ifndef PRERELEASE
void CmbCalc::runOneScalarK(double k,Cosmos* cosmos,const ControlPanel& control,
                            bool checkInitialConditions, std::vector<double>* steps, double prec) {
  Perturbation *pert = PerturbationFactory::perturbation(gauge,cosmos);

  double precision = 1e-4; 
  if (cosmos->omegaH2_b() < 0.01) {
    precision = 1e-5;
  }
  if (prec!=0) {
    precision=prec;
  }
  pert->k = k;
  pert->getReady(control);


   //  Begin timestep loop. 
   //
   //   d contains the sources for the anisotropy and dp 
   //   for the polarization. t means tensor. 
   // Unfortunately, taustart may not be the same for all modes
   // so we can't keep the very first (deep RD) initial perturbation here.

   vector<double> tauSteps;
   if (steps) {
     using namespace std;
     merge(cosmos->TimeStep.begin(), cosmos->TimeStep.end(), steps->begin(), steps->end(),
                                                                   back_inserter(tauSteps));
     unique(tauSteps.begin(), tauSteps.end());
   } else {
     tauSteps = cosmos->TimeStep;
   }

   double taustart = .001 / pert->k;
   taustart = min(tauSteps.front(), taustart);
   taustart = min(taustart, 1e-3*cosmos->tau_equ());

   cout << ".";   // indicate progress
   cout.flush();

   cout.flush();
   pert->initialScalarPerturbations(control,taustart);

   Perturbation *pert2=0;
   if (checkInitialConditions) {
     pert2 = PerturbationFactory::perturbation(gauge,cosmos);
     pert2->k = k;
     pert2->getReady(control);
     pert2->initialScalarPerturbations(control,taustart);
   }

   double tau = taustart;

   int timeSteps = tauSteps.size();
   for (int n = 0; n < timeSteps; n ++) {
     double tauend = tauSteps[n];
     if (tauend < taustart)
       continue;
     cout << "propagateScalar: from " << tau << " to " << tauend << "     \r";
     pert->propagateScalar(&tau,tauend,precision);
     if (pert2) {
       pert2->initialScalarPerturbations(control,tau);
     }
     double t,e,te;
     pert->scalarSources(tau, &t, &e  , &te );
     if (pert2)
       pert2->initialScalarPerturbations(control,tau);
   }
   cout << endl;

   delete pert; pert=0;
   delete pert2; pert2=0;
}
#endif // ndef PRERELEASE

void CmbCalc::oneKTransfer(int kNo, Cosmos* cosmos,const ControlPanel& control) {
  double precision = 1e-4;
  Perturbation *pert = PerturbationFactory::perturbation(gauge,cosmos);
  pert->k = transferK[kNo];

  if (cosmos->omegaH2_b() < 0.01) precision = 1e-5;
  //cout << "_"; cout.flush();

  if (kNo>1) {
    pert->setWarnOnEarlySources(false);
  }


//X   cout << " doing = " << pert->k << endl;
//X   return;
  pert->getReady(control);

  //    Begin when wave is far outside horizon. 
  //    Conformal time (in Mpc) in the radiation era, for photons plus 3 species
  //     of relativistic neutrinos. 
  double taustart = .001 / pert->k;
  //    Make sure to start early in the radiation era. 
  taustart = min(taustart,.1);
  //     Start when massive neutrinos are strongly relativistic. 
  if (cosmos->haveMassiveNeutrinos()) {
    arel = .001 / cosmos->amnu;
    taustart = min(taustart,  arel / cosmos->tau2adot(taustart));
  }
  pert->initialScalarPerturbations(control,taustart);

  double tau = taustart;

  //cout << "taustart: " << taustart<< " z=" << cosmos->tau2z(taustart) << endl;
  int TimeSteps = mTimeSteps.size();
  int n = 0;
  PowerSplines powerSplineSaver(*cosmos);
#ifdef _OPENMP
  powerSplineSaver.setLock(&mPowerSplinesLock);
#endif
  PowerSplines::PowerMaps pMaps = powerSplineSaver.maps(*cosmos);
  while (mTimeSteps[n]<tau) ++n;
  for (; n < TimeSteps; ++n) {
    // cout << " before prpscalar tau: " << tau << endl;
    double tauend = mTimeSteps[n];
    pert->propagateScalar(&tau, tauend, precision);
    powerSplineSaver.fillPowerMaps(pMaps, *pert, tauend);
  }
  if (control.power_cdm)
    powerSplineSaver.powerMapsToSplineWebs(*cosmos, pMaps, pert->k);
  if (pert->k > control.StopPert) throw Bad_Error("k-breakpoint reached");  
  delete pert; pert =0;
}

#define FOURPI (4.0*M_PI)
#define FOURPI2 (FOURPI*FOURPI)
/*!
  Perform folding of bessel function with scalar source to get C_l's.
  At the end of this routine, the Cl's are in fact l*(l+1)*C_l.
  In AnalyzeThis::quickWMAPNormalize(), the spectrum will then
  be normalized to fit the data with an additional factor of 2*Pi.
*/
void CmbCalc::oneCl(Cosmos* cosmos, int j,CL& cl,const ControlPanel& control) {
  //cout << "*"; cout.flush();
  double l=j2l[j];  // needs to be double, cause product ctnorm blows int bounds  
  // do the integration this is the most costy thing in the program :-)
  ClReturn clReturn = ScalarInt->integrate((int)l,bessel[(int) l],control);
  double ctnorm = l * (l - 1) * (l + 1) * (l + 2);     // conventional l factors 
  double l2 =  FOURPI2*l * (l+ 1);  // (4*pi)^2 * l*(l+1)
  // loop over all spectral indices 
  // The FOURPI2 is needed from the integration to get the C_l's
#ifdef WITH_OMP
#pragma omp critical
#endif
  for (unsigned int n = 0; n < cosmos->InitialPower.size(); n++) {
    cl.ts[n]->setForce(l,  clReturn.clts[n] * l2);      
    cl.es[n]->set(clReturn.cles[n]* ctnorm * l2);
    cl.cs[n]->set(clReturn.clcs[n]* sqrt(ctnorm) * l2);
    cl.kk[n]->set(0);
    cl.tk[n]->set(0);
  }
}

void CmbCalc::oneTensorCl(Cosmos* cosmos, int j,CL& cl, const ControlPanel &control) {
  //cout << "#"; cout.flush();
  int l=j2l[j];
  // do the integration this is the most costy thing in the program :-)
  ClReturn clReturn = TensorInt->integrate(l, bessel[l],control);
   double ctnorm = l * (l - 1.) * (l + 1) * (l + 2);
   double l2 = l * (l + 1)*FOURPI2;
   // loop over all spectral indices
   for (unsigned int n = 0; n < cosmos->InitialPower.size(); n++) {		   
     cl.tt[n]->setForce(l,ctnorm*l2*clReturn.cltt[n]);
     cl.et[n]->setForce(l,l2*clReturn.clet[n]);
     cl.bt[n]->setForce(l,l2*clReturn.clbt[n]);
     cl.ct[n]->setForce(l,sqrt(ctnorm)*l2*clReturn.clct[n]);
   } 
}


void CmbCalc::dumpTransfer(Cosmos *cosmos, ControlPanel &control, map<int,string>& transferFile) {
  Anchor anker;
  for (unsigned int i = 0; i < control.transferZ.size(); i++) {
    ofstream t(transferFile[i].c_str());
    t.setf(ios::scientific);
    double tau = cosmos->z2tau(control.transferZ[i]);
    Spline *c = cosmos->power_cdm()->createAlongX(tau,"cdm",&anker);
    Spline *b = cosmos->power_baryon()->createAlongX(tau,"baryon",&anker);
    Spline *g = cosmos->power_gamma()->createAlongX(tau,"photon",&anker);
    Spline *nu =  cosmos->power_nu()->createAlongX(tau,"neutrinos",&anker);

    // from v4.0 on, we correct for an additional factor of k
    // that is no longer stored in the power webs (see also scalarSpectrum() 
    // and fillPower() in cosmos )
    for (int i = 0; i < c->size(); i++) {
      double kh = c->x(i); // the power splines store k/h as argument
      double k = kh*cosmos->h();
      double tc = sqrt(exp(c->y(i))/k);
      double tb = sqrt(exp(b->y(i))/k);
      double tg = sqrt(exp(g->y(i))/k);
      double tn = sqrt(exp(nu->y(i))/k);
      if (cosmos->nuNR() != 0) {
        ; // should write massnu here
      }
        t <<  kh << " " << tc << " " << tb << " " << tg  << " " << tn << endl; 
    }
  }
}

void  CmbCalc::dumpCl(const vector<double>& InitialPower, const CL &cl, const ControlPanel &control,
                      const string fname, string tname, bool dumpbs) const {
  bool tensors = control.tensor;
  if (tname.empty()) tensors = false;
  ofstream scalarFile, tensorFile;
  if (control.scalar) scalarFile.open(fname.c_str());
  if (tensors) tensorFile.open(tname.c_str());

  for (unsigned int n = 0; n < InitialPower.size(); ++n) {
    bool wasArmed = cl.ts[n]->isArmed();
    if (!wasArmed)
      cl.ts[n]->arm(Spline::thoseReady);
    for (int l = 2; l <= cl.ts[n]->stop() -300; ++l) {
      if (control.scalar) {
        scalarFile.width(4);
        scalarFile << l << " ";
        scalarFile.precision(7);
        scalarFile.width(13); scalarFile << (*cl.ts[n])(l) << " ";
        scalarFile.width(13); scalarFile << (*cl.es[n])(l) << " ";
        if(dumpbs) {scalarFile.width(13); scalarFile << (*cl.bs[n])(l) << " ";}
        scalarFile.width(13); scalarFile << (*cl.cs[n])(l);
        scalarFile << endl;
      }
      if (tensors) {
	tensorFile << l << " " << (*cl.tt[n])(l);
	tensorFile << " " << (*cl.et[n])(l) << " " << (*cl.bt[n])(l) << " " << (*cl.ct[n])(l) << endl; 
      }
    }
    if (!wasArmed)
      cl.ts[n]->disarm(Spline::thoseReady);
    if (control.scalar) scalarFile << endl << endl;
    if (tensors) tensorFile << endl << endl;
  }
  if (control.scalar) scalarFile.close();
  if (tensors)  tensorFile.close();
}


/*!
  Reads a file of bessel-functions. The file should previously be computed
  using the program jlgen. 
  In contrast to CMBFAST, both the multipoles and the bessel functions 
  are entirely encoded in this file. A multipole of -1e10 means end of file.
  Other than that, you are free to use any number an value of multipoles.
  For example, one may use jlgen to calculate multipoles up to l=4000 and
  save this to some file (usually called "jlgen.dat"). 
  When running CMBEASY, one can the specify the MaxMultipole to
  which the spectrum is calculated. For instance MaxMultipole=1500 will
  output the spectrum up to l=1500.
  
  The MaxMultipole feature didn't work as promised in versions < v2.0, sorry for that.

  The ability of Splines to initialize from disk is used.

  ** Attention: the maximum argument that will be used is set to be no more than
  2 * MaxMultipole. As the argument is x = k*tau, this means that the maximum wave
  number that will later be calculated in oneK() etc, is at most 2 * MaxMultipole / tau_0.

  This is necessary, because you can create a very large bessel function table on
  disk and use only some low multipoles (say up to l = 1000). If three is no file on
  disk, initjl will create one for l = 5000 and x up tp 10000. 
  Yet, if you only want l = 1000, it will only use the  first bessel functions up to 1000,
  but in order to prevent the calculation of very high k (which have no contribution),
  we limit maxBesselX to 2*Multipole.

*/

void CmbCalc::initjl(string filename, int MaxMultipole,bool recursive) {
  //static float x,y;
  BesselAnchor.kill();
  ifstream data(filename.c_str());
  if (!data) { 
    if (recursive)  throw Bad_Error("CmbCalc::initijl():  Bessel-File not found. Presumed filename was: " + filename + "\nI tried to pre-calculate a file, but I failed doing so.");
    cout << "** initjl did not find the pre-calculated bessel functions ** " << endl;
    cout << "** it will therefore create a file with multipoles up to 5000 ** " << endl;
    precalculateBesselFunctions(filename,5000,10000); // large and mighty
    initjl(filename,MaxMultipole,true);  // call this routine once more
    return;   // for this run, return
  } 
  maxBesselX = 0;  // init this 
  int l,count=-1;
  while (data.read((char*)&l,sizeof(int))) {
    if (l > MaxMultipole + 300) break;
    j2l[++count] = l;  // better for bookkeeping. In principle, duplicated information
    bessel[l] = new Spline(data,"besselfunction",&BesselAnchor);
      bessel[l]->arm();
      maxBesselX = max(maxBesselX, bessel[l]->stop());
  }
  // we limit here the maximum x = k*tau_0 
  // and hence the k_max that is caculated later on
  // to MaxMultipole * 2
  maxBesselX = min(maxBesselX, MaxMultipole*2.0);
  //  cout << "max bessel x: " << maxBesselX << endl;
  data.close();  
  if (bessel.rbegin()->first < MaxMultipole +300) 
    {
      cerr << bessel.rbegin()->first  << " :: " << MaxMultipole +300 << endl;
    throw Bad_Error("CmbCalc::initjl()  The Besselfunction file does not have bessel functions with high enaugh multipoles\nUse jlgen to precalculate these and\ngive the name of the data file generated by jlgen to me"); 
    }
} 

vector<int> CmbCalc::lval() {
  vector<int> v(bessel.size());
  int j=-1;
  for (map<int,Spline*>::iterator i=bessel.begin(); i != bessel.end(); i++) 
    v[++j] = i->first;
  return v;
} 

/*!
  The routine formerly known as initlval :-)
  
  Given a MaxMultipole, return a  list<int> that
  contains all multipole numbers for which bessel
  functions should be pre-calculated.

  If you want denser sampling of l's, this is the
  routine!

*/
list<int> CmbCalc::generateMultipoleValues(int MaxMultipole) {
  list<int> multipoles;
  for (int l=2; l <= 10; l++) multipoles.push_back(l);
  multipoles.push_back(12);
  multipoles.push_back(15);
  for (int l = 20; l <= 90; l+=10) multipoles.push_back(l);
  multipoles.push_back(110);
  multipoles.push_back(130);
  for (int l = 150; l <= MaxMultipole + 300; l += 50) multipoles.push_back(l);
  return multipoles;
}

/*!
  Calculate the spherical bessel functions and store on disc. Actually,
  with today's computing power reading from disc and calculating on the fly
  is not so terribly different, at least for moderate l.
  
  Please note that the maximum l does not any longer determine the
  number of multipoles used in the Cl calculation. Nor does keta determine
  the maximum wave number for the actual calculation. The reason for this
  is the fact that you can specify the maximum multipole to use for the
  calculation in initjl().

  The values given here as arguments are  however the limiting maximum values,
  i.e. initjl() cannot access a higher multipole than one pre-calculated here.
 
  \param filename Usually "jlgen.dat"
  \param MaxMultipole The highest multipole for which j_l's will be tabulated
  \param keta The name is legacy, in our code it should be called ktau, cause this is the maximum k times conformal time tau.
*/
void CmbCalc::precalculateBesselFunctions(string filename,int MaxMultipole, int keta) {
  ofstream data(filename.c_str());
  list<int> multipoles = generateMultipoleValues(MaxMultipole);
  int kmax = keta + 126;
  
  for (list<int>::iterator i = multipoles.begin(); i != multipoles.end();i++) {
    int l = *i;
    Model::write<int>(data, l); // write multipole to file
    Spline *TmpBessel = new Spline(20000,"Temporary bessel spline"); 
    float x,xlim,ajl;
    for (int k = 1; k <= kmax; ++k) {
      if (k <= 151) {
	if (k <= 51) x = (k - 1) / 10.0;  else x = (k - 51) / 5.0 + 5.0;
      } else   x = k - 151 + 25.0;
      
      xlim = l / 20.0;
      xlim = max(xlim,(float)10.0);
      xlim = l - xlim;
      
      if (x > xlim) {
	// here we go and get it ...
	ajl = Bessel::bjl(l, x);	    
	// in principle, the Numerical Recipes' version
	// would also work. However, the results differ
	// and I am inclined to favour the original
	// jlgen's version.
	// ajl = Bessel::sphJ(lvalues1_1.l[j - 1],x+1e-8); //1e-8 to avoid 1/0	
	
	if (l == 3 && x < 0.2) ajl = 0;
	if (l > 3 && x < 0.6) ajl = 0;
	if ( l >= 5 && x < 1) ajl = 0;
	TmpBessel->set(x,ajl); // we store x and y to reduce program dependence
      }
    }
    TmpBessel->save(data);
    delete TmpBessel;
    
  }
}

