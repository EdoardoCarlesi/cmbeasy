#include "mcmodel.h"

#include "gauge.h"
#include "lensing.h"
#include "analyzethis.h"
#include "speedydeinvariant.h"
#include "quintcosmos.h"
#include "ratra.h"
#include "fastslowstepper.h"
#include "nestedsampler.h"

//#include <iostream>

class QuintModel: public McCustomModel<QuintCosmos>
{
  void compute();
  void info();
};

void McModel::initialize()
{
  // tell the McSettings class which model to compute
  // change this to match the name of your model
  cfg().setModel(new QuintModel());

  /*
  NestedSampler* master = new NestedSampler();
  cfg().setMaster(master);
  master->setLivePointsCount(400);
  master->setSamplingMethod(NestedSampler::EllipsoidalSampling);
  master->setEllipsoidalEnlargementFactor(1.15);
  master->setMaxRemainingLogEvidenceFraction(0.1);
  */

  //cfg().useDataSet("Astier");
  //cfg().useDataSet("Riess06");
  //cfg().useDataSet("AcbarData");
  //cfg().useDataSet("LyaMcDonald");
  //cfg().useDataSet("Constitution09");
  //cfg().useDataSet("Union21");
  //cfg().useDataSet("SDSS");
  //cfg().useDataSet("SDSSLRG");
  //cfg().useDataSet("SDSSLRG_10");
  //cfg().useDataSet("SdssBAO");
  cfg().useDataSet("WMAP7Data");

  // addParameter(name, lowerBound, upperBound, initialSigma)
  McTaskInfo::addMcParameter("omega_mh2", 0.05, 0.24, 0.01);
  McTaskInfo::addMcParameter("omega_bh2", 0.016, 0.03, 0.001);
  McTaskInfo::addMcParameter("h", 0.50, 0.85, 0.032);
  McTaskInfo::addMcParameter("alpha", 0., 3., 0.01);
  McTaskInfo::addMcParameter("amplitude", 1.e-8, 1.e-6, 1.e-7);

  // add additional things you would like to keep track of
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "age in years");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "omega_q");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_q");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_m");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_b");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_cdm");
	cfg().controlPanel().cmb = true;
	cfg().controlPanel().power_cdm = false;
  // if any of the datasets need cmb or transfer functions,
  // add these paramters and keep track of sigma8
  if (cfg().controlPanel().cmb || cfg().controlPanel().power_cdm) {
    McTaskInfo::addMcParameter("optdlss", 0.0, 0.3, 0.03);
    McTaskInfo::addMcParameter("n", 0.8, 1.4, 0.02);
    McTaskInfo::addMcParameter("ln (10^10 A_s) - 2\\tau", 2.5, 3.2, 0.03);
    McTaskInfo::parameterInfo("ln (10^10 A_s) - 2\\tau").isFast=true;
  }
  if (cfg().controlPanel().power_cdm) {
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "sigma8");
  }
}

void QuintModel::compute()
{
  McTaskInfo& param = *mParams;

  mCosmos->reset();
  bool computeCmb = mControlPanel->cmb;
  bool computePowerCdm = mControlPanel->power_cdm;
  bool computePerturbations = computeCmb || computePowerCdm;

  //mCmbCalc->setGauge(Gauge::speedyInvariant);
  //mCmbCalc->setGauge(Gauge::synchronous);
  //if you are computing quintessence models, the above line should be:
  mCmbCalc->setGauge(Gauge::speedyDEInvariant);

  // always set hubble first, because setting omega_b h^2 etc will depend on h()
  mCosmos->seth(param("h"));
 
 // fixed parameters (can be converted to variable parameters of course)
  mCosmos->setT_cmb(2.725); //Temperature of the CMB
  mCosmos->setY_he(0.24);  // Primordial He-Abundance
  mCosmos->setNuR(3.04); // Number of (massless) relativistic neutrios
  mCosmos->setNuNR(0);  // Number of non-relativistic neutrinos
  mCosmos->setOmega_nuNR(0.00); // Contribution of non-relativistic neutrinos
  if (param.hasEntry("optdlss")) {
    mCosmos->setOptDistanceLss(param("optdlss"));
  } else {
    mCosmos->setOptDistanceLss(0.08);
  }

cout << "McUncoupledQuintessence::params() " << param("h") << endl;

  double omega_m_h2 = param("omega_mh2");
  double omega_b_h2 = param("omega_bh2");
  mCosmos->setOmegaH2_cdm(omega_m_h2-omega_b_h2);
  mCosmos->setOmegaH2_b(omega_b_h2);
  mCosmos->setOmega_quintessence_flat(); // for quintessence

  /* **********************************************************************************
   ** if you would like to run quintessence instead of a cosmological constant,
   ** replace Comsos by QuintCosmos, uncomment the next two lines
   ** and make sure to modify and uncomment the two lines below
   ** that set the quintessence parameters and call tuneQuintessence()
   ************************************************************************************ */
   mCosmos->setQuintessence(Quintessence::ipl);   // the model you like
  double A = mCosmos->M_p(2)*param("amplitude");
  double alpha = param("alpha");
  Ratra* rat = dynamic_cast<Ratra*> (mCosmos->quintessence());
  rat->setalpha(alpha);
  rat->setA(A);
cout << "alpha: " << alpha << " amplitude: " << A << endl;

  mControlPanel->highPrecisionTransfer = true;  // if at all cdm, high precision ?
  mControlPanel->transferMaxK=5*mCosmos->h();  // maximal k for cdm
  mControlPanel->transferPerLog=5;  // k-values per log interval

  if (computePerturbations) {
    mCosmos->setInitialPower(param("n"));
    //  set tensor spectral index to scalar index - 1
    mCosmos->InitialTensorPower[0] = mCosmos->InitialPower[0] - 1.0;
  }

  if (mControlPanel->scalar) {
    mControlPanel->setInitialConditions(ControlPanel::adiabatic);
  }

  mCl->clear();
  mCl->resize(0);

  mCosmos->reset();

  // For quintessence. You might also have to call tuneQuintessnce() depending on the model
  // mCosmos->setQParameters(task[6],task[7],task[8],pow(10.0,task[9]));
  // mControlPanel->setPhantomCrossing(true);  // is the model potentially crossing the w=-1 line ?

  if (computePerturbations) {
    mCmbCalc->cmbflat(mCosmos, "", *mControlPanel, *mCl);
  } else {
    mCosmos->history();
  }


  if (computePerturbations) {

    AnalyzeThis ai;
    if (computeCmb) {
      mCl->ts[0]->setChildrensN();
      ai.scaleCls(*mCl, 0, pow(mCosmos->T_cmb()*1e6, 2)); // analyzethis needs all  Cl's in units of muK^2
    }

    vector<double> A_s, A_t;  // the Amplitudes. For each spectral index an Amplitude A_s and A_t

    ai.fiducialAmplitudes(*mCosmos, A_s, A_t); // Initialize the vectors A_s and A_t (convenience)
    // For our single spectral index, we choose A_s
    A_s[0] = exp(param("ln (10^10 A_s) - 2\\tau") + 2.*param("optdlss"))*1e-10;
    // Apply the A_t  = -8 n_t A_s inflationary formula (if you want this, if not, comment out)
    ai.applyInflationaryTensorRatio(*mCosmos, A_s, A_t);
    // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
    // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
    // automatically]
    ai.rescaleSpectra(*mCosmos, *mControlPanel, *mCl, A_s, A_t,
        true /*don't warn that cosmos has already been normalized*/);

    if (computeCmb) {
      mCl->ts[0]->arm(Spline::all);      // arm all output splines
      if (mControlPanel->isLensing()) {
        Lensing lens(*mCosmos, *mControlPanel, *mCl, *mCmbCalc);
        mCl = lens.lensedCls(); // get the lensed cl's
        mCl->ts[0]->arm(Spline::all); // arm the lensed ones
      }
    }

    if (computePowerCdm) {
      param.setEntryValue("sigma8",  mCosmos->sigma8[0]);
      cout << " sigma_8: " << mCosmos->sigma8[0] << endl; 
   }
  }

  // add more things here if you want to track the values of other quantities,
  // don't forget to add the corresponding entry in McModel::initializeParameters()
//info();
  double age = mCosmos->mpc2year()*mCosmos->tau2t(mCosmos->tau_0());
  param("age in years") = age;
  param("omega_q") = mCosmos->omega_q(false);
  param("true_omega_q") =  mCosmos->omega_q(true);
  param("true_omega_m") =  mCosmos->omega_m(true);
  param("true_omega_b") =  mCosmos->omega_b(true);
  param("true_omega_cdm") = mCosmos->omega_cdm(true);
//  param("lambda") = mCosmos->omega_v();
}


void QuintModel::info(){
cout << " ****************** " << endl;
cout << " Omega_cdm: " << mCosmos->omega_cdm(true) << endl; 
cout << " Omega_b: " << mCosmos->omega_b(true) << endl; 
cout << " Omega_q: " << mCosmos->omega_q(true) << endl; 
cout << " ****************** " << endl;
}

