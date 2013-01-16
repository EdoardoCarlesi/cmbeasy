#include "mcmodel.h"

#include "gauge.h"
#include "lensing.h"
#include "analyzethis.h"
#include "quintcosmos.h"

#include "fastslowstepper.h"
#include "nestedsampler.h"

//#include "couplingtuner.h"
#include "coupledquintcosmos.h"
#include "speedycoupledinvariant.h"
#include "coupledinvariant.h"
#include "exponential.h"
#include "exponentialcoupling.h"
#include "ratra.h"
#include <math.h>
#include <iostream>

class CoupledDEDMModel: public McCustomModel<CoupledQuintCosmos>
//class CoupledDEDMModel: public McCustomModel<QuintCosmos>
{
  void compute();
  void info();
};

void McModel::initialize()
{
  // tell the McSettings class which model to compute
  // change this to match the name of your model
  cfg().setModel(new CoupledDEDMModel());

  /*  NestedSampler* master = new NestedSampler();
  cfg().setMaster(master);
  master->setLivePointsCount(400);
  master->setSamplingMethod(NestedSampler::EllipsoidalSampling);
  master->setEllipsoidalEnlargementFactor(1.15);
  master->setMaxRemainingLogEvidenceFraction(0.1); */

  //cfg().useDataSet("Astier");
  //cfg().useDataSet("Riess06");
  //cfg().useDataSet("UnionSNeData08");
  //cfg().useDataSet("Constitution09");
  //cfg().useDataSet("Union21");

  //cfg().useDataSet("AcbarData");
  //cfg().useDataSet("LyaMcDonald");

  //cfg().useDataSet("SDSS");
  //cfg().useDataSet("SDSSLRG");
  //cfg().useDataSet("SdssBAO");

  cfg().useDataSet("WMAP7Data");

	cfg().controlPanel().cmb       = true;
	//cfg().controlPanel().cmb       = false;
	cfg().controlPanel().power_cdm = false;
//cfg().controlPanel().power_cdm = true;

  //  addParameter(name, lowerBound, upperBound, initialSigma)
  McTaskInfo::addMcParameter("omega_cdm", 0.15, 0.33, 0.075);
  McTaskInfo::addMcParameter("omega_b", 0.03, 0.08, 0.05);
  McTaskInfo::addMcParameter("omega_nu", 0.01, 0.05, 0.01);
  McTaskInfo::addMcParameter("h", 0.50, 0.85, 0.032);
  McTaskInfo::addMcParameter("alpha", 0.1, 0.5, 0.05);
  McTaskInfo::addMcParameter("beta", 0., 0.5, 0.05);
  //McTaskInfo::addMcParameter("amplitude", 1.e-7, 1.e-5, 1.e-7);

  // add additional things you would like to keep track of
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "omega_q");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_q");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_m");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_b");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_cdm");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "age in years");
  //McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_nu");

  // if any of the datasets need cmb or transfer functions,
  // add these paramters and keep track of sigma8
  if (cfg().controlPanel().cmb || cfg().controlPanel().power_cdm) {
 cout << "Perturbations enabled\n"; 
    //McTaskInfo::addMcParameter("optdlss", 0.0, 0.1, 0.03);
    McTaskInfo::addMcParameter("n", 0.8, 1.2, 0.02);
    //McTaskInfo::addMcParameter("ln (10^10 A_s) - 2\\tau", 2.5, 3.2, 0.03);
    //McTaskInfo::parameterInfo("ln (10^10 A_s) - 2\\tau").isFast=true;
  }
  if (cfg().controlPanel().power_cdm) {
cout << "Keeping track of sigma8\n";
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "sigma8");
  }
}

void CoupledDEDMModel::compute()
{
  McTaskInfo& param = *mParams;

  bool computeCmb = mControlPanel->cmb;
  bool computePowerCdm = mControlPanel->power_cdm;
  bool computePerturbations = computeCmb || computePowerCdm;

  //mCosmos->reset();
//mCmbCalc->setGauge(Gauge::CoupledInvariant);
mCmbCalc->setGauge(Gauge::synchronous);
//mCmbCalc->setGauge(Gauge::speedyDEInvariant);
//mCmbCalc->setGauge(Gauge::speedyCoupledInvariant);

  double h = param("h");
  const double h2 = h*h;
  mCosmos->seth(h);
  mCosmos->setT_cmb(2.725); //Temperature of the CMB
  mCosmos->setY_he(0.24);  // Primordial He-Abundance
  mCosmos->setNuR(3.);
  mCosmos->setNuNR(0.); 

  if (param.hasEntry("optdlss")) {
    mCosmos->setOptDistanceLss(param("optdlss"));
  } else {
    mCosmos->setOptDistanceLss(0.08);
  }

  mControlPanel->highPrecisionTransfer = true;  // if at all cdm, high precision ?
  mControlPanel->transferMaxK=5*mCosmos->h();  // maximal k for cdm
  mControlPanel->transferPerLog=5;  // k-values per log interval

	  double omega_cdm=param("omega_cdm");
	  double omega_b=param("omega_b");
  double omega_nu = param("omega_nu");

  // Set potential parameter
    double alpha=param("alpha");
    //double amplitude=param("amplitude");
    double amplitude=1.e-7; 

  mCosmos->setOmega_vacuum(0.); // for cosmological constant    
  mCosmos->setOmega_cdm(omega_cdm);
  mCosmos->setOmega_b(omega_b);
  mCosmos->setOmega_nuNR(0.); //omega_nu);
  mCosmos->setOmega_nu(omega_nu);
  mCosmos->setOmega_quintessence_flat();

  if (fabs(mCosmos->omega_k())>1e-3) {
    throw Bad_Error("curvature");
  }
	// Set quintessence type
  mCosmos->setQuintessence(Quintessence::ipl);
  Ratra* rat = dynamic_cast<Ratra*> (mCosmos->quintessence());
  rat->setalpha(alpha);
  rat->setA(amplitude*mCosmos->M_p(2));
  double Q=rat->initialQ(1.e-12);

 	// Set coupling stuff
    vector<double> v(1);
    ExponentialCoupling c;
    c.setQuintessence(mCosmos->quintessence());
    v[0] = param("beta");
    c.setCouplingParameters(v);
    mCosmos->setCoupling(&c);

	cout << " InitialQ : " << Q         << " alpha  : " << alpha   << " beta   : " << v[0] << endl;
	cout << " Omega_cdm: " << omega_cdm << " Omega_b: " << omega_b << " Omega_q: " << mCosmos->omega_q(false) << endl;

 // double n=0.953;
 double n=param("n");

  if (computePerturbations) {
    mCosmos->setInitialPower(n);
    mCosmos->InitialTensorPower[0] = mCosmos->InitialPower[0] - 1.0;
  }

  if (mControlPanel->scalar) {
    mControlPanel->setInitialConditions(ControlPanel::adiabatic);
  }

  mCl->clear();
  //mCl->resize(n);
  mCl->resize(0);

  mCosmos->reset();

  if (computePerturbations) {
/*
cout << "mCosmos->history()\n"; 
   mCosmos->history(false);
cout << "mCosmos->history(). Done.\n"; 
   mCosmos->reset();
cout << "mCosmos->history()-2\n"; 
   mCosmos->history(false);
cout << "mCosmos->history()-2. Done.\n"; 
*/
   mCmbCalc->cmbflat(mCosmos, "42", *mControlPanel, *mCl);
  } else {
   mCosmos->history(false);
  }

  if (computePerturbations) {

    AnalyzeThis ai;
    if (computeCmb) {
      mCl->ts[0]->setChildrensN();
      ai.scaleCls(*mCl, 0, pow(mCosmos->T_cmb()*1e6, 2)); // analyzethis needs all  Cl's in units of muK^2
    }

    vector<double> A_s, A_t;  

    ai.fiducialAmplitudes(*mCosmos, A_s, A_t); 
    A_s[0] = std::exp( 2.92  +2.*0.08 /*param("ln (10^10 A_s) - 2\\tau") + 2.*param("optdlss")*/)*1e-10;
//    A_s[0] = exp(param("ln (10^10 A_s) - 2\\tau") + 2.*param("optdlss"))*1e-10;
    ai.applyInflationaryTensorRatio(*mCosmos, A_s, A_t);
    ai.rescaleSpectra(*mCosmos, *mControlPanel, *mCl, A_s, A_t,true);

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
    }
  }

  // add more things here if you want to track the values of other quantities,
  // don't forget to add the corresponding entry in McModel::initializeParameters()
//	info();
//	mCosmos->printStatus();
  double age = mCosmos->mpc2year()*mCosmos->tau2t(mCosmos->tau_0());
  param("omega_q") = mCosmos->omega_q(false);
  param("age in years") = age;
  param("true_omega_q") =  mCosmos->omega_q(true);
  param("true_omega_m") =  mCosmos->omega_m(true);
  param("true_omega_b") =  mCosmos->omega_b(true);
  param("true_omega_cdm") = mCosmos->omega_cdm(true);
  //param("true_omega_nu") =  mCosmos->omega_nu(true);
  //param("sigma8") = mCosmos->sigma8[0];
  //param("omega_q_diff") =  mCosmos->omega_q(false)-mCosmos->omega_q(true);
  //param("omega_nuNR_diff") =  mCosmos->omega_nuNR(false)-mCosmos->omega_nuNR(true);
  //param("omega_m_diff") =  omega_m-mCosmos->omega_m(); 
  //param("omega_b_diff") =  mCosmos->omega_b()-omega_b;
}

void CoupledDEDMModel::info(){
cout << " ****************** " << endl;
cout << " Omega_cdm: " << mCosmos->omega_cdm(true) << endl; 
cout << " Omega_b: " << mCosmos->omega_b(true) << endl; 
cout << " Omega_q: " << mCosmos->omega_q(true) << endl; 
cout << " ****************** " << endl;
}
