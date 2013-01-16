#include "mcmodel.h"

#include "gauge.h"
#include "lensing.h"
#include "analyzethis.h"
#include "quintcosmos.h"

#include "fastslowstepper.h"
#include "nestedsampler.h"

//#include "couplingtuner.h"
#include "decayingdedmcosmos.h"
#include "exponential.h"
#include "exponentialcoupling.h"
#include "ratra.h"
#include <math.h>
#include <iostream>

class DecayingDEDMModel: public McCustomModel<DecayingDEDMCosmos>
{
  void compute();
};

void McModel::initialize()
{
  // tell the McSettings class which model to compute
  // change this to match the name of your model
  cfg().setModel(new DecayingDEDMModel());

  /*  NestedSampler* master = new NestedSampler();
  cfg().setMaster(master);
  master->setLivePointsCount(400);
  master->setSamplingMethod(NestedSampler::EllipsoidalSampling);
  master->setEllipsoidalEnlargementFactor(1.15);
  master->setMaxRemainingLogEvidenceFraction(0.1); */

  //cfg().useDataSet("Astier");
  //cfg().useDataSet("Riess06");
  //cfg().useDataSet("UnionSNeData08");
  cfg().useDataSet("ConstitutionSNeData09");

  //cfg().useDataSet("AcbarData");
  //cfg().useDataSet("LyaMcDonald");

  //cfg().useDataSet("SDSS");
  //cfg().useDataSet("SDSSLRG");
  //cfg().useDataSet("SdssBAO");

  //cfg().useDataSet("WMAP5Data");
  //cfg().useDataSet("WMAP7Data");

  //  addParameter(name, lowerBound, upperBound, initialSigma)
  McTaskInfo::addMcParameter("omega_cdm", 0.1, 0.5, 0.1);
  //McTaskInfo::addMcParameter("omega_cdm1", 0.01, 0.15, 0.01);
  McTaskInfo::addMcParameter("omega_b", 0.01, 0.15, 0.01);
  McTaskInfo::addMcParameter("h", 0.50, 0.85, 0.032);
  McTaskInfo::addMcParameter("tau", 0.50, 0.85, 0.032);
  McTaskInfo::addMcParameter("lambda", 0.01, 3, 0.005);
  McTaskInfo::addMcParameter("beta", 0.01, 1, 0.005);

  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_q");
 // McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "delta_nu_max");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "omega_q");
  //McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "sigma8");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_m");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_nuNR");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_b");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_cdm");
  //McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "Omega_e");
  //McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "h");

  // add additional things you would like to keep track of
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "age in years");

  // if any of the datasets need cmb or transfer functions,
  // add these paramters and keep track of sigma8
  if (cfg().controlPanel().cmb || cfg().controlPanel().power_cdm) {
  //  McTaskInfo::addMcParameter("optdlss", 0.0, 0.3, 0.03);
  //  McTaskInfo::addMcParameter("n", 0.8, 1.4, 0.02);
  //  McTaskInfo::addMcParameter("ln (10^10 A_s) - 2\\tau", 2.5, 3.2, 0.03);
  //  McTaskInfo::parameterInfo("ln (10^10 A_s) - 2\\tau").isFast=true;
  }
  if (cfg().controlPanel().power_cdm) {
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "sigma8");
  }
}

void DecayingDEDMModel::compute()
{
  McTaskInfo& param = *mParams;

  mCosmos->reset();

  bool computeCmb = mControlPanel->cmb;
  bool computePowerCdm = mControlPanel->power_cdm;
  bool computePerturbations = computeCmb || computePowerCdm;

  mCmbCalc->setGauge(Gauge::speedyCoupledInvariant);

  double h = param("h");
  mCosmos->seth(h);
  const double h2 = h*h;
  mCosmos->setT_cmb(2.725); //Temperature of the CMB
  mCosmos->setY_he(0.24);  // Primordial He-Abundance
  mCosmos->setNuR(0.);
  mCosmos->setNuNR(3); 

  if (param.hasEntry("optdlss")) {
    mCosmos->setOptDistanceLss(param("optdlss"));
  } else {
    mCosmos->setOptDistanceLss(0.08);
  }

  mControlPanel->highPrecisionTransfer = true;  // if at all cdm, high precision ?
  mControlPanel->transferMaxK=5*mCosmos->h();  // maximal k for cdm
  mControlPanel->transferPerLog=5;  // k-values per log interval

  double Omega_e=0.06;
	  double omega_cdm=param("omega_cdm");
	  double omega_b=param("omega_b");
  double omega_nu = param("omega_nu");
  double tau = param("tau");
    double lambda=param("lambda");
    double beta=param("lambda");

  mCosmos->setOmega_vacuum(0.); 
  mCosmos->setOmega_cdm(omega_cdm);
  mCosmos->setOmega_b(omega_b);
  mCosmos->setOmegaH2_nuNR(omega_nu*h*h);
  mCosmos->setOmega_quintessence_flat();


  if (fabs(mCosmos->omega_k())>1e-3) {
    throw Bad_Error("curvature");
  }
	// Set quintessence type
	Quintessence::Type quinttype;
	quinttype = cosmos.quintessence()->type();
	
	// Exponential potential
	if(quinttype==Quintessence::exponential){
	Exponential* ex = dynamic_cast<Exponential*>(cosmos.quintessence());
    ex->setInitialQDot(0.);
    ex->setLambda(lambda);
    ex->setV0(5.e-7/cosmos.M_p(2));
    ex->setInitialQ(ex->initialQ(a));
    //cout << " set ExponentialPotential" << endl;
		}

  	// Set coupling stuff
    vector<double> v(1);
    ExponentialCoupling c;
    c.setQuintessence(mCosmos->quintessence());
    v[0] = beta;
    c.setCouplingParameters(v);
    mCosmos->setCoupling(&c);

    mCosmos->setTauDM(tau);

  double n=0.953;

  if (computePerturbations) {
    mCosmos->setInitialPower(/*param("n")*/n);
    mCosmos->InitialTensorPower[0] = mCosmos->InitialPower[0] - 1.0;
  }

  if (mControlPanel->scalar) {
    mControlPanel->setInitialConditions(ControlPanel::adiabatic);
  }

  mCl->clear();
  mCl->resize(0);

  mCosmos->reset();

  if (computePerturbations) {
    mCmbCalc->cmbflat(mCosmos, "", *mControlPanel, *mCl);
  } else {
	  //cout << " History started... " << endl;
    mCosmos->history();
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

  double age = mCosmos->mpc2year()*mCosmos->tau2t(mCosmos->tau_0());
  param("omega_q") = mCosmos->omega_q(false);
  param("age in years") = age;
  param("true_omega_q") =  mCosmos->omega_q(true);
  param("true_omega_m") =  mCosmos->omega_m(true);
  param("tau") = tau;
  param("true_omega_b") =  mCosmos->omega_b(true);
  param("true_omega_cdm_1") = mCosmos->omega_cdm(true);
  param("true_omega_cdm_2") = mCosmos->omega_cdm(true);
  //param("sigma8") = mCosmos->sigma8[0];
  //param("omega_q_diff") =  mCosmos->omega_q(false)-mCosmos->omega_q(true);
  //param("omega_nuNR_diff") =  mCosmos->omega_nuNR(false)-mCosmos->omega_nuNR(true);
  //param("omega_m_diff") =  omega_m-mCosmos->omega_m(); 
  //param("omega_b_diff") =  mCosmos->omega_b()-omega_b;
}
