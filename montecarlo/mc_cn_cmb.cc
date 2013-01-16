#include "mcmodel.h"

#include "gauge.h"
#include "lensing.h"
#include "analyzethis.h"
#include "quintcosmos.h"

#include "fastslowstepper.h"
#include "nestedsampler.h"

#include "cncosmos.h"
#include "cninvariant.h"
#include "exponential.h"
#include <math.h>
#include <iostream>

class CoupledNeutrinoModel: public McCustomModel<CnCosmos>
{
  void compute();
};

void McModel::initialize()
{
  // tell the McSettings class which model to compute
  // change this to match the name of your model
  cfg().setModel(new CoupledNeutrinoModel());

  /*  NestedSampler* master = new NestedSampler();
  cfg().setMaster(master);
  master->setLivePointsCount(400);
  master->setSamplingMethod(NestedSampler::EllipsoidalSampling);
  master->setEllipsoidalEnlargementFactor(1.15);
  master->setMaxRemainingLogEvidenceFraction(0.1); */

  //cfg().useDataSet("Astier");
  //cfg().useDataSet("Riess06");
  //cfg().useDataSet("UnionSNeData08");
  //cfg().useDataSet("ConstitutionSNeData09");

  //cfg().useDataSet("AcbarData");
  //cfg().useDataSet("LyaMcDonald");

  //cfg().useDataSet("SDSS");
  //cfg().useDataSet("SDSSLRG");
  //cfg().useDataSet("SdssBAO");

  //cfg().useDataSet("WMAP5Data");
  cfg().useDataSet("WMAP7Data");

  //  addParameter(name, lowerBound, upperBound, initialSigma)
  //  McTaskInfo::addMcParameter("omega_cdm", 0.01, 0.5, 0.01);
  //  McTaskInfo::addMcParameter("omega_b", 0.001, 0.2, 0.001);
  McTaskInfo::addMcParameter("omega_nu", 0.001, 0.05, 0.001);
  //  McTaskInfo::addMcParameter("h", 0.50, 0.85, 0.032);
  //  McTaskInfo::addMcParameter("Omega_e", -0.1, 0.1, 0.03);
  McTaskInfo::addMcParameter("lambda", 10, 30, 1);
  McTaskInfo::addMcParameter("beta", 20, 50, 5);
  //McTaskInfo::addMcParameter("log_beta",1,2,0.1);
  //McTaskInfo::addMcParameter("log_delta", -4, 0, 0.5);
  McTaskInfo::addMcParameter("delta_nu_max", 1.e-4, 1.e-3, 5.e-4);
  //McTaskInfo::addMcParameter("stopIndex", -2.5, -1.5, 0.25);
  //McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "lambda");

  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_q");
  McTaskInfo::addentry(McTaskInfo::ExtraInfo, "delta_nu_max");
	  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_m");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_nuNR");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_b");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "massNu_eV_today");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "true_omega_cdm");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "Omega_e");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "h");

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

void CoupledNeutrinoModel::compute()
{
  McTaskInfo& param = *mParams;

  mCosmos->reset();

  bool computeCmb = mControlPanel->cmb;
  bool computePowerCdm = mControlPanel->power_cdm;
  bool computePerturbations = computeCmb || computePowerCdm;

  mCmbCalc->setGauge(Gauge::cnInvariant);
  CnInvariant::earliestStopZ=15;
  //CnInvariant::maxGravPotential=param("stopFrac");
  //CnInvariant::cutoffIndex=param("stopIndex");
  //double delta_nu_max = pow(10, param("log_delta"));
  double intUnits = 1.; //pow(10,4.5);
	CnInvariant::maxDelta=param("delta_nu_max"); // multiplicative factor to convert 
					       	 // to internal units used by cnInvariant
	
  double h = 0.71; //param("h");
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
  double omega_cdm = 0.26; 
//double omega_cdm=fabs(param("omega_cdm"));
  double omega_b = 0.046; 
//double omega_b=fabs(param("omega_b"));
  double omega_nu = fabs(param("omega_nu"));
  double lambda=param("lambda");

  mCosmos->setOmega_vacuum(0.); // for cosmological constant    
  mCosmos->setOmega_cdm(omega_cdm);
  mCosmos->setOmega_b(omega_b);
  mCosmos->setOmegaH2_nuNR(omega_nu*h*h);
  mCosmos->setOmega_quintessence_flat();
  if (fabs(mCosmos->omega_k())>1e-3) {
    throw Bad_Error("curvature");
  }

  mCosmos->setQuintessence(Quintessence::exponential);
  Exponential* expon = dynamic_cast<Exponential*> (mCosmos->quintessence());

  //  const double Omega_early = fabs(param("Omega_e"));
  //  const double lambda = sqrt(4.0/Omega_early);
  //  double lambda=param("lambda");

  expon->setLambda(lambda);
  //  exp->setV0(mCosmos->rho_0()*mCosmos->omega_q(false)/mCosmos->M_p(4));

  expon->setV0(5.e-7/mCosmos->M_p(2));
    double a=1.e-12;
    double Q=expon->initialQ(a);
    expon->setInitialQ(Q);
	double beta=pow(10,param("log_beta"));

  mCosmos->tuneQuintessence();
  mCosmos->setNeutrinoCoupling(beta);

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
	  cout << " History started... " << endl;
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
  //param("delta_nu_max") = delta_nu_max;
  param("age in years") = age;
  //param("lambda") =  lambda;
  //param("beta") = beta;
  param("true_omega_q") =  mCosmos->omega_q(true);
  param("true_omega_m") =  mCosmos->omega_m();
  param("true_omega_nuNR") =  mCosmos->omega_nuNR(true);
  param("true_omega_b") =  mCosmos->omega_b(true);
  param("massNu_eV_today") = mCosmos->mass_nu_eV(mCosmos->tau2phi(mCosmos->tau_0()));
  param("true_omega_cdm") = mCosmos->omega_cdm(true);
  //param("omega_q_diff") =  mCosmos->omega_q(false)-mCosmos->omega_q(true);
  //param("omega_nuNR_diff") =  mCosmos->omega_nuNR(false)-mCosmos->omega_nuNR(true);
  //param("omega_m_diff") =  omega_m-mCosmos->omega_m(); 
  //param("omega_b_diff") =  mCosmos->omega_b()-omega_b;
}
