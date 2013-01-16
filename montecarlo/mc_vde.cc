#include "mcmodel.h"

#include "gauge.h"
#include "lensing.h"
#include "analyzethis.h"
#include "../vde/vdecosmos.h"
#include "../vde/grid.h"

#include "fastslowstepper.h"
#include "nestedsampler.h"

#include <iostream>
#include <string>

class VdeCosmosModel: public McCustomModel<VdeCosmos>
{ 
  private: 
  GridXYZ *g;

  public:
  VdeCosmosModel(GridXYZ*);
  void compute();
};

VdeCosmosModel::VdeCosmosModel(GridXYZ *gr) {
g=gr;
}

void McModel::initialize()
{ 
  cout << "initialize()" << endl;

  string fileW="/home/edoardo/devel/cmbeasy/resources/wtable.dat";
  double w0 = 1./3.;
  int rep = 720;
  GridXYZ *g; 
  g = new GridXYZ;
  g->set_file_name(fileW);
  g->set_repeat(rep);
  g->set_w0(w0);
  g->initialize_grid();

  // tell the McSettings class which model to compute
  // change this to match the name of your model
  cfg().setModel(new VdeCosmosModel(g));
  //cfg().useDataSet("Astier");
  //cfg().useDataSet("Riess06");
  //cfg().useDataSet("Constitution09");
  //cfg().useDataSet("Union2");
  //cfg().useDataSet("LegacySN");
  //cfg().useDataSet("AcbarData");
  //cfg().useDataSet("LyaMcDonald");
  cfg().useDataSet("SDSS");
  //cfg().useDataSet("SDSSLRG");
  //cfg().useDataSet("SdssBAO");
  //cfg().useDataSet("AngularBAO");
  cfg().useDataSet("WMAP7Data");

  // addParameter(name, lowerBound, upperBound, initialSigma)
  //McTaskInfo::addMcParameter("omega_mh2", 0.05, 0.24, 0.01);
  //McTaskInfo::addMcParameter("omega_bh2", 0.01, 0.05, 0.002);
  McTaskInfo::addMcParameter("omega_c", 0.2, 0.6, 0.05);
  McTaskInfo::addMcParameter("omega_b", 0.02, 0.1, 0.0075);
  McTaskInfo::addMcParameter("h", 0.5, 0.8, 0.05);

  // add additional things you would like to keep track of
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "age_in_years");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "eos_today");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "initial_A0");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "omega_vde");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "omega_m");

  // if any of the datasets need cmb or transfer functions,
  // add these paramters and keep track of sigma8
  if (cfg().controlPanel().cmb || cfg().controlPanel().power_cdm) {
    McTaskInfo::addMcParameter("optdlss", 0.0, 0.3, 0.03);
    McTaskInfo::addMcParameter("n", 0.8, 1.4, 0.02);
    McTaskInfo::addMcParameter("ln (10^10 A_s) - 2\\tau", 2.5, 3.2, 0.03);
    McTaskInfo::parameterInfo("ln (10^10 A_s) - 2\\tau").isFast=true;
    //McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "sigma8");
  }
  if (cfg().controlPanel().power_cdm) {
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "sigma8");
  }
}

void VdeCosmosModel::compute()
{
  McTaskInfo& param = *mParams;
 // mControlPanel->power_cdm=true;

  mCosmos->reset();
  cout << "VdeCosmosModel::compute()" << endl;
  bool computeCmb = mControlPanel->cmb;
  bool computePowerCdm = mControlPanel->power_cdm;
  bool computePerturbations = computeCmb || computePowerCdm;

  mCmbCalc->setGauge(Gauge::speedyInvariant);
  //mCmbCalc->setGauge(Gauge::synchronous);
	mCosmos->initVde();
  // always set hubble first, because setting omega_b h^2 etc will depend on h()
  mCosmos->seth(param("h"));
  // fixed parameters (can be converted to variable parameters of course)
  mCosmos->setT_cmb(2.725); //Temperature of the CMB
  mCosmos->setY_he(0.24);  // Primordial He-Abundance
  mCosmos->setNuR(3.04); // Number of (massless) relativistic neutrios
  mCosmos->setNuNR(0);  // Number of non-relativistic neutrinos
  mCosmos->setOmega_nuNR(0.00); // Contribution of non-relativistic neutrinos

  mCosmos->vde->setGrid(g);
  mCosmos->z_pk=0.103;
// Normalize perturbations to this value
  mCosmos->Amplitude=1.; 

  if (param.hasEntry("optdlss")) {
    mCosmos->setOptDistanceLss(param("optdlss"));
  } else {
    mCosmos->setOptDistanceLss(0.08);
  }

  //cout << "VdeCosmosModel::compute()" << endl;
  double omega_c = param("omega_c");
  double omega_b = param("omega_b");
  double omega_m = omega_b + omega_c;
  double h = param("h");

cout << "VdeCosmosModel::compute(). Omega_c, Omega_b, h: " << omega_c << " " << omega_b << " " << h << endl;

  mCosmos->setOmega_cdm(omega_c);
  mCosmos->setOmega_b(omega_b);
  mCosmos->setOmegaH2_b(omega_b*h*h);
  mCosmos->setOmegaH2_cdm(omega_c*h*h);
  mCosmos->setOmega_vde_flat();

  mControlPanel->highPrecisionTransfer = false;  // if at all cdm, high precision ?
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

	//cout << "VdeCosmosModel::compute_flat()" << endl;
 mCosmos->reset();

  if (computePerturbations) {
	//cout << "VdeCosmosModel::compute_perturbations(). history" << endl;
    mCmbCalc->cmbflat(mCosmos, "", *mControlPanel, *mCl);
	//cout << "VdeCosmosModel::compute(). perturbations " << endl;
  } else {
    mCosmos->history();
  }

  double w_0 = mCosmos->vde->w(1.);
  cout << "w0 today:  " << w_0 << endl;
  param("eos_today") = w_0;
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
	//cout << "VdeCosmosModel::ai_compute_cmb()" << endl;
    // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
    // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
    // automatically]

    ai.rescaleSpectra(*mCosmos, *mControlPanel, *mCl, A_s, A_t, true);

    if (computeCmb) {
      mCl->ts[0]->arm(Spline::all);      // arm all output splines
      if (mControlPanel->isLensing()) {
        Lensing lens(*mCosmos, *mControlPanel, *mCl, *mCmbCalc);
        mCl = lens.lensedCls(); // get the lensed cl's
        mCl->ts[0]->arm(Spline::all); // arm the lensed ones
      }
    }

    if (computePowerCdm) {
//cout << "sigma8: " << mCosmos->sigma8[0] << endl;
//mCosmos->k2Pk->makeProper();
param.setEntryValue("sigma8", mCosmos->sigma8_z0);
//cout << "sigma8: " << mCosmos->sigma8_z0 << endl;
    }
  }

  // add more things here if you want to track the values of other quantities,
  // don't forget to add the corresponding entry in McModel::initializeParameters()
  double age = mCosmos->mpc2year()*mCosmos->tau2t(mCosmos->tau_0());
  param("age_in_years") = age*1.e-9;
  
  double A_0 = mCosmos->vde->getInitialA0();
  param("initial_A0") = A_0;
 
  param("omega_m") = omega_m;
  
  double omega_vde = mCosmos->vde->get_OmegaA0();
  param("omega_vde") = omega_vde;

//delete mCosmos;
	//cout << "Done task : " << ThisTask << endl;
}
