#include "quintcosmos.h"
#include "global.h"
#include <list>
#include <fstream>
#include <string>


#include "cleanvector.h"
#include "safevector.h"
#include "lensing.h"
#include "controlpanel.h"
#include "cl.h"
#include "cmbcalc.h"
#include "spline.h"
#include "DataManager.h"
#include "analyzethis.h"




/*
  Extended configuration file controlled driver written by Christian M. Mueller.
  You can either use this driver and its configuration file format or the original
  (and less capable) version.
*/
int main(int argc, char* argv[]) {
  QuintCosmos cosmos;
  CmbCalc cmbcalc;
  ControlPanel control;
  try {
    DataManager dataManager;    // this class is for reading the data file & holds all the needes values
    if (argc == 2) dataManager.setConfigFileName(argv[1]); else 
      dataManager.setConfigFileName("configuration.cfg"); // sets the file where the configuration is read from 
    dataManager.synchronize(); // read all the info from the file
    dataManager.printInfo();

    string base = dataManager.getOutputBase();

    cmbcalc.setGauge(dataManager.getGauge());

    static double riflag;
    cosmos.seth(dataManager.getHubbleH());

    string filename;  //! Used frequently for storing the names of various files
    string scalarFileName, tensorFileName, lensedFileName;
    string filejlens; // the filename of the lensing bessels

    map<int,string> transferFile;

    control.cmb = dataManager.getCMB();   // want cl - spectrum  ?
    control.power_cdm = dataManager.getCDMPowerSpectrum(); // want cdm power-spectrum ?
   
    if (control.power_cdm) { 
      control.highPrecisionTransfer = dataManager.getHighPrecisionTransfer();  // if at all cdm, high precision ?
      control.transferMaxK=dataManager.getTransferMaxK();  // maximal k for cdm  
      control.transferPerLog=dataManager.getTransferPerLog();  // k-values per log interval
      
      control.transferZ[0] = dataManager.getCDMRedshift();   // cdm-power output redshift
      transferFile[0] = base + dataManager.getTransferFileName(); // cdm-power output filename
    }
      
    // here, you specify the quintessence type you like
    // please note that whenever you choose Omega_quintessence=0
    // this setting becomes irrelevant
    cosmos.setQuintessence(dataManager.getQuintessenceType());
   
    cosmos.setT_cmb(dataManager.getCMBTemperature());
    cosmos.setY_he(dataManager.getHeAbundance());
    cosmos.setNuR(dataManager.getRelativisticNeutrinoSpecies()); 
    cosmos.setNuNR(dataManager.getNonRelativisticNeutrinoSpecies());
    
    control.transferMaxK *= cosmos.h();    
    control.recombination = dataManager.getRecombination();  // recombination

    
    // reionization: 0 means no reionization
    riflag = dataManager.getReionization();
    // reionization: optical depth to lss surface
    if (riflag == 1) cosmos.setOptDistanceLss(dataManager.getOptDistanceLss());
    // reionization: redshift and inonization fraction
    if (riflag == 2) {
      double z = dataManager.getReionizationZ();
      cosmos.setReionizationZ(z);
      cosmos.setReionizationFraction(dataManager.getReionizationFraction());
      cosmos.setReionizationZStop(max(z *  .07 - 1,0.0));
    }
    
    control.scalar = dataManager.getScalar();  // scalar pert. 
    control.tensor = dataManager.getTensor(); // tensor per.
     
    // some lines for the spectal indices
    if (control.scalar)  cosmos.setInitialPower(dataManager.getInitialPower()); // sets to canonical : one index with n = 1;
    // this means that InitialPower.size() will return 1 so all the for loops are only traversed once
    if (control.scalar && control.tensor) { //if (tensor_1.itflag == 1) {
      //  set tensor spectral index to scalar index - 1
      for (unsigned int in = 0; in < cosmos.InitialPower.size(); ++in) {
	cosmos.InitialTensorPower[in] = cosmos.InitialPower[in] - 1.0; 
      }	
    }
    if (control.tensor && (!control.scalar) ) {
      int n =1;
      cosmos.InitialTensorPower.clear();
      for (int in = 0; in < n; ++in) cosmos.InitialTensorPower[in]=-0.03;  
    }
    
    if (control.scalar && control.power_cdm ) {
      //cout << "Enter (0) unlensed Cls only" << endl;
      //cout << "Enter (1) lensed Cls, linear evolution" << endl;
      //cout << "Enter (2) lensed Cls, non-linear evolution" << endl; 
      int lensflag = dataManager.getLensflag();
      control.setLensing(lensflag);
    }
    
    if (control.isLensing()) control.setAllSkyLensing(dataManager.getAllSkyLensing());
	
      if (control.isLensing()   && (!control.power_cdm)) {
	cout << "You did not request the transfer function" << endl;
	cout << "calculation needed to do the lensing" << endl;
	cout << "you will have to start again" << endl;
	throw Bad_Error("no transfer functions requested, needed for lensing");
      }	  
      
      if (control.scalar) {
	scalarFileName = dataManager.getScalarFileName();  // unlensed output file name
	if (control.isLensing()) {
	  //  cout << "Enter output filename for LENSED SCALAR cl" << endl;
	  lensedFileName = dataManager.getLensedFileName(); // output file name for lensed spectra
	  filejlens = dataManager.getBesselFileName();  // file with bessel functions for lensing
	}
      }
      if (control.tensor) {
	// cout << "Enter output filename for TENSOR cl" << endl;
	tensorFileName=dataManager.getTensorFileName();
      }	
      
      
      if  (control.scalar)
	control.setInitialConditions(dataManager.getInitialConditions());  // initial conditions
    
      // for mixed initial conditions, we need
      // to find the relative contributions of each type.
    control.adiabaticContribution=dataManager.getAdiabaticContribution();
    control.isoCDMContribution=dataManager.getIsoCDMContribution();
    control.isoBaryonContribution=dataManager.getIsoBaryonContribution();
    control.isoNeutrinoContribution=dataManager.getIsoNeutrinoContribution();
    
    if (fabs(cosmos.omega_k()) <= .001) {
      if (control.cmb) {
 	string cmbeasydir = ControlPanel::cmbeasyDir();
	filename = cmbeasydir + "/resources/jlgen.dat";  // bessel function filename
	cmbcalc.initjl(filename, 2000 ); // maximum nr of l's
      }
    } else   throw Bad_Error("cmbeasy only supports flat universes");
   
    
    int n = cosmos.InitialPower.size();
    CL cl;
    cl.resize(n);
    
    cosmos.reset(); //reset is always good, epsecially, if you would like to run a loop of models
    
    cosmos.seth(dataManager.getHubbleH());   
    cosmos.setOmegaH2_b(dataManager.getOmegaH2_b());
    cosmos.setOmega_vacuum(dataManager.getOmegaVacuum()); 
    cosmos.setOmega_quintessence(dataManager.getOmegaQuintessence());
	   
    cosmos.setOmega_nuNR(dataManager.getOmegaNonRelativisticNeutrinos());
    cosmos.setOmega_cdm(1.0  - cosmos.omega_v() - cosmos.omega_b() - cosmos.omega_q()-cosmos.omega_nuNR());
	   
    
    if (dataManager.getOmegaQuintessence()!=0){
      if (dataManager.getQuintessenceType()==Quintessence::leaping) { cosmos.setQParameters(dataManager.getQuintessenceParameters()); 
      cosmos.tuneQuintessence(dataManager.getLeapingTune()); }
      else if (dataManager.getQuintessenceType()==Quintessence::ipl || dataManager.getQuintessenceType()==Quintessence::crossover ||
          dataManager.getQuintessenceType()==Quintessence::celestine ||
	       dataManager.getQuintessenceType()==Quintessence::crossoverfield) {
       	cosmos.setQParameters(dataManager.getQuintessenceParameters()); 
	//	cosmos.setQParameters(20.0,0.01,0.0,0.0,280);
	cosmos.tuneQuintessence();
      }
    }	

    //cosmos.reset();
    //cosmos.history(true);
    //cosmos.reset();
    //control.StopPert=0.01;

    string gaga("42");
    try {
      cmbcalc.cmbflat(&cosmos ,gaga, control, cl);
    } catch (Bad_Error x) {
      cout << "leandriver bad error:" << endl;
      cout << x.s << endl;
      return 0;
    }	  
    cout << endl;
    cosmos.printStatus();
    // we now set all splines to have the same number of data points as mother cl.ts
    // this is just needed for all the splines that have not been added data points,
    // for example tensor splines, if only scalar is wanted. 
    // for them to return zero when accessed (and even before, to be armed by
    // a call to clts ->arm( childrenToo) 

    AnalyzeThis ai;
    if (control.cmb) {
      cl.ts[0]->setChildrensN();
      ai.scaleCls(cl,0,pow(2.725e6,2));  // bring the cl's to muK^2 
    }
    // New in cmbeasy v4.0, you can set the initial amplitudes for A_s and A_t
    // or use the inflationary result 
    // and apply or not apply the wmap normalization
    vector<double> A_s, A_t;  // the Amplitudes. For each spectral index an Amplitude A_s and A_t
    ai.fiducialAmplitudes(cosmos,A_s,A_t); // Initialize the vectors A_s and A_t (convenience)
    // For our single spectral index, you can choose A_s & A_t (not needed if WMAP normalized later)
    if (!dataManager.getWmapNormalize()){
        A_s[0] = dataManager.getA_s();
        A_t[0] = dataManager.getA_t();
    }
    //cout << "got A_s and A_t: " << A_s[0] << " and " << A_t[0] << endl;
    // Apply the A_t  = -8 n_t A_s inflationary formula
    if (dataManager.getScalarTensorRatio()) ai.applyInflationaryTensorRatio(cosmos,A_s,A_t);

    // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
    // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
    // automatically]
    ai.rescaleSpectra(cosmos,control,cl,A_s,A_t);
    // finally WMAP 1 year quick normalize
    if (control.cmb && dataManager.getWmapNormalize())
        ai.quickWMAPNormalize(cosmos,control,cl); // WMAP normalize 
    if (control.cmb) cl.ts[0]->arm(Spline::all);      // arm all output splines

    // SN1a
    bool Riess04 = dataManager.getRiess04();
    bool HST = dataManager.getHST();
    bool Tonry = dataManager.getTonry();
    bool Astier05 = dataManager.getAstier05();
    bool AnySNe = Riess04 | HST | Tonry | Astier05;
    if (AnySNe) {
      // SN1a
      cout << endl;
      cout << "Supernovae\n";
      if (Riess04) cout << " Riess04: " << fixed <<
	ai.SNIaRiess04(cosmos) << endl;
      //Use best fit value of alpha from Knop03 for this set
      if (HST) cout << " HST: " << fixed <<
	ai.SNIaKnop03(cosmos,1.18) << endl;
      if (Tonry) cout << " Tonry: " << fixed << 
	ai.SNIaTonry03(cosmos) << endl;
      //For simplicity, use best fit alpha/beta from Astier 05
      if (Astier05) cout << " Astier05: " << fixed << 
	ai.SNIaAstier05(cosmos,1.52466,1.5690) << endl;
      cout << endl << endl;
      Spline lum(100, "Model::lum");
      for (double z = 0; z <= 2; z+= 0.05) lum.set(z, cosmos.luminosityDistance(z)); 
      lum.arm();
      cout << " old riess 06: " << ai.Sn1aRiess06(lum) << endl;
    }
      
    //output cl spectrum
    if (control.cmb) cmbcalc.dumpCl(cosmos.InitialPower,cl,control,base+scalarFileName,base+tensorFileName);
    // output tranfer (i.e. cdm-power-spectrum)

    if (control.isLensing()){
      Lensing lens(cosmos,control,cl,cmbcalc);
      lens.setBesselFileName(filejlens);
      lens.dumpLensedCls( base+lensedFileName, Lensing::Automatic );
    }
  
    if (control.power_cdm) { 
      // output power spectrum and transferfiles
      // the power spectrum is "ready to use" as P(k) over k/h 
      cout << "dumping cdm: " << base+"cdm" << endl;
      cosmos.dumpPower(0, "Pk_cdm_4",cosmos.power_cdm(),cosmos.z2tau(4));//cosmos.z2tau(dataManager.getCDMRedshift()));
//X       cosmos.dumpPower(0, "cdm",cosmos.power_cdm(),cosmos.z2tau(dataManager.getCDMRedshift()));
      //cosmos.dumpGrowthFactor();
      transferFile[0] = "trans.dat"; // cdm-power output filename
      cmbcalc.dumpTransfer(&cosmos,control,transferFile);
      cout << "sigma_8 = " << cosmos.sigma8[0] << endl;
    }

  cout << endl; 
  //cosmos.Tau2T->dump("tau2t");
  } catch (Bad_Error x) {
    cout << "\n\n******** BAD ERROR OCCURED *****\n\n";
    cout << x.s << "\n\n";
    cout << "******************************************\n";
  } catch (SafeVectorOutOfRange) {
    cout << "\n\n*** SafeVector out of range *** \n\n";
  }
  return 0;

}

int  MAIN__() {
  cout << "MAIN"<< endl;
  return 0;
}

