#include "cosmos.h"
#include "cncosmos.h"
#include "cmbcalc.h"
#include "exponential.h"
#include "splinetools.h"
#include <iostream>
#include <sstream>
#include "stdio.h"
#include "cninvariant.h"
#include "cleanvector.h"
#include "safevector.h"
#include "lensing.h"
#include "cl.h"
#include "spline.h"
#include "DataManager.h"
#include "analyzethis.h"

string convert(double d) {
 std::ostringstream oss;
 oss << d;
 std::string value = oss.str();
 return value;
}


int main ()
{ 
    double Q;
    CnCosmos cosmos;

    // declare CMB perturbations variables
    string base = ControlPanel::cmbeasyDir("/output/");  // base for output
    DataManager dataManager;
    string filename;  //! Used frequently for storing the names of various files
    string scalarFileName, tensorFileName, lensedFileName;
    string filejlens; // the filename of the lensing bessels

    CmbCalc cmbcalc;
    ControlPanel control;

    // setting all cosmos values to default
    cmbcalc.setGauge(Gauge::cnInvariant);
    cosmos.setQuintessence(Quintessence::exponential);
    cosmos.setT_cmb(2.725);
    cosmos.setY_he(0.24);
    cosmos.setNuR(0.0); 
    cosmos.setNuNR(3.04);
    
    // check quantities to be printed
    control.transferMaxK *= cosmos.h();   
    control.power_cdm=false; 
    control.cmb=true;
    
    cosmos.setOptDistanceLss(0.088); // default: 0.088
    control.scalar = true;  // scalar pert. 
    control.tensor = false; // tensor per.     
    cosmos.setInitialPower(0.951); // default: 0.951 //  sets to canonical : one index with n = 1;

    // cosmological parameters (SNIa)
    double omega_cdm=22.;
    double omega_b=4.6;
    double omega_m=omega_cdm + omega_b;
    double omega_nuNR=0.1;
    double omega_q=100. - omega_b - omega_cdm - omega_nuNR;
    double hub=71.;
    double a=1.e-12;

    double lambda=205.16;
    double beta=50.0;
    double index=-1;
    double intUnits = pow(2,-0.5)*pow(10,4.5); // factor to convert to internal units used by cninvariant class
    					       // which is the inverse sqrt of A_s = 20*10^-100
    
   double mDelta=0.15*intUnits; 
   double mEta=0.01*intUnits; 
    
    double norm=20e-10;
    double omega_v=0.;

    CnInvariant::maxEta=mEta;
    CnInvariant::maxDelta=mDelta;
  
    omega_cdm /= 100;
    omega_b /= 100;
    hub     /= 100;
    omega_q /= 100;
    omega_nuNR /= 100;
    omega_v /=100;

    cosmos.seth(hub);
    cosmos.setOmega_cdm(omega_cdm);
    cosmos.setOmega_b(omega_b);
    cosmos.setOmega_nuNR(omega_nuNR);
    cosmos.setOmega_quintessence_flat();
    	
    	cout << "This gives omega_q: " << cosmos.omega_q(false) << endl;
 	cout << "Omega_cdm(): " << cosmos.omega_cdm() << endl;
	cout << "Omega_nu(): " << cosmos.omega_nu() << endl;
	cout << "Omega_nuNR(): " << cosmos.omega_nuNR()<<endl;
     	
	Exponential* exp = dynamic_cast<Exponential*>(cosmos.quintessence());
    if (!exp) {
      cout << "main():  expected exponential quintessence." << endl;
      return 0;
    }

    exp->setInitialQDot(0.);
    exp->setLambda(lambda);
    exp->setV0(5.e-7/cosmos.M_p(2));

    cosmos.setNeutrinoCoupling(beta);
    Q=exp->initialQ(a);
    exp->setInitialQ(Q);
    cout << "set initial Q: " << Q << endl;
    cosmos.reset();
    cosmos.tuneQuintessence();
    
	cosmos.history();

    if (control.cmb) cmbcalc.initjl(ControlPanel::cmbeasyDir("/resources/jlgen.dat"),1500);  // max nr of l's 
    if (control.power_cdm) {
      control.highPrecisionTransfer = true;  // if at all cdm, high precision ?
      control.transferMaxK=5*cosmos.h();  // maximal k for cdm  
      control.transferPerLog=5;  // k-values per log interval
    }

    if (fabs(cosmos.omega_k())  > .001) throw Bad_Error("cmbeasy only supports flat universes");
  
    int n= cosmos.InitialPower.size();
    CL cl;
    cl.resize(n);
    
       //cosmos.reset();  //reset is always good, especially, if you would like to run a loop of models
    
       bool accept=true;     

/*	  Spline lum(100, "luminosity");
for (double z = 0.05; z <= 2; z+= 0.05) lum.set(z, 7 +5*log10(cosmos.luminosityDistance(z)));// *cosmos.H_0_cpm())); 
	  lum.arm();
	  lum.dump();
*/
	cosmos.reset();
    // start computing the whole evolution
    try {
        cmbcalc.cmbflat(&cosmos ,"42", control, cl);
    } catch (Bad_Error x) {
      cout << "useFile()  bad error:" << endl;
      cout << x.s << endl;
      cout << "It happened in cmbflat()" << endl;
      accept = false;
    }	 

    AnalyzeThis ai;
    if (control.cmb) {
      cl.ts[0]->setChildrensN();
      ai.scaleCls(cl,0,pow(2.725e6,2));  // bring the cl's to muK^2 
    }
    // New in cmbeasy v4.0, you can set the initial amplitudes for A_s and A_t
    // or use the inflationary result 
    // and apply or not apply the wmap normalization
/*	cosmos.reset();

	ofstream lumDist;
	string lum = "luminosity.dat";
	lumDist.open(lum.c_str(), std::ios_base::app);
	for(int zi; zi<20; zi++ ){
	double z = zi/10;
	double dist = cosmos.luminosityDistance(z);
	lumDist << z << " " << dist << endl;
	}
	lumDist.close();
*/
    vector<double> A_s, A_t;  // the Amplitudes. For each smixedpectral index an Amplitude A_s and A_t
    ai.fiducialAmplitudes(cosmos,A_s,A_t); // Initialize the vectors A_s and A_t (convenience)
    // For our single spectral index, you can choose A_s & A_t (not needed if WMAP normalized later)

       A_s[0] = norm;//dataManager.getA_s();
        A_t[0] = norm;//dataManager.getA_t();

    cout << "got A_s and A_t: " << A_s[0] << " and " << A_t[0] << endl;
    // Apply the A_t  = -8 n_t A_s inflationary formula
    //    if (dataManager.getScalarTensorRatio()) 
    ai.applyInflationaryTensorRatio(cosmos,A_s,A_t);

    // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
    // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
    // automatically]
    ai.rescaleSpectra(cosmos,control,cl,A_s,A_t);

    if (control.cmb) cl.ts[0]->arm(Spline::all);      // arm all output splines

 string scal= "cmb_spectra.dat" ;
 string power_sp = "power_sp_nu_scalar_.dat";
	 
    if (control.cmb) cmbcalc.dumpCl(cosmos.InitialPower,cl,control,scal,"tensor.dat"); //base+scalarFileName,base+tensorFileName);
    if (control.power_cdm) cosmos.dumpPower(0,power_sp);
	cosmos.printStatus();    
	cout << " EXP COUPLING() : " << cosmos.neutrinoCoupling() << endl;
     dumpOmegas(cosmos);
     dumpDensities(cosmos);

  return 0;
}
