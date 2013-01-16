#include "cosmos.h"
#include "decayingdedmcosmos.h"
#include "cmbcalc.h"
#include "coupling.h"
#include "exponential.h"
#include "exponentialcoupling.h"
#include "splinetools.h"

#include "stdio.h"
#include "ratra.h"
#include "speedycoupledinvariant.h"
#include "analyzethis.h"

#include <iostream>
#include <sstream>

int main ()
{ 
    DecayingDEDMCosmos cosmos;

    string base = ControlPanel::cmbeasyDir("/output/");  // base for output
    string filename;  //! Used frequently for storing the names of various files
    string scalarFileName, tensorFileName, lensedFileName;
    string filejlens; // the filename of the lensing bessels
	
    // Perturbation variables
    CmbCalc cmbcalc;
    ControlPanel control;

    // Setting all cosmos values to default
    cmbcalc.setGauge(Gauge::speedyCoupledInvariant);
    
    // Set Quintessence type
    //cosmos.setQuintessence(Quintessence::ipl);
    cosmos.setQuintessence(Quintessence::exponential);

    cosmos.setT_cmb(2.725);
    cosmos.setY_he(0.24);
    cosmos.setNuR(3.0); 
    cosmos.setNuNR(0.0);

    // Perturbation calculations
    control.transferMaxK *= cosmos.h();   
    control.power_cdm=false; 
    control.cmb=false;
    control.setInitialConditions(ControlPanel::adiabatic);

    cosmos.setOptDistanceLss(0.088);
    control.scalar = false; 
    control.tensor = false;    
    cosmos.setInitialPower(0.951); 

    // Amplitude normalizations
    double norm1=1.0e-10;
    double norm=20.;

    // Background cosmological parameters (SNIa)
    double omega_cdm=38.43;
    double omega_b=3.6;
    double omega_m=omega_cdm + omega_b;
    double omega_nuNR=0.;
    double omega_nu=0.;
    double omega_q=100. - omega_b - omega_cdm - omega_nuNR;
    double hub=71.;
    double a=1.e-12;
	
    // Parameter for EXP potential
    double lambda=5.; 
    double alpha=0.143; 
    double beta=-1.1;//.1; //1.; //500.20;

    // Decaying constant
    // tau_dm=10,000 ---> 15 BNY
    cosmos.setTauDM(10000.);

    // Vacuum and Quintessence stuff
    double initialQ;
    double omega_v=0.;
    double yes=0;
    
    omega_cdm /= 100;
    omega_b /= 100;
    hub     /= 100;
    omega_q /= 100;
    omega_nuNR /= 100;
    omega_nu /= 100;
    omega_v /=100;

    cosmos.seth(hub);
    cosmos.setOmega_nu(omega_nu);
    //cosmos.setOmega_quintessence(omega_q);
    cosmos.setOmega_cdm(omega_cdm);
    cosmos.setOmega_b(omega_b);
    cosmos.setOmega_nuNR(omega_nuNR);
    cosmos.setOmega_quintessence_flat();
	 cout << "omega_q: " << omega_q << endl;

 	 // check wether we use Exponential or IPL potential
	Quintessence::Type quinttype;
	quinttype = cosmos.quintessence()->type();
	
	// Exponential potential
	if(quinttype==Quintessence::exponential){
	Exponential* ex = dynamic_cast<Exponential*>(cosmos.quintessence());
    ex->setInitialQDot(0.);
    ex->setLambda(lambda);
    ex->setV0(5.e-7*cosmos.M_p(0));
    ex->setInitialQ(ex->initialQFromRadiationAttractor(a)*pow(cosmos.M_p(),-1));
    //ex->setInitialQ(1e54);
    cout << " set ExponentialPotential" << endl;
    cout << " initial Q: " << ex->initialQFromRadiationAttractor(a)/cosmos.M_p() <<  endl;
	ex->printStatus();
	}

		// Ratra-Peebles IPL potential
		else if(quinttype==Quintessence::ipl){
		Ratra* rat = dynamic_cast<Ratra*>(cosmos.quintessence());
		double A = 1.e-7*cosmos.M_p(2);
		rat->setalpha(alpha);
		rat->setA(A);
		cout << " setting InversePowerLaw potential...  " << endl;
	}

	// Exponential coupling stuff
    vector<double> v(1);
    ExponentialCoupling c;
    c.setQuintessence(cosmos.quintessence());
	    v[0] = beta;
	    c.setCouplingParameters(v);
    cosmos.setCoupling(&c);
	
    if (control.cmb) cmbcalc.initjl(ControlPanel::cmbeasyDir("/resources/jlgen.dat"),1500);  // max nr of l's 

    if (control.power_cdm) {
      control.highPrecisionTransfer = false;  // if at all cdm, high precision ?
      control.transferMaxK=10*cosmos.h();  // maximal k for cdm  
      control.transferPerLog=5;  // k-values per log interval
    }

    if (fabs(cosmos.omega_k())  > .001) throw Bad_Error("cmbeasy only supports flat universes");

    int n= cosmos.InitialPower.size();
    CL cl;
    cl.resize(n);
    
    cosmos.reset();

    bool accept=true;     
       //bool lumd=true;
       	// Calculate and print luminosity distance 
/*
       if(lumd == true) {
       Spline lum(100, "luminosity_distance");
for (double z = 0.05; z <= 2; z+= 0.05) lum.set(z, 25 + 5*log10(cosmos.luminosityDistance(z)));
	  lum.arm();
	  lum.dump();
	cosmos.reset();
}*/
    // Start computing the whole evolution
    try {
   cosmos.history(false);
        //cmbcalc.cmbflat(&cosmos ,"42", control, cl);
	//	runOneK(&cosmos, cmbcalc, control, 1.7e-3, "k_001_");
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

    vector<double> A_s, A_t;  // the Amplitudes. For each smixedpectral index an Amplitude A_s and A_t
    ai.fiducialAmplitudes(cosmos,A_s,A_t); // Initialize the vectors A_s and A_t (convenience)
    // For our single spectral index, you can choose A_s & A_t (not needed if WMAP normalized later)

        A_s[0] = norm*norm1;
        A_t[0] = norm*norm1;

    cout << "got A_s and A_t: " << A_s[0] << " and " << A_t[0] << endl;

    // Apply the A_t  = -8 n_t A_s inflationary formula
        //if (dataManager.getScalarTensorRatio()) 
    ai.applyInflationaryTensorRatio(cosmos,A_s,A_t);

    // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
    // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
    // automatically]
    ai.rescaleSpectra(cosmos,control,cl,A_s,A_t);
	
    if (control.cmb) {
	    cl.ts[0]->arm(Spline::all);      // arm all output splines
    cmbcalc.dumpCl(cosmos.InitialPower,cl,control,"scalar_cmb.dat","tensor_cmb.dat");
    }
	// Print stuff
    cosmos.dumpAll();
	//dumpOmegas(cosmos);
    //dumpDensities(cosmos);

    if(control.power_cdm){
    //cosmos.dumpGrowthFactor(0.1);
    cosmos.dumpPower(0, "k100_RP5_power_spectrum_z60",cosmos.power_cdm(),cosmos.z2tau(60));
    }

  return 0;
}
