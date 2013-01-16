#include "cosmos.h"
#include "coupledquintcosmos.h"
#include "mass_function.h"
#include "cmbcalc.h"
#include "coupling.h"
#include "exponential.h"
#include "exponentialcoupling.h"
#include "splinetools.h"
#include <iostream>
#include <sstream>
#include "stdio.h"
#include "ratra.h"
#include "speedycoupledinvariant.h"
#include "analyzethis.h"
#include "cl.h"
#include "spline.h"
#include "DataManager.h"

//#define LCDM

int main ()
{ 
	int stop=1;
#ifdef LCDM
	cout << " * * * LCDM * * * " << endl;
    Cosmos cosmos;
    //CoupledQuintCosmos cosmos;
#else
    CoupledQuintCosmos cosmos;
    //QuintCosmos cosmos;
#endif

    string base = ControlPanel::cmbeasyDir("/output/");  // base for output
    string filename;  //! Used frequently for storing the names of various files
    string scalarFileName, tensorFileName, lensedFileName;
    string filejlens; // the filename of the lensing bessels
	
    // Perturbation variables
    CmbCalc cmbcalc;
    ControlPanel control;

    // Setting all cosmos values to default
#ifdef LCDM
    cmbcalc.setGauge(Gauge::synchronous);
#else
    cmbcalc.setGauge(Gauge::speedyCoupledInvariant);
    //cmbcalc.setGauge(Gauge::speedyInvariant);
#endif
    
    // Set Quintessence type
    //cosmos.setQuintessence(Quintessence::exponential);
#ifndef LCDM
    cosmos.setQuintessence(Quintessence::ipl);
#endif

    cosmos.setT_cmb(2.725);
    cosmos.setY_he(0.24);
    cosmos.setNuR(3.0); 
    cosmos.setNuNR(0.0);

    // Perturbation calculations
    control.power_cdm=false; 
    control.cmb=false;
    control.scalar = false; 

if(true)
{
    control.power_cdm=true; 
    control.cmb=true;
    control.scalar = true; 
} 
   control.tensor = false;    
    control.setInitialConditions(ControlPanel::adiabatic);
	
	double index=0.951;
	double optdlss=0.088;

    cosmos.setOptDistanceLss(optdlss);
    cosmos.setInitialPower(index); 

    // Background cosmological parameters tuning
    // RP5 cdm = 44.75; ampli = 1.35; beta=0.2
    // RP4 cdm = 36.38; ampli = 1.28; beta=0.16
    // RP3 cdm = 30.89; ampli = 1.21; beta=0.12
    // RP2 cdm = 27.85; ampli = 1.14; beta=0.08
    // RP1 cdm = 26.05; ampli = 1.07; beta=0.04
    
// Final Quint: 0.73 DM: 0.225 b: 0.046
	// uDE ::cdm=21.3   b=4.4  beta=0.00
	// cDE1::cdm=20.017 b=3.9  beta=0.05 
	// cDE2::cdm=20.7   b=3.75 beta=0.10
	// cDE3::cdm=23     b=3.53 beta=0.15
	// RP5::cdm=38.3    b=4.57 beta=0.2

// Final Quint: 0.73 DM: 0.225 b: 0.046
// 	// LCDM  ::cdm=
	// cDE000::cdm=21.9   b=4.53 beta=0.000
	// cDE033::cdm=21.    b=4.29 beta=0.033 
	// cDE066::cdm=20.34  b=3.96 beta=0.066
	// cDE099::cdm=20.65  b=3.77 beta=0.099
double kkk=0.22; // Growth factor scale

#ifndef LCDM
    double omega_cdm=20.65; 
    double omega_b=3.77;
    double omega_m=omega_cdm + omega_b;
    double omega_nuNR=0.;
    double omega_nu=0.;
    double omega_q=100. - omega_b - omega_cdm - omega_nu;
    double hub=70.;
    double a=1.e-12;
#else
    double omega_cdm=21.9; 
    double omega_b=4.53;
    double omega_m=omega_cdm + omega_b;
    double omega_nuNR=0.;
    double omega_nu=0.;
    double omega_q=100. - omega_b - omega_cdm - omega_nu;
    double hub=70.;
    double a=1.e-12;
#endif
    // Coupling parameter (both EXP and IPL)
   double beta=0.099;
    // Parameter for EXP potential
    double lambda=20.; 
    
    // Parameter for IPL potential
    double alpha=0.134; 
    double amplitude=1.00e-07;

    // Amplitude normalizations
   double norm1=1.0e-10;
   double norm=20.;

    // Vacuum and Quintessence stuff
    double initialQ;
    double omega_v=0.;
    double yes=0;
    double A=0;
  
    omega_cdm /= 100;
    omega_b /= 100;
    hub     /= 100;
    omega_q /= 100;
    omega_nuNR /= 100;
    omega_nu /= 100;
    omega_v /=100;

    cosmos.seth(hub);
    cosmos.setOmega_nu(omega_nu);
    cosmos.setOmega_cdm(omega_cdm);
    cosmos.setOmega_b(omega_b);
    cosmos.setOmega_nuNR(omega_nuNR);

#ifdef LCDM
    cosmos.setOmega_vacuum_flat();
  //  cosmos.setOmega_quintessence(0);
#else
    cosmos.setOmega_vacuum(0);
    cosmos.setOmega_quintessence_flat();
#endif

	// Return input values (check ?) 
    	//cout << "This gives omega_q: " << cosmos.omega_q(false) << endl;
    	cout << "This gives omega_v: " << cosmos.omega_v() << endl;
 	cout << "Omega_cdm(): " << cosmos.omega_cdm() << endl;
 	cout << "Omega_nu(): " << cosmos.omega_nu() << endl;
	cout << "Omega_nuNR(): " << cosmos.omega_nuNR()<<endl;

#ifndef LCDM
	// check wether we use Exponential or IPL potential
	Quintessence::Type quinttype;
	quinttype = cosmos.quintessence()->type();
	
	// Exponential potential
	if(quinttype==Quintessence::exponential){
	Exponential* ex = dynamic_cast<Exponential*>(cosmos.quintessence());
    ex->setInitialQDot(0.);
    ex->setLambda(lambda);
	 A = amplitude*cosmos.M_p(0);
    ex->setV0(A); // *cosmos.M_p(2));
    ex->setInitialQ(ex->initialQ(a));
    cout << "set ExponentialPotential" << endl;
	}
		// Ratra-Peebles IPL potential
		else if(quinttype==Quintessence::ipl){
		Ratra* rat = dynamic_cast<Ratra*>(cosmos.quintessence());
	        A = amplitude*cosmos.M_p(2);
		rat->setalpha(alpha);
		rat->setA(A);
		//rat->tuneQuintessence();
		cout << " setting InversePowerLaw potential...  " << endl;
		cout << " InitialQ: " << rat->initialQ(a) << endl;
	//	rat->printStatus();
	//	cout << " Quintessence tuned to: " << rat->A() << endl;  
	}
	// Exponential coupling stuff
    vector<double> v(1);
    ExponentialCoupling c;
    c.setQuintessence(cosmos.quintessence());
	    v[0] = beta;
	    c.setCouplingParameters(v);
	cosmos.setCoupling(&c);
    // Calculate cosmological background
#endif

	cosmos.reset();
	//cosmos.tuneQuintessence();

    if (control.cmb) cmbcalc.initjl(ControlPanel::cmbeasyDir("/resources/jlgen.dat"),1500);  // max nr of l's 

    if (control.power_cdm) {
      control.highPrecisionTransfer = true;  // if at all cdm, high precision ?
      control.transferMaxK=50*cosmos.h();  // maximal k for cdm  
      control.transferPerLog=5;  // k-values per log interval
    }

    if (fabs(cosmos.omega_k())  > .001) throw Bad_Error("cmbeasy only supports flat universes");

    int n= cosmos.InitialPower.size();
    CL cl;
    ///cl.resize(n);
    cl.resize(0);
    
       bool accept=true;     
       bool lumd=true;

 /*   // Calculate and print luminosity distance 
      // if(lumd == true) {
      Spline lum(100, "luminosity_distance");
for (double z = 0.05; z <= 2; z+= 0.05) lum.set(z, 25 + 5*log10(cosmos.luminosityDistance(z)));
	  lum.arm();
	  lum.dump();
  cosmos.reset();
}*/

cout << "************** \n" ;
    // Start computing the whole evolution
       try {
     cmbcalc.cmbflat(&cosmos ,"42", control, cl);
    } catch (Bad_Error x) {
      cout << "useFile()  bad error:" << endl;
      cout << x.s << endl;
      cout << "It happened in cmbflat()" << endl;
      accept = false;
    }

double freqYr = 1.5e+8;

double z0=0; double z1=5;
	
//cosmos.printOutputList("output_list.txt", freqYr, z0, z1);	 

cout << "Evolution calculated.\n";
cout << "************** \n" ;

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
    ai.applyInflationaryTensorRatio(cosmos,A_s,A_t);

    // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
    // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
    // automatically]
    ai.rescaleSpectra(cosmos,control,cl,A_s,A_t);
	
    if (control.cmb) {
	    cl.ts[0]->arm(Spline::all);      // arm all output splines
   cmbcalc.dumpCl(cosmos.InitialPower,cl,control,"scalar_cmb.dat","tensor_cmb.dat");
    }

	cosmos.printStatus();

#ifndef LCDM
	cosmos.dumpExtraSplines();
	cosmos.normExtraSplines();
#endif

#ifdef LCDM

	//double hub, a, tau; 
	//tau = cosmos.tau2a);
	
	double H0 = cosmos.Z2H_spline->fastY(0.);
	int size = cosmos.Z2H_spline->size();

	for(int j=0; j<size; j++) 
{
	double oldv = cosmos.Z2H_spline->ydat[j];
	cosmos.Z2H_spline->setY(j, oldv/H0);
	//cout << "spline Z2H, oldv: " << oldv << " H0: " << H0 << endl;
}
	cosmos.Z2H_spline->dump();

#endif

    /* Print stuff */
    dumpOmegas(cosmos);
    dumpDensities(cosmos);

    if(control.power_cdm){
    cosmos.dumpGrowthFactor(kkk);
	//cosmos.k2Pk->dump("k2pk"); 
   cosmos.dumpPower(0, "power_spectrum",cosmos.power_cdm(),cosmos.z2tau(0));
   cosmos.dumpMassFunction("tinker",0,true,1.); 
    }

  return 0;
}
