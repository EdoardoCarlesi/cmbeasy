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
#include "analyzethis.h"
#include <string>

#include "vdecosmos.h"
#include "vectorde.h"

using namespace std;
int main() {
  VdeCosmos cosmos;
  CmbCalc cmbcalc;
  ControlPanel control;

  // set your gauge here
  cmbcalc.setGauge(Gauge::speedyInvariant);
  //cmbcalc.setGauge(Gauge::synchronous);
  
  static double riflag;
  bool accept=true;
    double Amp = 1.; //2.75e-6; 

  string filename;  //! Used frequently for storing the names of various files
  string scalarFileName, tensorFileName, lensedFileName, pkFileName;
  string filejlens; // the filename of the lensing bessels
   
	scalarFileName = "cmb_spectrum_scalar.dat";  // unlensed output file name
	tensorFileName="plottensor.dat";
	lensedFileName = "lensed.dat"; // output file name for lensed spectra
	pkFileName = "matter_pk";

  map<int,string> transferFile;
  try {
    control.cmb = true;   // want cl - spectrum  ?
    control.power_cdm = true; // want cdm power-spectrum ?
  
  double h=0.62;
	double b_h2 = 0.06*h*h;
	double cdm = 0.28; 
    //cosmos.z_pk=0.103;
    //cosmos.Amplitude=1.;
    cosmos.setT_cmb(2.726);
    cosmos.setY_he(0.24);
    cosmos.setNuR(0.0);   // relativistic neutrinos species
    cosmos.setNuNR(0.0); // massive neutrinos species
    
    cosmos.initVde();
    cosmos.reset(); //reset is always good, especially, if you would like to run a loop of models
    cosmos.seth(h);   // hubble parameter
    cosmos.setOmegaH2_b(b_h2);  // omega_b * h^2
    cosmos.setOmega_cdm(cdm);
    cosmos.setOmega_nuNR(0.00); // Omega of nu non-relativistic
    cosmos.setOmega_vde_flat();  // quintessence 
    cosmos.setOptDistanceLss(0.0);  // few lines below, you can see the case of reionization redshift etc.
	//cout << "driver. set" << endl; 

    /* VDE initialization and parameters */
        string wfile="/home/edoardo/devel/cmbeasy/resources/wtable.dat";
	double w0=1./3.;
	int rep=720;
    cosmos.vde->initialize_grid(wfile,rep,w0);
//	cout << "driver. initialized " << endl; 
    //cout << "main, OM0: " << cosmos.vde->get_OmegaM0() << endl; 
   cout << "VDE initialized. " << endl; 
	control.recombination = Recombination::Recfast ;  // recombination (switches between recfast + recfastalpha)
	
    if (control.power_cdm) { 
      control.highPrecisionTransfer = true;  // if at all cdm, high precision ?
      control.transferMaxK=5*cosmos.h();  // maximal k for cdm  
      control.transferPerLog=5;  // k-values per log interval
      control.transferZ[0] = 0;   // cdm-power output redshift
      transferFile[0] = "trans.dat"; // cdm-power output filename
    }

    // reionization: 0 means no reionization
    riflag = 0;
    // reionization: optical depth to lss surface
    if (riflag == 1) cosmos.setOptDistanceLss(0.1); //take your pick ;-)
    // reionization: redshift and inonization fraction
    if (riflag == 2) {
      double z = 6.2;
      cosmos.setReionizationZ(6.2);
      cosmos.setReionizationFraction(0.1);
      cosmos.setReionizationZStop(max(z *  .07 - 1,0.0));
    }
 
    if (control.cmb) {	
      control.scalar = true;  // scalar spectrum  ? 
      control.tensor = false; // tensor spectrum ?
    
      // If you want full control, manipulate InitialTensorPower[] directly.
      // check out e.g. the the arguments (0.8,1.2,5) for setInitialPower(), these will give 5 indices, 0.8,0.9,1.0,1.1,1.2 	
      if (control.scalar) cosmos.setInitialPower(0.95); // one index with n_s = 0.95
      //  set tensor spectral index to scalar index - 1
       
      if (control.scalar && control.tensor) { 
	// use inflationary result n_t = n_s - 1
	for (uint in = 0; in < cosmos.InitialPower.size(); ++in) cosmos.InitialTensorPower[in] = cosmos.InitialPower[in] - 1.0; 
      }
      // if no scalars are wanted, but tensors are, 
      // you may explicitly specify n_t here
      if (control.tensor && (!control.scalar) ) {
	int n =1;
	cosmos.InitialTensorPower.clear();
	for (int in = 0; in < n; ++in)   cosmos.InitialTensorPower[in]=-0.03;
      }
          
      if (control.scalar && control.power_cdm ) {
	cout << "Enter (0) unlensed Cls only" << endl;
	cout << "Enter (1) lensed Cls, linear evolution" << endl;
	cout << "Enter (2) lensed Cls, non-linear evolution" << endl;
	int lensflag = 0;
control.setLensing(lensflag);
      }

      if (control.isLensing()   && (!control.power_cdm)) {
	cout << "You did not request the transfer function" << endl;
	cout << "calculation needed to do the lensing" << endl;
	cout << "you will have to start again" << endl;
	throw Bad_Error("no transfer functions requested, needed for lensing");
      }	  
      
      if (control.scalar) {
	// cout << "Enter output filename for SCALAR cl" << endl;
	if (control.isLensing()) {
	  cout << "If lensing was requested this will" << endl;
	  cout << "be the unlensed power spectrum" << endl;
	}
	cout << "scalar spectrum file name: " << scalarFileName << endl;

	if (control.isLensing()) {
	  cout << "Enter output filename for LENSED SCALAR cl" << endl;
	  filejlens = ControlPanel::cmbeasyDir("/resources/jlens.dat");  // file with bessel functions for lensing
	}
      }
      if (control.tensor) {
	cout << "Enter output filename for TENSOR cl" << endl;
	cout << "tensor spectrum file name: " << tensorFileName << endl;
      }	
    } else {
      control.scalar=true; control.tensor=false;
      control.setLensing(0);
    }
      
    if  (control.scalar)
      control.setInitialConditions(ControlPanel::adiabatic);  // initial conditions
    
    if (fabs(cosmos.omega_k()) >= .001) throw Bad_Error("cmbeasy only supports flat universes");
   
    if (control.cmb) cmbcalc.initjl(ControlPanel::cmbeasyDir("/resources/jlgen.dat"),1500);  
    
    int n = cosmos.InitialPower.size();
    CL cl;
    cl.resize(n);
	
    string gaga("42");
    try {
	//cosmos.history();
//cosmos.reset(); 
  	cmbcalc.cmbflat(&cosmos ,gaga, control, cl);
   
	cosmos.integrateComovingVolume(0.025,0.25);
 } catch (Bad_Error x) {
      cout << endl << "###########################################" << endl;
      cout << "leandriver bad error:" << endl;
      cout << x.s << endl;
      cout << "###########################################" << endl;
      accept = false;
    }	  
    //
    // we now set all splines to have the same number of data points as mother cl.ts
    // this is just needed for all the splines that have not been added data points,
    // for example tensor splines, if only scalar is wanted. 
    // for them to return zero when accessed (and even befor, to be armed by
    // a call to clts ->arm( childrenToo) 
    
    if (control.cmb && accept) {
      cl.tt[0]->dump("pretensor");
      cl.ts[0]->setChildrensN();
       AnalyzeThis ai;
      vector<double> A_s, A_t;  // the Amplitudes. For each spectral index an Amplitude A_s and A_t
      ai.scaleCls(cl,0,pow(2.725e6,2));  // bring the cl's to muK^2 
      // New in cmbeasy v4.0, you can set the initial amplitudes for A_s and A_t
      // or use the inflationary result 
      // and apply or not apply the wmap normalization
      ai.fiducialAmplitudes(cosmos,A_s,A_t); // Initialize the vectors A_s and A_t (convenience)
      // Apply the A_t  = -8 n_t A_s inflationary formula (if you want this, if not, comment out)
      ai.applyInflationaryTensorRatio(cosmos,A_s,A_t); 
      // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
     A_s[0] = 20e-10; A_t[0] = 20e-10; 
      // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
      // automatically]
      ai.rescaleSpectra(cosmos,control,cl,A_s,A_t);
      // finally WMAP normalize (if you like to, if not, comment out);
      ai.quickWMAPNormalize(cosmos,control,cl); // WMAP normalize 
            cl.ts[0]->arm(Spline::all);      // arm all output splines
         //output cl spectrum
      if (control.cmb) cmbcalc.dumpCl(cosmos.InitialPower,cl,control,scalarFileName,tensorFileName); 
      // output tranfer (and cdm-power-spectrum)
      
      // just for fun, you may want to know SNe Ia likelihoods
      // add more of the analyzethis stuff here, if you want, e.g. WMAP likelihood etc.
      //double Knop03Chisq = ai.SNIaKnop03(cosmos,AnalyzeThis::K03Lowe);
      //cout <<  "Knop03 (Low Extinction subsample): " << Knop03Chisq <<endl; 
    } 
    if (control.power_cdm && accept) {
      cout << "Printing Power Spectrum" << endl;
      //cmbcalc.dumpTransfer(&cosmos,control,transferFile);
      //  cosmos.dumpPower(0,pkFileName);  // Spline dumps it, hence filename is cdm.dat
	//cosmos.k2Pk->makeProper();
	//cosmos.k2Pk = cosmos.correctPower(cosmos.Amplitude, 0, "power_cdm", cosmos.k2Pk, 0, cosmos.z2tau(cosmos.z_pk));
	//cosmos.k2Pk->arm();
//for(int i=0; i<100; i++) {cosmos.k2Pk->mul(i, 1./Amp);} // cout << cosmos.k2Pk->fastY(i) << " amp: " << 1/Amp <<endl;}
	//cosmos.k2Pk->dump();
    }
    if (control.isLensing()){
        Lensing lens(cosmos,control,cl,cmbcalc);
        lens.setBesselFileName(filejlens);
        lens.dumpLensedCls(lensedFileName);
    }

  } catch (Bad_Error x) {
    cout << "\n\n******** BAD ERROR OCCURED *****\n\n";
    cout << x.s << "\n\n";
    cout << "******************************************\n";
  } catch (SafeVectorOutOfRange) {
    cout << "\n\n*** SafeVector out of range *** \n\n";
  }
    cosmos.printStatus();
    cosmos.dumpSplines();
    cout << "done" << endl;
return 0;
}
