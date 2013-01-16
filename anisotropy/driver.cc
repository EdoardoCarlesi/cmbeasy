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
#include "analyzethis.h"
#include <string>

using namespace std;

/*! small sub routine.
  a map<string, double> and a key to search
  for, it will return the value corresponding to
  the key in the map. If the key is not found,
  it will return the default value given
*/
double lookUp(map<string, double>& m, string key, double deflt) {
  if (m.find(key) == m.end()) return deflt;
  return m[key];
}

/*!
  Same as above, just for string values. Too easy for template :-)
*/
string lookUp(map<string, string>& m, string key, string deflt) {
  if (m.find(key) == m.end()) return deflt;
  return m[key];
}

/*!
  Return true, if key has been specified in map<string, double>
*/
bool testKey(map<string, double>& m, string key) {
  if (m.find(key) == m.end()) return false;
  return true;
}

string killspace(string s) {
  string copy;
  for (unsigned int i = 0; i < s.length(); i++) {
    char c = s.at(i);
    if (isalnum(c)) copy += tolower(c);
    if (c == '_') copy += tolower(c);
    if (c=='.') copy += c;
  }
  return copy;
}



/*!
  Given a filename, this routine will
  read in parameter values and run cmbeasy.
  Please note that not all functionality is
  available using these instruction files.
  Feel free to extend it (and if you like, send
  me your version)
*/

void useFile(const char* file) {
  map<string, double> keys;
  map<string,string> keys2;
  cout << "USEFILE: " << file << endl;
  try {
    // First, we have to read in the  file line by line
    // and determine the key, value pairs
    ifstream in(file);
    if (!file) throw Bad_Error("useFile(): File [" + string(file) + "] not found");
    char buffer[1000];
    buffer[0]=0;
    long pos = in.tellg();
    while (in.getline(buffer,999)) {
      bool take = true;
      if (buffer[0] == '#') take = false;
      int count =0;
      while (buffer[count++] != 0);
      if (count < 5) take = false;
      if (take) {
	in.seekg(pos); // rewind
	in.getline(buffer,999,'=');
	string key(buffer);
	key = killspace(key);

	bool isstring = false;
	if (key == "scalar" ||  key == "tensor" || key == "cmb" || key == "power") isstring = true;
	if (key == "scalarfile" || key == "tensorfile" || key == "powerfile" || key == "gauge") isstring = true;
     	if (key == "transfile" || key == "initialcond" || key == "verbose")  isstring = true;
	
	if (isstring) {
	  string value;
	  in >> value;
	  keys2[key] = killspace(value);  // store string , string pair 
	} else {
	  double value = -1;
	  in >> value;
	  keys[key] = value;  // store string , double pair
	}
	in.getline(buffer,999);	
      }
      pos = in.tellg();
      buffer[0]=0;
    }
    cout << "==========================" << endl;
    for (map<string,string >::iterator i = keys2.begin(); i != keys2.end(); i++)
      cout << i->first << "    ::::   " << i->second << endl;
    for (map<string,double>::iterator i = keys.begin(); i != keys.end(); i++)
      cout << i->first << "    ::::   " << i->second << endl;
    cout << "===========================" << endl;
  
    
    
    // Here we go
    string scalarFileName, tensorFileName, powerFileName;

    Cosmos cosmos;
    CmbCalc cmbcalc;
    ControlPanel control;
    
    // the gauge
    string gauge = lookUp(keys2,"gauge","synchronous");
    bool unknown = true;
    if (gauge == "synchronous")   { cmbcalc.setGauge(Gauge::synchronous); unknown =false;}
    if (gauge == "speedyinvariant") { cmbcalc.setGauge(Gauge::speedyInvariant); unknown = false;}
    if (unknown)  throw Bad_Error("useFile(): I don't recognize the gauge in parameter file");
    
   
    // what do we want to calculate ?
    control.scalar= lookUp(keys2, "scalar", "yes") == "yes";
    control.tensor = lookUp(keys2,"tensor","no") == "yes";
    control.cmb = lookUp(keys2,"cmb","yes") == "yes";
    control.power_cdm =  lookUp(keys2, "power", "no") == "yes";
    //control.power_cdm =  true; 
    control.verbose= lookUp(keys2,"verbose","no") == "yes";
    // hubble and content
    double h = lookUp(keys,"h",0.65);
    cosmos.seth(h);
    // for baryons and cdm, you have the choice of Omega_b or Omega_b h^2:
    if (testKey(keys,"omega_b"))  cosmos.setOmega_b(lookUp(keys,"omega_b",0.05));
    else cosmos.setOmegaH2_b(lookUp(keys,"omega_bh2",0.023));
    
    if (testKey(keys,"omega_cdm")) cosmos.setOmega_cdm(lookUp(keys,"omega_cdm",0.35));
    else {
      // if matter*h^2  and not cdm has been specified, 
      // we set omega_cdm * h^2 to omega_matter*h^2 - omega_baryon * h2
      cosmos.setOmegaH2_cdm(lookUp(keys,"omega_mh2",0.15) - cosmos.omega_b()*h*h); 
    }

    cosmos.setNuR(lookUp(keys,"nur",3.04)); 
    cosmos.setY_he(lookUp(keys,"yhe",0.24));
    
    // if no vacuum energy is specified, we set omega_vac such that omega_total = 1
    if (testKey(keys, "omega_vac")) cosmos.setOmega_vacuum(lookUp(keys,"omega_vac",0.6)); 
    else cosmos.setOmega_vacuum(1.0  - cosmos.omega_cdm() - cosmos.omega_b() -cosmos.omega_nuNR());
    
    // transfer redshift
    control.transferZ[0] = lookUp(keys,"transz",0);
    if (control.power_cdm) {
      control.highPrecisionTransfer = true;  // if at all cdm, high precision ?
      control.transferMaxK=5*cosmos.h();  // maximal k for cdm  
      control.transferPerLog=5;  // k-values per log interval
    }
    
    // spectral index
    int n_number = (int)lookUp(keys, "n_number",1);
    double n_start =  lookUp(keys, "n_start",1);
    double n_stop =  lookUp(keys, "n_stop",1);
    double dnsdlnk = lookUp(keys,"n_running",0);
    

    //    control.peebles=false;  // recombination (switches between peebles and recfast)

    if (control.scalar)  {
      cosmos.setInitialPower(n_start, n_stop, n_number); 
      // running of index
      for (unsigned int i = 0; i < cosmos.InitialPower.size(); i++) 
	cosmos.InitialPower_dnsdlnk[i] = dnsdlnk;
    }	
    if (control.scalar && control.tensor) { 
	//  set tensor spectral index to scalar index - 1
      for (uint in = 0; in < cosmos.InitialPower.size(); ++in) cosmos.InitialTensorPower[in] = cosmos.InitialPower[in] - 1.0;  
    }
    if (control.tensor && (!control.scalar) ) {
      cosmos.InitialTensorPower.clear();
      cosmos.InitialTensorPower[0]=0.0;
    }

    control.setLensing(ControlPanel::none);  // no lensing
    
    // Output files
    map<int,string> transferFile;
    scalarFileName = lookUp(keys2,"scalarfile","plot.dat");
    tensorFileName = lookUp(keys2,"tensorfile","plot.dat");
    powerFileName = lookUp(keys2,"powerfile","cdm");  // as the power spectrum is dumped by a spline, this will be called cdm.dat
    transferFile[0] = lookUp(keys2,"transfile","trans.dat");

    // initial conditions
    string initialcond = lookUp(keys2,"initialcond","adiabatic");
    unknown = true;
    if (initialcond == "adiabatic")  { control.setInitialConditions(ControlPanel::adiabatic); unknown = false;}
    if (initialcond == "isocdm")   { control.setInitialConditions(ControlPanel::isoCDM); unknown =false;}
    if (unknown)  throw Bad_Error("useFile(): I don't recognize the initial conditions in parameter file"); 
      
    // set the jlgen.dat filename and initialize it
    if (control.cmb) cmbcalc.initjl(ControlPanel::cmbeasyDir("/resources/jlgen.dat"),1500);  // max nr of l's 

    if (fabs(cosmos.omega_k()) > .001) throw Bad_Error("cmbeasy only supports flat universes");

    // optical depth
    cosmos.setOptDistanceLss(lookUp(keys,"tau",0));
    
    //
    // now, do the calculation
    //
    
    int n = cosmos.InitialPower.size();
    CL cl;
    cl.resize(n);
    
    cosmos.reset();  //reset is always good, epsecially, if you would like to run a loop of models
    
    bool accept=true;
    try {
      cmbcalc.cmbflat(&cosmos ,"42", control, cl);
    } catch (Bad_Error x) {
      cout << "useFile()  bad error:" << endl;
      cout << x.s << endl;
      cout << "It happened in cmbflat()" << endl;
      accept = false;
    }	  
    cout << endl;
    cosmos.printStatus();
    
    // we now set all splines to have the same number of data points as mother cl.ts
     // this is just needed for all the splines that have not been added data points,
     // for example tensor splines, if only scalar is wanted. 
    // for them to return zero when accessed (and even before), to be armed by
    // a call to clts ->arm( childrenToo) 
    
    if (control.cmb && accept) {
      cl.ts[0]->setChildrensN();
      AnalyzeThis ai;
      ai.scaleCls(cl,0,pow(2.725e6,2));  // bring the cl's to muK^2 
      // New in cmbeasy v4.0, you can set the initial amplitudes for A_s and A_t
      // or use the inflationary result 
      // and apply or not apply the wmap normalization
      vector<double> A_s, A_t;  // the Amplitudes. For each spectral index an Amplitude A_s and A_t
      ai.fiducialAmplitudes(cosmos,A_s,A_t); // Initialize the vectors A_s and A_t (convenience)
      // For our single spectral index, you can choose A_s & A_t (not needed if WMAP normalized later)
      A_s[0] = 20e-10; A_t[0] = 20e-10; 
      // Apply the A_t  = -8 n_t A_s inflationary formula (if you want this, if not, comment out)
      ai.applyInflationaryTensorRatio(cosmos,A_s,A_t); 
      // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
      // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
      // automatically]
      ai.rescaleSpectra(cosmos,control,cl,A_s,A_t);
      // finally WMAP normalize (if you like to, if not, comment out);
      ai.quickWMAPNormalize(cosmos,control,cl); // WMAP normalize 

      //If you like to get wmap likelihoods, uncomment the next few lines 
      /* ai.initWMAPCommon(); 
       vector<WMAPNorm> like = ai.WMAPNormalize(cosmos,control,cl);
	 
       cout << "WMAP chi2's: " << endl;
       cout << like[0].chi2_tt << endl;
       cout << like[0].chi2_te << endl;
       cout << "total: " <<  like[0].chi2_tt  + like[0].chi2_te << endl;
      */
      cl.ts[0]->arm(Spline::all);      // arm all output splines
      
      //output cl spectrum
      if (control.cmb) cmbcalc.dumpCl(cosmos.InitialPower,cl, control,scalarFileName,tensorFileName); 

      //just for fun, calculate the SNe Ia likelihoods of Knop03
      double Knop03Chisq = ai.SNIaKnop03(cosmos,AnalyzeThis::K03Lowe);
      cout <<  "Knop03 (Low Extinction subsample): " << Knop03Chisq <<endl; 
    } 
    // The power spectrum will only be dumped for the first spectral index
    if (control.power_cdm && accept)  {
	powerFileName="power_cdm";
      cosmos.dumpPower(0,powerFileName);
      cmbcalc.dumpTransfer(&cosmos,control,transferFile);
      cout << "sigma_8 = " << cosmos.sigma8[0] << endl;
    }
    cout  << endl << "done" << endl;
  } catch (Bad_Error x) {
    cout << "useFile() bad error:" << endl;
    cout << x.s << endl;
  } catch (SafeVectorOutOfRange) {
    cout << "\n\n*** SafeVector out of range *** \n\n";
  }
  cout << endl;
}

// ====================== END OF useFile() =========================

/*!
  If you like full control, go call cmb without argument.
  detailed() will then be executed. To change the settings,
  however, you will need to recompile the program (personally,
  I don't find this uncomfortable). To run many programs in a loop,
  you can of course write a different driver. Or copy this routine
  and modify it a bit.
*/
void detailed() {
  // We start off by choosing a cosmos object.
  // If you do not intend to run quintessence models,
  // simply switch substitue Cosmos for  QuintCosmos.
  // Some of the functions we use in this driver, however, are
  // quintessence specific. For instance, Cosmos has no
  // setOmega_quintessence() function. You need to get rid
  // of these (there aren't many), for the program to compile then.
  // SPEED: Cosmos is considerably faster than QuintCosmos
  // mainly, because the perturbation evolution is different.
  // so it is worth switching, if you don't need it  
  QuintCosmos cosmos;
  CmbCalc cmbcalc;
  ControlPanel control;
   
  // set your gauge here
  // please note that if you have Omega_quintessence not equal to zero,
  // you need to take a quintessence perturbation class. 
  //cmbcalc.setGauge(Gauge::speedyDEInvariant);
  cmbcalc.setGauge(Gauge::synchronous);
  
  static double riflag;
  bool accept=true;

  string filename;  //! Used frequently for storing the names of various files
  string scalarFileName, tensorFileName, lensedFileName;
  string filejlens; // the filename of the lensing bessels
   
  map<int,string> transferFile;
  try {
    control.cmb = true;   // want cl - spectrum  ?
    control.power_cdm = true; // want cdm power-spectrum ?
       
    // here, you specify the quintessence type you like
    // please note that whenever you choose Omega_quintessence=0
    // this setting becomes irrelevant
    //cosmos.setQuintessence(Quintessence::ipl);
    cosmos.setQuintessence(Quintessence::none);
  
    cosmos.setT_cmb(2.726);
    cosmos.setY_he(0.24);
    cosmos.setNuR(3.04);   // relativistic neutrinos species
    cosmos.setNuNR(0); // massive neutrinos species
    
    cosmos.seth(0.7);   // hubble parameter
    //cosmos.seth(1.);  // hubble parameter
    cosmos.setOmegaH2_b(0.022);  // omega_b * h^2
    cosmos.setOmega_vacuum(0.73);   // vacuum energy
    cosmos.setOmega_quintessence(0.0);  // quintessence 
    
    cosmos.setOmega_nuNR(0.00); // Omega of nu non-relativistic
    cosmos.setOmega_cdm(1.0  - cosmos.omega_v() - cosmos.omega_b() - cosmos.omega_q()-cosmos.omega_nuNR());	   
    cosmos.setOptDistanceLss(0.0);  // few lines below, you can see the case of reionization redshift etc.
    
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
	scalarFileName = "plot.dat";  // unlensed output file name
	cout << "scalar spectrum file name: " << scalarFileName << endl;

	if (control.isLensing()) {
	  cout << "Enter output filename for LENSED SCALAR cl" << endl;
	  lensedFileName = "lensed.dat"; // output file name for lensed spectra
	  filejlens = ControlPanel::cmbeasyDir("/resources/jlens.dat");  // file with bessel functions for lensing
	}
      }
      if (control.tensor) {
	cout << "Enter output filename for TENSOR cl" << endl;
	tensorFileName="plottensor.dat";
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
    
    cosmos.reset(); //reset is always good, especially, if you would like to run a loop of models
    //  For the inverse - power - law, we chose here, there is only
    // one parameter that you may specify: the power law
    // Other models have more parameters, some may have none.
    cosmos.setQParameters(3.0); 
    cosmos.tuneQuintessence();
	
    string gaga("42");
    try {
      cmbcalc.cmbflat(&cosmos ,gaga, control, cl);
    } catch (Bad_Error x) {
      cout << endl << "###########################################" << endl;
      cout << "leandriver bad error:" << endl;
      cout << x.s << endl;
      cout << "###########################################" << endl;
      accept = false;
    }	  
    cout << endl;
    cosmos.printStatus();
    // we now set all splines to have the same number of data points as mother cl.ts
    // this is just needed for all the splines that have not been added data points,
    // for example tensor splines, if only scalar is wanted. 
    // for them to return zero when accessed (and even befor, to be armed by
    // a call to clts ->arm( childrenToo) 
    
    if (control.cmb && accept) {
      cl.tt[0]->dump("pretensor");
      cl.ts[0]->setChildrensN();
       AnalyzeThis ai;
      ai.scaleCls(cl,0,pow(2.725e6,2));  // bring the cl's to muK^2 
      // New in cmbeasy v4.0, you can set the initial amplitudes for A_s and A_t
      // or use the inflationary result 
      // and apply or not apply the wmap normalization
      vector<double> A_s, A_t;  // the Amplitudes. For each spectral index an Amplitude A_s and A_t
      ai.fiducialAmplitudes(cosmos,A_s,A_t); // Initialize the vectors A_s and A_t (convenience)
      // For our single spectral index, you can choose A_s & A_t (not needed if WMAP normalized later)
      A_s[0] = 20e-10; A_t[0] = 20e-10; 
      // Apply the A_t  = -8 n_t A_s inflationary formula (if you want this, if not, comment out)
      ai.applyInflationaryTensorRatio(cosmos,A_s,A_t); 
      // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
      // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
      // automatically]
      ai.rescaleSpectra(cosmos,control,cl,A_s,A_t);
      // finally WMAP normalize (if you like to, if not, comment out);
      ai.quickWMAPNormalize(cosmos,control,cl); // WMAP normalize 
      
      /* If you like to get wmap likelihoods, uncomment the next few lines
       ai.initWMAPCommon(); 
       vector<WMAPNorm> like = ai.WMAPNormalize(cosmos,control,cl);
	 
       cout << "WMAP chi2's: " << endl;
       cout << like[0].chi2_tt << endl;
       cout << like[0].chi2_te << endl;
       cout << "total: " <<  like[0].chi2_tt  + like[0].chi2_te << endl;
      */

      cl.ts[0]->arm(Spline::all);      // arm all output splines
      
      //output cl spectrum
      if (control.cmb) cmbcalc.dumpCl(cosmos.InitialPower,cl, control,scalarFileName,tensorFileName); 
      // output tranfer (and cdm-power-spectrum)
      
      // just for fun, you may want to know SNe Ia likelihoods
      // add more of the analyzethis stuff here, if you want, e.g. WMAP likelihood etc.
    //  double Knop03Chisq = ai.SNIaKnop03(cosmos,AnalyzeThis::K03Lowe);
    //  cout <<  "Knop03 (Low Extinction subsample): " << Knop03Chisq <<endl; 
    } 
    if (control.power_cdm && accept) {
      cout << "POWER WANTED" << endl;
      cmbcalc.dumpTransfer(&cosmos,control,transferFile);
	double zz = 0; bool print=true;
	//double Mass=1.e15*cosmos.h();
	double Mass=1.e15;
	cosmos.setPkBaryon(true);
	//cosmos.setPkBaryon(false);
	//char red[5]; sprintf(red, "%lf", zz);
	//string *namepower = "cdm" + red;
      //cosmos.dumpPower(0,"cdm");  // Spline dumps it, hence filename is cdm.dat
      string type = "tinker";
      //string type = "tinker_z";
  //  	cosmos.dumpPower(0,"cdm", 0, zz);  
//cosmos.growthFactor()->dump("growth_factor_lcdm");	
//cosmos.growthIndex()->dump("growth_index_lcdm");	
cosmos.dumpMassFunction(type, zz, print, 1.);  
cosmos.integrateComovingVolume(0.025,0.25);    
//cosmos.dumpNumberDensity("number_density", Mass);  
    	double z0 = 1; double z1 = 3;
	double th0 = 0; double th1 = 180;
	double ph0 = 0; double ph1 = 360;
	//double prob=cosmos.integrateComovingVolume(z0,z1);
	//double prob=cosmos.integrateOverComovingVolume(type, Mass, z0, z1, th0, th1, ph0, ph1);  
//cout << "Between redshift z0: " << z0 << " and z1: " << z1 << endl;
	//cout << "Above mass " << Mass/cosmos.h() << " h^-1 S.M., the total expected number of objects is: " << prob << endl;  
	//cout << "Above mass " << Mass << " h^-1 S.M., the total expected number of objects is: " << prob << endl;  
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
  cout << endl;
  cout << "done" << endl;
}

// ===================== END OF detailed() 


 

int main(int argc, char* argv[]) {
  switch (argc) {
  case 2:
    useFile(argv[1]);
    break;
  default:
    detailed();
  }
  return 0;
}
