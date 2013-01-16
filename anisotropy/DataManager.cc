/* file DataManager.cc
   function definitions for DataManager class
   written by Christian M. Mueller
   15.10.2003
*/
#include "DataManager.h"
#include <string>
#include <fstream>
#include <iostream>
#include <exception>
#include <stdlib.h>
#include <vector>
#include "gauge.h"
#include "quintcosmos.h"
#include "controlpanel.h"

/*! The only constructor of this class. All data values and the configuration file name are unspecified until set by 
  the appropriate functions
  \sa setConfigFileName(string name), synchronize()
*/
DataManager::DataManager() {
  _Riess04=_HST=_Tonry=_Astier05=false;
  _outputBase=ControlPanel::cmbeasyDir("/output/");
  _configFileName="";
  _leapingParameters.resize(3);
  _celestineParameters.resize(2);
  _crossoverParameters.resize(3);
  _crossoverfieldParameters.resize(5);
  _expParameters.resize(4);
}
 
/*! Sets the configuration file name to \a name. The file will not be opened by this function.
  \param name the name of the configuration file
 */
void DataManager::setConfigFileName(string name){  _configFileName=name;}

/*! This prints a formatted table of the data values stored in the private class variables. Will output random values if no 
  configuration file name has been set and \a synchronize has not been called. This does not change any variables.
  \sa setConfigFileName(string name), synchronize()
*/
void DataManager::printInfo(){
  cout << "----------------------------------------" << endl;
  cout << "----Begin configuration entries---------" << endl;
  cout << "----------------------------------------" << endl;

  cout << "Configuration file: " << _configFileName << endl;

  cout << "Gauge: ";
  if (_Gauge == Gauge::invariant) cout << "invariant" << endl;
  else if (_Gauge == Gauge::synchronous) cout << "synchronous" << endl;
  else if (_Gauge == Gauge::invariant2)  cout << "invariant2" << endl;
  else if (_Gauge == Gauge::quintSynchronous) cout << "quintSynchronous" << endl;
  else if (_Gauge == Gauge::quintInvariant) cout << "quintInvariant" << endl;
  else if (_Gauge == Gauge::speedyInvariant) cout << "speedyInvariant" << endl;
  else if (_Gauge == Gauge::speedyDEInvariant) cout << "speedyDEInvariant" << endl;
  else if (_Gauge == Gauge::coupledInvariant) cout << "coupledInvariant" << endl;


  cout << "InitialConditions:  " ;
  switch(_InitialConditions){
  case ControlPanel::adiabatic:
    cout << "adiabatic" << endl;
    break;
  case ControlPanel::isoCDM:
    cout << "isoCDM" << endl;
    break;
  case ControlPanel::isoBaryon:
    cout << "isoBaryon" << endl;
    break;
  case ControlPanel::isoNeutrino:
    cout << "isoNeutrino" << endl;
    break;
  case ControlPanel::mixed:
    cout << "mixed" << endl;
    cout << "---Contributions:" << endl;
    cout << "     Adiabatic:" << _AdiabaticContribution << endl;
    cout << "     IsoCDM:" << _IsoCDMContribution << endl;
    cout << "     IsoBaryon:" << _IsoBaryonContribution << endl;
    cout << "     IsoNeutrino:" << _IsoNeutrinoContribution << endl;
    break;
  default: 
    cout << "!!!Not specified!!!" << endl;
    break;
  }

  cout << "QuintessenceType: " ;
  switch(_QuintessenceType){
  case Quintessence::none:
    cout << "none" << endl;
    break;

  case Quintessence::exponential:  // Adding the possibility to set the value of the alpha parameter in the exponential from the config file
    cout << "exponential" << endl;
    cout << "--Parameters:" << endl;
    cout << "       Alpha:" << _expParameters[0] << endl ;
    cout << "   Q_Initial:" << _expParameters[1] << endl;
    cout << " Q_DotInital:" << _expParameters[2] << endl;
    cout << "          V0:" << _expParameters[3] << endl; 
   break;

  case Quintessence::leaping:
    cout << "leaping" << endl;
    cout << "--Parameters:" << endl;
    cout << "      k_min: " << _leapingParameters[0] << endl;
    cout << "      phi_0: " << _leapingParameters[1] << endl;
    cout << "      alpha: " << _leapingParameters[2] << endl;
    cout << "      omegals: " << _LeapingTune << endl;
    break;
  case Quintessence::ipl:
    cout << "ipl" << endl;
    cout << "--Parameters:" << endl;
    cout << "      alpha: " << _iplParameter << endl;
    break;
  case Quintessence::celestine:
    cout << "celestine" << endl;
    cout << "--Parameters:" << endl;
    cout << "      w0: " << _celestineParameters[0] << endl;
    cout << "      w1: " << _celestineParameters[1] << endl;
    break;
  case Quintessence::crossover:
    cout << "crossover" << endl;
    cout << "--Parameters:" << endl;
    cout << "      w0: " << _crossoverParameters[0] << endl;
    cout << "      wlsbar: " << _crossoverParameters[1] << endl;
    cout << "      A: "   << _crossoverParameters[2] << endl;
    break;
 case Quintessence::crossoverfield:
    cout << "crossoverfield" << endl;
    cout << "--Parameters:" << endl;
    cout << "      E: " << _crossoverfieldParameters[0] << endl;
    cout << "      J: " << _crossoverfieldParameters[1] << endl;
    cout << "      C: "   << _crossoverfieldParameters[2] << endl; 
    cout << "      D: "   << _crossoverfieldParameters[3] << endl; 
    cout << "      phi_crit: "   << _crossoverfieldParameters[4] << endl;
    break;   
  default: 
    cout << "!!!Not specified!!!" << endl;
    break;
  }

  cout << "-----Omega values today------" << endl;
  cout << "OmegaH2_b: " << _OmegaH2_b << endl;
  cout << "OmegaVacuum: " << _OmegaVacuum << endl;
  cout << "OmegaQuintessence: " << _OmegaQuintessence << endl;
  cout << "OmegaNonRelativisticNeutrinos: " << _OmegaNonRelativisticNeutrinos << endl;
  cout << "OmegaCDM (calculated from the other Omegas): " << 1.0 - _OmegaH2_b/(_HubbleH*_HubbleH) - _OmegaVacuum - _OmegaQuintessence - _OmegaNonRelativisticNeutrinos << endl;

  
  cout << "-----other cosmological parameters---"  << endl;
  cout << "HubbleH: " << _HubbleH << endl;
  cout << "CMBTemperature: " << _CMBTemperature << endl;
  cout << "HeAbundance: " << _HeAbundance << endl;
  cout << "RelativisticNeutrinoSpecies: " << _RelativisticNeutrinoSpecies << endl;
  cout << "NonRelativisticNeutrinoSpecies: " << _NonRelativisticNeutrinoSpecies << endl;

  // Coupled Neutrino

  cout << "-----flags and other stuff----" << endl;
  
  cout << "CMB: ";
  if(_CMB) cout << "true" << endl;
  else cout << "false" << endl; 
  
  cout << "CDMPowerSpectrum: ";
  if(_CDMPowerSpectrum) cout << "true" << endl;
  else cout << "false" << endl;
  
  cout << "TransferMaxK: " << _TransferMaxK << endl;
  
  cout << "HighPrecisionTransfer: ";
  if(_HighPrecisionTransfer) cout << "true" << endl;
  else cout << "false" << endl;
  
  cout << "TransferPerLog: " << _TransferPerLog << endl;
  cout << "CDMRedshift: " << _CDMRedshift << endl;


  cout << "Recombination: ";
  if(_Recombination == Recombination::Recfast) cout << "Recfast" << endl;
  else cout << "RecfastAlpha" << endl;

  cout << "Reionization:  " ;
  switch(_Reionization){
  case 0:
    cout << "0" << endl;
    break;
  case 1:
    cout << "1" << endl;
    break;
  case 2:
    cout << "2" << endl;
    break;
  default: 
    cout << "!!!Not specified!!!" << endl;
    break;
  }
  
  cout << "ReionizationZ: " << _ReionizationZ << endl;
  cout << "ReionizationFraction: " << _ReionizationFraction << endl;

  cout << "OptDistanceLss: " << _OptDistanceLss << endl;

  cout << "Scalar: ";
  if(_Scalar) cout << "true" << endl;
  else cout << "false" << endl;

  cout << "Tensor: ";
  if(_Tensor) cout << "true" << endl;
  else cout << "false" << endl;

  cout << "InitialPower: " << _InitialPower << endl;

  cout << "Normalizing to WMAP 1: " << (_wmapNormalize?"true":"false") << endl;
  cout << "Apply scalar/tensor ratio: " << (_ScalarTensorRatio?"true":"false") << endl;

  if (!_wmapNormalize) cout << "A_s: " << _A_s << endl;
  if (!_wmapNormalize && !_ScalarTensorRatio ) cout << "A_t: " << _A_t << endl;

  cout << "Lensflag: " ;
  switch(_Lensflag){
  case ControlPanel::none:
    cout << "none" << endl;
    break;
  case ControlPanel::linear:
    cout << "linear" << endl;
    break;
  case ControlPanel::nonelinear:
    cout << "nonlinear" << endl;
    break;
  default: 
    cout << "!!!Not specified!!!" << endl;
    break;
  }

  bool AnySNe = _Riess04 | _HST | _Tonry | _Astier05;
  if (AnySNe) {
    cout << "Supernova Data: true" << endl;
    if (_Riess04) cout << " Riess04 SNe\n";
    if (_HST) cout << " HST SNe\n";
    if (_Tonry) cout << " Tonry SNe\n";
    if (_Astier05) cout << " Astier05 (SNLS 1st year)\n";
  } else cout << "Supernova Data: false" << endl;

  cout << "AllSkyLensing: ";
  if (_AllSkyLensing) cout << "true" << endl;
  else cout << "false" << endl;

  cout << "-----file names--------" << endl;
  cout << "Output dir: " << _outputBase << endl;
  cout << "LensedFileName: " << _LensedFileName << endl;
  cout << "ScalarFileName: " << _ScalarFileName << endl;
  cout << "TensorFileName: " << _TensorFileName << endl;
  cout << "TransferFileName: " << _TransferFileName << endl; 
  cout << "BesselFileName: " << _BesselFileName << endl;


  cout << "--------------------------------------------" << endl;
  cout << "--------End of configuration entries--------" << endl;
  cout << "--------------------------------------------" << endl;
}  

/*!
  This opens the configuration file, reads the data and performs error checking on the obtained values. The usual sequence is 
  \a setConfigFileName(string name) and then \a synchronize()
 */
void DataManager::synchronize(){ // open configuration file and read it
    if (_configFileName ==  "") {
      cout << "No configuration file has been specified" <<endl;
      throw exception();   
    }              
    readConfigFile();
}

/*!
  Helper function called by \a synchronize(). Reads all lines of the configuration file and enters the values in the variables
  with error checking.
*/
void DataManager::readConfigFile(){
  ifstream inFile( _configFileName.c_str()); // The .c_str() returns a char * needed by ifstream constructor
    if (!inFile) {
      cout << "The configuration file could not be opened (wrong filename?)" << endl;
      throw exception();
    }

  string line;
  while(!inFile.eof()){
    getline(inFile,line);                               //read all lines until end of configuration file is reached
    if ( line[0]=='#' || line[0]==' ' || line[0]=='\n' ) {} /*Ignore lines beginning with empty space, # 
		  					    or newline character. */
    else{
      if (line.find(" ")== string::npos) {        
	                                             // if there are no empty spaces, we are done
      }
      else {
	line=line.substr(0,line.find_first_of(" "));  //remove everything after first empty space
      }
      
      enterLine(&line);                  //use helper function to sort into respective variables
    }
  }
}

/*!
  This function is called by \a readConfigFile() as a helper application
*/

void DataManager::enterLine(string *line){

  string pre=line->substr(0,line->find("="));   // part of the string in front of "="
  string post=line->substr(line->find("=")+1);  // part of the string after "="

// Unfortunately, switch does not work for strings, so we will use "else if"....

  if (pre == "Gauge") {	
    if (post == "invariant") _Gauge=Gauge::invariant;
    else if(post == "synchronous") _Gauge=Gauge::synchronous;
    else if(post == "invariant2") _Gauge=Gauge::invariant2;
    else if(post == "quintSynchronous") _Gauge=Gauge::quintSynchronous;
    else if(post == "quintInvariant")  _Gauge=Gauge::quintInvariant;    
    else if(post == "speedyInvariant")  _Gauge=Gauge::speedyInvariant;    
    else if(post == "speedyDEInvariant")  _Gauge=Gauge::speedyDEInvariant;    
#ifndef PRERELEASE
    else if(post == "coupledInvariant")  _Gauge=Gauge::coupledInvariant;    
#endif
    else {
      cout << "Invalid setting for " << pre << " in file " << _configFileName << endl;
      throw exception();
    }
  }

  else if(pre == "CMB") {
    if(post =="true" || post =="TRUE" || post =="True") _CMB=true;
    else if(post =="false" || post =="FALSE" || post =="False") _CMB=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }


  else if(pre == "CDMPowerSpectrum"){
    if(post =="true" || post =="TRUE" || post =="True") _CDMPowerSpectrum=true;
    else if(post =="false" || post =="FALSE" || post =="False") _CDMPowerSpectrum=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }

 else if(pre == "HighPrecisionTransfer"){
    if(post =="true" || post =="TRUE" || post =="True") _HighPrecisionTransfer=true;
    else if(post =="false" || post =="FALSE" || post =="False") _HighPrecisionTransfer=false; 
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }

  else if(pre == "TransferMaxK") _TransferMaxK=atoi(post.c_str());

  else if(pre == "TransferPerLog") _TransferPerLog=atoi(post.c_str());

  else if(pre == "CDMRedshift") _CDMRedshift=strtod(post.c_str(),'\0');

  else if(pre == "OutputBase") _outputBase=post;

  else if(pre == "TransferFileName") _TransferFileName=post;

  else if (pre == "Riess04") {
    if (post =="true" || post =="TRUE" || post =="True") _Riess04=true;
    else if (post =="false" || post =="FALSE" || post =="False") 
      _Riess04=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << 
	_configFileName << endl;
      throw exception();
    }
  }
  else if (pre == "HST") {
    if (post =="true" || post =="TRUE" || post =="True") _HST=true;
    else if (post =="false" || post =="FALSE" || post =="False") _HST=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << 
	_configFileName << endl;
      throw exception();
    }
  }
  else if (pre == "Tonry") {
    if (post =="true" || post =="TRUE" || post =="True") _Tonry=true;
    else if (post =="false" || post =="FALSE" || post =="False") _Tonry=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << 
	_configFileName << endl;
      throw exception();
    }
  }
  else if (pre == "Astier05") {
    if (post =="true" || post =="TRUE" || post =="True") _Astier05=true;
    else if (post =="false" || post =="FALSE" || post =="False") _Astier05=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << 
	_configFileName << endl;
      throw exception();
    }
  }

  else if(pre == "QuintessenceType") {
    if (post == "none") _QuintessenceType=Quintessence::none;
    else if (post == "exponential") _QuintessenceType=Quintessence::exponential;
    else if (post == "leaping") _QuintessenceType=Quintessence::leaping;
    else if (post == "ipl") _QuintessenceType=Quintessence::ipl;
    else if (post == "celestine") _QuintessenceType=Quintessence::celestine;
    else if (post == "crossover") _QuintessenceType=Quintessence::crossover;
    else if (post == "crossoverfield") _QuintessenceType=Quintessence::crossoverfield;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }

  else if(pre == "leaping:k_min" ) _leapingParameters[0]=strtod(post.c_str(),'\0');

  else if(pre == "leaping:phi_0" ) _leapingParameters[1]=strtod(post.c_str(),'\0');

  else if(pre == "leaping:alpha" ) _leapingParameters[2]=strtod(post.c_str(),'\0');
 
  else if(pre == "leaping:omegals" ) _LeapingTune=strtod(post.c_str(),'\0');

  else if (pre == "exponential:alpha_exp") _expParameters[0]=strtod(post.c_str(),'\0'); // exponential alpha parameter, V0 to be added 
  else if (pre == "exponential:qInitial") _expParameters[1]=strtod(post.c_str(),'\0');
  else if (pre == "exponential:qInitialDot") _expParameters[2]=strtod(post.c_str(),'\0');
  else if (pre == "exponential:V0") _expParameters[3]=strtod(post.c_str(),'\0');

  else if(pre == "ipl:alpha") _iplParameter=strtod(post.c_str(), '\0');

  else if(pre == "celestine:w0" ) _celestineParameters[0]=strtod(post.c_str(),'\0');
  else if(pre == "celestine:w1" ) _celestineParameters[1]=strtod(post.c_str(),'\0');

  else if(pre == "crossover:w0" ) _crossoverParameters[0]=strtod(post.c_str(),'\0');

  else if(pre == "crossover:wlsbar" ) _crossoverParameters[1]=strtod(post.c_str(),'\0');

  else if(pre == "crossover:A" ) _crossoverParameters[2]=strtod(post.c_str(),'\0');

  else if(pre == "crossoverfield:E" ) _crossoverfieldParameters[0]=strtod(post.c_str(),'\0');

  else if(pre == "crossoverfield:J" ) _crossoverfieldParameters[1]=strtod(post.c_str(),'\0');

  else if(pre == "crossoverfield:C" ) _crossoverfieldParameters[2]=strtod(post.c_str(),'\0');

  else if(pre == "crossoverfield:D" ) _crossoverfieldParameters[3]=strtod(post.c_str(),'\0');

  else if(pre == "crossoverfield:phi_crit" ) _crossoverfieldParameters[4]=strtod(post.c_str(),'\0');

  else if (pre == "CMBTemperature") _CMBTemperature=strtod(post.c_str(),'\0');

  else if (pre == "HeAbundance") _HeAbundance=strtod(post.c_str(),'\0');
  
  else if (pre == "RelativisticNeutrinoSpecies") _RelativisticNeutrinoSpecies=strtod(post.c_str(),'\0');

  else if (pre == "NonRelativisticNeutrinoSpecies") _NonRelativisticNeutrinoSpecies=strtod(post.c_str(),'\0');

  else if (pre == "BetaNeutrinoCoupling" ) _BetaNeutrinoCoupling=strtod(post.c_str(),'\0');

  else if (pre == "HubbleH") _HubbleH=strtod(post.c_str(),'\0' );

  else if(pre == "Recombination"){
    if (post == "recfast") _Recombination = Recombination::Recfast;
    else if (post == "recfastalpha") _Recombination = Recombination::RecfastAlpha;
   else {
     cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
     throw exception();
   }
  }

  else if (pre == "Reionization"){
    _Reionization=atoi(post.c_str());
    if (_Reionization>2 || _Reionization <0 ) {
      cout << "In configuration file: parameter ''Reionization'' can only take value 0,1 or 2" << endl;
      throw exception();
    }
  }

  else if (pre=="ReionizationZ")  _ReionizationZ=strtod(post.c_str(),'\0');

  else if (pre=="ReionizationFraction") _ReionizationFraction=strtod(post.c_str(),'\0');
  
  else if(pre == "Scalar"){
    if(post =="true" || post =="TRUE" || post =="True") _Scalar=true;
    else if(post =="false" || post =="FALSE" || post =="False") _Scalar=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }

  else if(pre == "Tensor"){
    if(post =="true" || post =="TRUE" || post =="True") _Tensor=true;
    else if(post =="false" || post =="FALSE" || post =="False") _Tensor=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }

  else if (pre == "InitialPower") _InitialPower=strtod(post.c_str(),'\0');

  else if(pre == "WmapNormalize") {
    if(post =="true" || post =="TRUE" || post =="True") _wmapNormalize=true;
    else if(post =="false" || post =="FALSE" || post =="False") _wmapNormalize=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }

  else if(pre == "ScalarTensorRatio") {
    if(post =="true" || post =="TRUE" || post =="True") _ScalarTensorRatio=true;
    else if(post =="false" || post =="FALSE" || post =="False") _ScalarTensorRatio=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }

  else if (pre == "A_s") _A_s=strtod(post.c_str(), '\0');

  else if (pre == "A_t") _A_t=strtod(post.c_str(), '\0');

  else if(pre == "Lensflag") {
    if (post == "none") _Lensflag=ControlPanel::none;
    else if (post == "linear") _Lensflag=ControlPanel::linear;
    else if (post == "nonlinear") _Lensflag=ControlPanel::nonelinear;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }

  else if(pre == "AllSkyLensing") {
    if(post =="true" || post =="TRUE" || post =="True") _AllSkyLensing=true;
    else if(post =="false" || post =="FALSE" || post =="False") _AllSkyLensing=false;
    else {
      cout << "Invalid setting for " << pre <<  " in file " << _configFileName << endl;
      throw exception();
    }
  }
  else if (pre == "BesselFileName") _BesselFileName=post;

  else if(pre == "ScalarFileName") _ScalarFileName=post;

  else if(pre == "TensorFileName") _TensorFileName=post;

  else if( pre== "LensedFileName") _LensedFileName=post;

  else if (pre == "InitialConditions") {
    if (post == "adiabatic") _InitialConditions=ControlPanel::adiabatic;
    else if (post =="isoCDM") _InitialConditions=ControlPanel::isoCDM;
    else if (post =="isoNeutrino") _InitialConditions=ControlPanel::isoNeutrino;
    else if (post =="mixed") _InitialConditions=ControlPanel::mixed;
    else if (post =="isoBaryon") _InitialConditions=ControlPanel::isoBaryon;
    else {
      cout << "Invalid setting for " << pre << " in file " << _configFileName << endl;
      throw exception();
    }
  }

  else if(pre == "AdiabaticContribution") _AdiabaticContribution=strtod(post.c_str(),'\0');
  else if(pre == "IsoCDMContribution") _IsoCDMContribution=strtod(post.c_str(),'\0');
  else if(pre == "IsoBaryonContribution") _IsoBaryonContribution=strtod(post.c_str(),'\0');
  else if(pre == "IsoNeutrinoContribution") _IsoNeutrinoContribution=strtod(post.c_str(),'\0');

  else if(pre == "OmegaH2_b") _OmegaH2_b=strtod(post.c_str(),'\0');

  else if(pre == "OmegaVacuum") _OmegaVacuum=strtod(post.c_str(),'\0');
 
  else if(pre == "OmegaQuintessence") _OmegaQuintessence=strtod(post.c_str(),'\0');
  
  else if(pre == "OmegaNonRelativisticNeutrinos") _OmegaNonRelativisticNeutrinos=strtod(post.c_str(),'\0');
  
  else if(pre == "OptDistanceLss") _OptDistanceLss=strtod(post.c_str(),'\0');
  
}


// the following functions merely give the corresponding values and are therefore quite simple

string DataManager::getConfigFileName() { return _configFileName;} 
string DataManager::getOutputBase() { return _outputBase;} 
Gauge::gauge DataManager::getGauge() { return _Gauge;}
bool DataManager::getCMB(){ return _CMB;}
bool DataManager::getCDMPowerSpectrum() {return _CDMPowerSpectrum; }
bool DataManager::getHighPrecisionTransfer(){ return _HighPrecisionTransfer;}
int DataManager::getTransferMaxK() {return _TransferMaxK;}
int DataManager::getTransferPerLog() {return _TransferPerLog;}
double DataManager::getCDMRedshift() {return _CDMRedshift;}
string DataManager::getTransferFileName() {return _TransferFileName;}
bool DataManager::getRiess04() const { return _Riess04; }
bool DataManager::getHST() const { return _HST; }
bool DataManager::getTonry() const { return _Tonry; }
bool DataManager::getAstier05() const { return _Astier05; }
Quintessence::Type DataManager::getQuintessenceType() {return _QuintessenceType;}

vector<double> DataManager::getQuintessenceParameters() {
  if (_QuintessenceType== Quintessence::leaping) {
    return _leapingParameters;
  }

  else if (_QuintessenceType==Quintessence::exponential) {  // EXP PARAMETERS RETURNED 
     return _expParameters; 
  }

  else if (_QuintessenceType==Quintessence::celestine){
    return _celestineParameters;
  }
  else if ( _QuintessenceType==Quintessence::ipl) {
    vector<double> paramVector(1);
    paramVector[0]=_iplParameter;
    return paramVector;
  }
  else if (_QuintessenceType==Quintessence::crossover){ 
    return _crossoverParameters;
  } 
  else if (_QuintessenceType==Quintessence::crossoverfield){ 
    return _crossoverfieldParameters;
  }
  else {
    vector<double> empty;
    return empty;
  }
}

double DataManager::getCMBTemperature() {return _CMBTemperature;}
double DataManager::getHeAbundance() {return _HeAbundance;}
double DataManager::getRelativisticNeutrinoSpecies() {return _RelativisticNeutrinoSpecies;}
double DataManager::getNonRelativisticNeutrinoSpecies() {return _NonRelativisticNeutrinoSpecies;}
double DataManager::getHubbleH() {return _HubbleH;}

double DataManager::getBetaCoupledNeutrino() {return _BetaNeutrinoCoupling;} // coupled neutrino get method implemented

Recombination::recombination DataManager::getRecombination() {return _Recombination;}
int DataManager::getReionization() {return _Reionization;}
double DataManager::getReionizationZ() { return _ReionizationZ;}
double DataManager::getReionizationFraction() {return _ReionizationFraction;}
bool DataManager::getScalar(){ return _Scalar;}
bool DataManager::getTensor() {return _Tensor;}
double DataManager::getInitialPower() {return _InitialPower;}
bool DataManager::getWmapNormalize() {return _wmapNormalize;}
bool DataManager::getScalarTensorRatio() {return _ScalarTensorRatio;}
double DataManager::getA_s() {return _A_s;}
double DataManager::getA_t() {return _A_t;}
int DataManager::getLensflag() {return _Lensflag;}
bool DataManager::getAllSkyLensing() {return _AllSkyLensing;}
string DataManager::getBesselFileName() {return _BesselFileName;}
string DataManager::getScalarFileName() {return _ScalarFileName;}
string DataManager::getTensorFileName() {return _TensorFileName;}
string DataManager::getLensedFileName() {return _LensedFileName;}
ControlPanel::icond DataManager::getInitialConditions() {return _InitialConditions;}
double DataManager::getAdiabaticContribution() {return _AdiabaticContribution;}
double DataManager::getIsoCDMContribution() {return _IsoCDMContribution;}
double DataManager::getIsoBaryonContribution(){return _IsoBaryonContribution;}
double DataManager::getIsoNeutrinoContribution(){return _IsoNeutrinoContribution;}
double DataManager::getOmegaH2_b() {return _OmegaH2_b;}
double DataManager::getOmegaVacuum() {return _OmegaVacuum;}
double DataManager::getOmegaQuintessence() {return _OmegaQuintessence; }
double DataManager::getOmegaNonRelativisticNeutrinos(){return _OmegaNonRelativisticNeutrinos;}
double DataManager::getOptDistanceLss(){ return _OptDistanceLss;}
double DataManager::getLeapingTune() {return _LeapingTune;}

