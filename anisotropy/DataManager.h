/* file DataManager.h
   Class for reading values from a configuration file
   written by Christian M. Mueller
   15.10.2003
*/
#ifndef DATAMANAGER_H
#define DATAMANAGER_H
#include <string>
#include <vector>
#include "gauge.h"
#include "quintcosmos.h"
#include "controlpanel.h"

/*! \brief Class for reading data from a configuration file.

  This Class reads data from a configuration file, converts it and
  performs some minimal error checking. Exceptions are thrown if a data 
  value could not be read. All getXXX functions merely return the stored value.
*/

class DataManager {
 public:
  DataManager();

  void setConfigFileName(string name); 
  void synchronize(); /* no synronization if filename has not been specified!
			 This reads the data values from the file into the 
			 private member variables */
  void printInfo();  // output of member value variables 
  
  string getConfigFileName(); 
  string getOutputBase(); 
  Gauge::gauge getGauge();
  bool getCMB();
  bool getCDMPowerSpectrum() ;
  bool getHighPrecisionTransfer();
  int getTransferMaxK();
  int getTransferPerLog();
  double getCDMRedshift() ;
  string getTransferFileName() ;
  bool getRiess04() const;
  bool getHST() const;
  bool getTonry() const;
  bool getAstier05() const;
  Quintessence::Type getQuintessenceType();
  vector<double> getQuintessenceParameters();
  double getCMBTemperature();
  double getHeAbundance();
  double getRelativisticNeutrinoSpecies();
  double getNonRelativisticNeutrinoSpecies();
  double getBetaCoupledNeutrino(); // get Beta parameter value 
  double getHubbleH();



  Recombination::recombination getRecombination() ;
  int getReionization() ;
  double getReionizationZ();
  double getReionizationFraction();
  bool getScalar() ;
  bool getTensor() ;
  double getInitialPower() ;
  int getLensflag();
  bool getAllSkyLensing();
  string getBesselFileName();
  string getScalarFileName() ;
  string getTensorFileName();
  string getLensedFileName();
  ControlPanel::icond getInitialConditions() ;
  double getAdiabaticContribution();
  double getIsoCDMContribution();
  double getIsoBaryonContribution();
  double getIsoNeutrinoContribution();
  double getOmegaH2_b() ;
  double getOmegaVacuum() ;
  double getOmegaQuintessence(); 
  double getOmegaNonRelativisticNeutrinos();
  double getOptDistanceLss() ;
  bool   getWmapNormalize();
  bool   getScalarTensorRatio();
  double getA_s();
  double getA_t();
  double getLeapingTune();

 private:
  void readConfigFile();
  void enterLine(string *line);

  string _configFileName;
  string _outputBase;
  Gauge::gauge _Gauge;
  bool _CMB;
  bool _CDMPowerSpectrum;
  bool _HighPrecisionTransfer;
  int _TransferMaxK;
  int _TransferPerLog;
  double _CDMRedshift;
  string _TransferFileName;

  //quintessence parameters

  Quintessence::Type _QuintessenceType;
  vector<double> _leapingParameters;
  vector<double> _crossoverParameters;
  vector<double> _crossoverfieldParameters;
  vector<double> _celestineParameters;
 
  vector<double> _expParameters; // setting the exponential parameters into vectordouble type, to allow multiple parameter cfg.file choice

  double _iplParameter;

  //---------------------//

  double _CMBTemperature;
  double _HeAbundance;
  double _RelativisticNeutrinoSpecies;
  double _NonRelativisticNeutrinoSpecies;
  double _BetaNeutrinoCoupling; // CN coupling
  double _HubbleH;



  Recombination::recombination  _Recombination;
  int _Reionization;
  double _ReionizationZ;
  double _ReionizationFraction;
  bool _Scalar;
  bool _Tensor;
  bool _wmapNormalize;
  bool _ScalarTensorRatio;
  double _A_s;
  double _A_t;
  double _InitialPower;
  bool _TensorToScalarQuadrupoleQ;
  int _Lensflag;
  bool _AllSkyLensing;
  string _ScalarFileName;
  string _BesselFileName;
  string _TensorFileName;
  string _LensedFileName;
  ControlPanel::icond _InitialConditions;
  double _OmegaH2_b;
  double _OmegaVacuum;
  double _OmegaQuintessence;
  double _OmegaNonRelativisticNeutrinos;
  double _OptDistanceLss;
  double _AdiabaticContribution;
  double _IsoCDMContribution;
  double _IsoBaryonContribution;
  double _IsoNeutrinoContribution;
  double _LeapingTune;

  bool _Riess04;
  bool _HST;
  bool _Tonry;
  bool _Astier05;

};


#endif
