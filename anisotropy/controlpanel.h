#ifndef CONTROLPANEL_H 
#define CONTROLPANEL_H

#include "safevector.h"
#include "recombination.h"
#include "recfast.h"

#include <string.h>

/*! small helper class to manipulate initial conditions
 *
 *  @warning: not yet implemented for all of the perturbation/gauge classes
 */
class InitialConditionFactors
{
  public:

    InitialConditionFactors()
      : DgpFactor(1), DgnFactor(1), DgbFactor(1),
        DgcFactor(1), VpFactor(1), VbFactor(1),
        VnFactor(1), VcFactor(1)
    {
    }

    double DgpFactor, DgnFactor, DgbFactor, DgcFactor,
     VpFactor, VbFactor, VnFactor, VcFactor;

    void apply(double& Dgp, double& Dgn, double& Dgb, double& Dgc,
               double& Vp,  double& Vn,  double& Vb, double& Vc) const {
      Dgp *= DgpFactor;
      Dgn *= DgnFactor;
      Dgb *= DgbFactor;
      Dgc *= DgcFactor;
      Vp *= VpFactor;
      Vb *= VpFactor;
      Vn *= VnFactor;
    }
};


/*!
  Miniature struct to tell cmbcalc and friends what it should 
  do, i.e. transfer, CMB, which recombination and 
  so on. replaces cryptic transfer.itf = 0,1,2 etc
  
  The logic structure of the flags cmb, power_cdm, scalar and tensors 
  is as follows:

  cmb: cmb anisotropies wanted ?
  power_cdm: cdm power spectrum / transfer functions wanted ?
  
  scalar: calculate scalar fluctuations. Needed for scalar Cl's and power_cdm
  tensor: calculate tensor fluctuations. Needed for tensor Cl's

  Hence, you will get:

  Scalar Cl's:   cmb = scalar = true
  Tensor Cl's:  cmb = tensor = true
  CDM only:  cmb = false,  power_cdm = scalar = true

*/
struct ControlPanel {
  enum icond { adiabatic, isoCDM, isoBaryon, isoSeed,isoNeutrino,mixed};
  InitialConditionFactors initialConditionFactors;
  enum lens  {none, linear, nonelinear};
  bool power_cdm;  //!< if true, fill cosmos->power_cdm SplineWeb with the powerspectrum
  bool cmb; //!< if true, calculate cmb spectrum
  bool verbose; //!< if true, dump several splines to harddisk and inform about cosmological quantities in terminal

  double transferMaxK;  //!< Maximum k-mode in h/Mpc  (akmaxt)
  double minK; //!< if set, gives minimum k (otherwise calculated automatically)
  int transferPerLog; //!< Number of transfer functions per logarithmic interval 
  bool IsPhantomCrossing; //!< true, if w of dark energy might possibly cross -1
  SafeVector<double> transferZ; //!< Holds the z values for each equested transfer output. This is legacy stuff :-) (ntf, ztf)

  bool scalar;  //!< scalar fluctuations on/off
  bool tensor; //!< tensor fluctuations on/off
  //  bool canonicalTensorRatio; //!< Tensor/Scalar ratio 7(1-n_s) yes / no ??? OBSOLETE
  //  bool originalTensorHandling; //!< Use cmbfast original Tensor/Scalar handling,  OBSOLETE

  bool highPrecisionTransfer; //!< true if high precision for transfer functions at small scales requested
  bool allSkyLensing; //!< true when using the all-sky correlation functions method from astro-ph/0502425
  icond initialCond; 
  double adiabaticContribution;
  double isoCDMContribution;
  double isoBaryonContribution;
  double isoNeutrinoContribution;
  //static vector<string> IcondHuman(5) =  {"none","adiabatic","isoCDM","isoBaryon","isoSeed"};
  Recombination::recombination recombination;  //!< Choose the kind of recombination. Currently Recfast or RecfastAlpha
  Recfast::HeliumRecombinationSwitch recfastHeliumRecombination; //!< integer option (0..6) for recfast HeI recombination; see recfast.h for details 
  unsigned int recfastHSwitch;   //!< either 0 or 1; Recfast modification for H correction: 0) no change from old Recfast, 1) include correction

  lens lensing;

  /*! This variable may be used to look at a specif k-mode, if for
    instance OUTPUTFILE is defined in Invariant our
    Synchronous. Another possibility would be to leave stop_pert alone
    and create a SplineWeb much the same as dWeb etc. To store the
    quantities you are interested in and then use SplineWeb->dump() to look
   at its evolution.*/
  double StopPert;  //!< k-value, at which CmbCalc should terminate. 

  /*! Creation and initialization of Controlpanel. Please observe that
    by default, StopPert is very large, thus CmbCalc will not interrupt its
    calculation, which is good so :-)
  */
  ControlPanel() : power_cdm(false), cmb(true), verbose(false), transferMaxK(5), minK(0), transferPerLog(5), IsPhantomCrossing(false), scalar(true), tensor(false) , highPrecisionTransfer(false), allSkyLensing(true), initialCond(adiabatic), 
     recombination(Recombination::Recfast),StopPert(1000) {
    transferZ[0] = 0;
    setLensing(none); //default 
    recfastHeliumRecombination=6;
    recfastHSwitch=1;
  }

  virtual ~ControlPanel() {}

  void setInitialConditions(icond i) {  initialCond = i; }
  void setInitialConditions(int i) { 
    switch (i) {
    case 2: 
       initialCond= isoCDM; break;
    case 3:
       initialCond= isoBaryon; break;
    case 4:
       initialCond= isoSeed; break;
    default:
       initialCond = adiabatic; break;
    }
  }
  void setLensing(int i) {
    switch (i) {
    case 1:
      lensing=linear; break;
    case 2:
      lensing=nonelinear;break;
    default:
      lensing=none;break;
    }
  }
  bool isLensing() { return lensing != none; }
  void setAllSkyLensing(bool b) { allSkyLensing = b; }
  bool isAllSkyLensing() { return allSkyLensing; }
  /*! Read CMBEASY environment variable and return
    path to this directory. Returns string="", if no variable is 
  defined. If throwerror is true, then a Bad_Error will be thrown,
  in case the environment variable is not defined. */
  static string cmbeasyDir(string addon="", bool throwerror=true) { 
    const char *cmbdir = getenv("CMBEASYDIR");
    if (!cmbdir) {
      if (! throwerror) {
        return "";
      }
      throw Bad_Error("No CMBEASYDIR environment variable defined. Please define one in your .bashrc or  the appropriate resource file, if you use a different shell");
    }
    return string(cmbdir) + addon;
  }

  bool isPhantomCrossing() const { return IsPhantomCrossing; } //!< might dark energy w crossing -1?
  void setPhantomCrossing(bool p) { IsPhantomCrossing=p;}  //!< set dark energy w crossing possible, i.e. treat as fluid with c_s^2 = 1 on all scales


  void printStatus() {
    cout << endl << "######### ControlPanel -- Status ######" << endl;
    cout << "Initial Condition: " <<initialCond << endl;

  }
};
#endif
