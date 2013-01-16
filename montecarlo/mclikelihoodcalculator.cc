#include "mclikelihoodcalculator.h"

#include "cosmos.h"
#include "quintcosmos.h"
#include "vdecosmos.h"

#include "mcsettings.h"
#include "mcerror.h"

#include "sdsslrg.h"
#include "sdsslrg10.h"
#include "cbi2.h"

#ifndef NOWMAP7
#include "WMAP_7yr_util.h"
#endif

#include <algorithm>


std::list<DataSet*> McLikelihoodCalculator::mDataSets;


void McLikelihoodCalculator::useData(DataSet* data)
{
  mDataSets.push_back(data);
  data->initialize();
}

void McLikelihoodCalculator::initialize()
{
}

double McLikelihoodCalculator::computeLogLike(const McModel& model, McTaskInfo& result)
{
  list<DataSet*>::iterator it, end;
  it = mDataSets.begin();
  end = mDataSets.end();
  for ( ; it != end; ++it) {
    (*it)->computeLogLike(model, result);
    if (result.skipFurtherLikelihoodComputations)
      break;
  }
  return result.totalLogLike();
}

class Astier: public DataSet
{
  public:
    void initialize() {
      McTaskInfo::addMcParameter("Astier-alpha", 1., 2.2, 0.15);
      McTaskInfo::parameterInfo("Astier-alpha").isFast = true;
      McTaskInfo::addMcParameter("Astier-beta", 1.1, 2.4, 0.15);
      McTaskInfo::parameterInfo("Astier-beta").isFast = true;
      McTaskInfo::addEntry(McTaskInfo::LogLike, "Astier-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      double astier = mAnalyzeThis->SNIaAstier05(*(model.cosmos()), result("Astier-alpha"),
                                                                    result("Astier-beta"));
      result("Astier-loglike") = -0.5*astier;
    }
};
REGISTER_DATASET(Astier)


class SDSSLRG: public DataSet
{
  public:
    void initialize() {
      needTransfer();
      McTaskInfo::addMcParameter("SDSSLRG-bias", 0.2, 5., 0.4);
      McTaskInfo::parameterInfo("SDSSLRG-bias").isFast=true;
      McTaskInfo::addEntry(McTaskInfo::LogLike, "SDSSLRG-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      Cosmos& c = *model.cosmos();
      Spline* cdm = c.createPower(0, "cdm",  c.power_cdm(), 0, c.z2tau(0));
      cdm->arm();

      SdssLrg s(c, ControlPanel::cmbeasyDir("/resources/sdss_lrg/"));
      double bias = result("SDSSLRG-bias");
      double logLike = -0.5*s.chi2(*cdm, bias);
      delete cdm;
      result.setEntryValue("SDSSLRG-loglike", logLike);
    }
};
REGISTER_DATASET(SDSSLRG)


class SDSSLRG_10: public DataSet
{
  public:
    void initialize() {
/*
      needTransfer();
      McTaskInfo::addMcParameter("SDSSLRG_10-bias", 0.2, 5., 0.4);
      McTaskInfo::parameterInfo("SDSSLRG_10-bias").isFast=true;
      McTaskInfo::addEntry(McTaskInfo::LogLike, "SDSSLRG_10-loglike");
*/
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
/*
      Cosmos& c = *model.cosmos();
      AnalyzeThis& ai = *mAnalyzeThis;
      Spline* cdm = c.createPower(0, "cdm",  c.power_cdm(), 0, c.z2tau(0));
      //Spline *cdm =  c.k2Pk;
 //     cdm->arm();
      cdm->makeProper();
      //cdm = c.correctPower(c.Amplitude, 0, "power_cdm", cdm, 0, c.z2tau(c.z_pk));
      c.sigma8_z0=ai.sigma8(cdm, c.h());
      cout << "McLikelihoodCalculator::SDSSLRG_10 sigma8(): " << c.sigma8_z0 << endl;
 
      SdssLrg10 s(c, ControlPanel::cmbeasyDir("/resources/sdss_lrg10/"));
      double bias = result("SDSSLRG_10-bias");
      double logLike = -0.5*s.chi2(*cdm, bias);
      delete cdm;
      result.setEntryValue("SDSSLRG_10-loglike", logLike);
*/
    }
};
REGISTER_DATASET(SDSSLRG_10)





class WMAP7Data: public DataSet
{
  public:
#ifdef NOWMAP7
    void initialize() { throw McError("WMAP7Data::initialize() - compiled with -DNOWMAP7, but trying to use WMAP7 data."); }
    void computeLogLike(const McModel& model, McTaskInfo& result) {}
#else
  void initialize() {
    needCmb();
    McTaskInfo::addEntry(McTaskInfo::LogLike, "wmap7-totalLoglike");
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "wmap7-TTlike");
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "wmap7-TTlowllike");
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "wmap7-TTlowldet");
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "wmap7-Beamlike");
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "wmap7-TElike");
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "wmap7-TEdet");
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "wmap7-Lowllike");
    McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "wmap7-Lowldet");
  }


  void computeLogLike(const McModel& model, McTaskInfo& result) {
    CL& cl = *model.cmbSpectra();
    AnalyzeThis::WMAP7Likelihood like = mAnalyzeThis->WMAP7computeLikelihood(cl);

    // very last minute sanity check
    bool ok = (like.total>0.);

    if (wmap_util::wmap_likelihood_ok() && ok) {
      result("wmap7-totalLoglike") = -0.5*like.total;
      result("wmap7-TTlike") = -0.5*like.TTlike;
      result("wmap7-TTlowllike") = -0.5*like.TTlowllike;
      result("wmap7-TTlowldet") = -0.5*like.TTlowldet;
      result("wmap7-Beamlike") = -0.5*like.Beamlike;
      result("wmap7-TElike") = -0.5*like.TElike;
      result("wmap7-TEdet") = -0.5*like.TEdet;
      result("wmap7-Lowllike") = -0.5*like.Lowllike;
      result("wmap7-Lowldet") = -0.5*like.Lowldet;
    } else {
      static const double nan = numeric_limits<double>::quiet_NaN();
      result("wmap7-totalLoglike") = nan;
      result("wmap7-TTlike") = nan;
      result("wmap7-TTlowllike") = nan;
      result("wmap7-TTlowldet") = nan;
      result("wmap7-Beamlike") = nan;
      result("wmap7-TElike") = nan;
      result("wmap7-TEdet") = nan;
      result("wmap7-Lowllike") = nan;
      result("wmap7-Lowldet") = nan;
      result.setErrorLogLike();
      ofstream errLog("wmap_errorfile.dat", ios::app);
      errLog << noPrettyPrint << result << endl;
    }
  }
#endif
};
REGISTER_DATASET(WMAP7Data)

/* *************************************************************************************
 ** CBI and ACBAR, VSA, BOOMERANG03. First, we deselect several bands. Particularily
 ** The very high l bands. Hence, it sufficies to calculate up to l=2000,
 ** as the window functions would give almost zero weight to these high l's
 ************************************************************************************* */

class CBIData: public DataSet
{
  public:
    void initialize() {
      needCmb();
      McTaskInfo::addEntry(McTaskInfo::LogLike, "CBI-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      CL& cl = *model.cmbSpectra();
      CBI2 cbi2(cl);
      for (int i = 0; i < 5; i++)
        cbi2.bandDeselect[i] = true;
      for (int i = 11; i <= 15; i++)
        cbi2.bandDeselect[i] = true;
      result.setEntryValue("CBI-loglike", cbi2.chi2WithCalibration());
    }
};
REGISTER_DATASET(CBIData)

class AcbarData: public DataSet
{
  public:
    void initialize() {
      needCmb();
      McTaskInfo::addEntry(McTaskInfo::LogLike, "Acbar-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      AnalyzeThis& ai = *mAnalyzeThis;
      for (int i = 1; i < 5; i ++)
        ai.acbarBandDeselect[i] = true;
      for (int i = 12; i  < 15; i++)
        ai.acbarBandDeselect[i] = true;

      CL& cl = *model.cmbSpectra();
      double lnL = -0.5*ai.ACBARChi2WithCalibration(*cl.ts[0]);
      result("Acbar-loglike") = lnL;
    }
};
REGISTER_DATASET(AcbarData)

class VSAData: public DataSet
{
  public:
    void initialize() {
      needCmb();
      McTaskInfo::addEntry(McTaskInfo::LogLike, "VSA-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      AnalyzeThis& ai = *mAnalyzeThis;
      for (int i = 1; i  < 11; i++)
        ai.VSABandDeselect[i] = true;

      CL& cl = *model.cmbSpectra();
      double lnL = -0.5*ai.VSAChi2WithCalibration(*cl.ts[0]);
      result.setEntryValue("VSA-loglike", lnL);
    }
};
REGISTER_DATASET(VSAData)

class Boomerang: public DataSet
{
  public:
    void initialize() {
      needCmb();
      McTaskInfo::addEntry(McTaskInfo::LogLike, "Boomerang-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      AnalyzeThis& ai = *mAnalyzeThis;
      CL& cl = *model.cmbSpectra();

      //B-Mode deselect:
      ai.BOOMERANGDeselect(false, false, false, true);

      //for combining with WMAP consider only l_eff > 924
      for ( int i = 1; i <= 17; ++i )
        ai.BOOMERANGDeselect(i);

      ai.BOOMERANGInit(cl);
      double lnL = -0.5*ai.BOOMERANGChi2WithCalibration();
      result.setEntryValue("Boomerang-loglike", lnL);
    }
};
REGISTER_DATASET(Boomerang)


class TwoDF: public DataSet
{
  public:
    void initialize() {
      needTransfer();
      McTaskInfo::addEntry(McTaskInfo::LogLike, "2df-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      Cosmos& c = *model.cosmos();

      // obtain a powerspectrum at the redshift of the 2df-Survey (z=0.17)
      Spline *pwr = c.createPower(0, "TwoDF-2dfPower", c.power_cdm(), 0, c.z2tau(0.17));
      double bestbias = mAnalyzeThis->TwoDF_bestBias(pwr);
      *pwr *= bestbias * bestbias;
      pwr->arm();
      double loglike2df = -0.5*mAnalyzeThis->TwoDF_convolutePowerSpectrum(pwr);
      delete pwr;
      result("2df-loglike") = loglike2df;
    }
};
REGISTER_DATASET(TwoDF)

class SDSS: public DataSet
{
  public:
    void initialize() {
      needTransfer();
      McTaskInfo::addEntry(McTaskInfo::LogLike, "sdss-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
/*      AnalyzeThis& ai = *mAnalyzeThis;
      Cosmos& c = *model.cosmos();

      // obtain a powerspectrum at the redshift of the SDSS survey (z=0.1)
      // TODO changed because of memory leaks in SplineWeb
      //Spline *sdsspwr =  c.createPower(0, "SDSS", c.power_cdm(), 0, c.z2tau(0.10));
      Spline *sdsspwr =  c.k2Pk;
      //sdsspwr->arm();
      sdsspwr->makeProper();
      sdsspwr = c.correctPower(c.Amplitude, 0, "power_cdm", sdsspwr, 0, c.z2tau(c.z_pk));
      c.sigma8_z0=ai.sigma8(sdsspwr, c.h());
      cout << "McLikelihoodCalculator::SDSS sigma8(): " << c.sigma8_z0 << endl;
      double sdssbias = ai.SDSS_bestBias(sdsspwr);
      result("sdss-loglike") = -0.5*ai.SDSS_chiSquared(sdsspwr, sdssbias);
      delete sdsspwr;
*/
    }
};
REGISTER_DATASET(SDSS)

class Riess06: public DataSet
{
  public:
    void initialize() {
      McTaskInfo::addEntry(McTaskInfo::LogLike, "riess06-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      Cosmos& c = *model.cosmos();

      Spline lum(100, "Model::lum");
      for (double z=0; z<=2; z+=0.05) {
        lum.set(z, c.luminosityDistance(z));
      }
      lum.arm();
      result.setEntryValue("riess06-loglike", -0.5*mAnalyzeThis->Sn1aRiess06(lum));
    }
};
REGISTER_DATASET(Riess06)


/*! ***********************************************************************************
 ** SDSS Baryon coustic osciallation likelihood
 ** *********************************************************************************** */
class SdssBAO: public DataSet
{
  public:
    void initialize() {
      McTaskInfo::addEntry(McTaskInfo::LogLike, "sdss-bao-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      Cosmos& c = *model.cosmos();

      double ns = c.InitialPower[0];
      double ode = 0;
      QuintCosmos* qc = dynamic_cast<QuintCosmos*>(&c);
      if (qc) {  // if there's quintessence, take early dark energy into account
        ode = qc->omebar();
      }

      result("sdss-bao-loglike") = -0.5*mAnalyzeThis->SDSS_BAP_A_chiSquared(c, ns, ode);
    }
};
REGISTER_DATASET(SdssBAO)

/*! ***********************************************************************************
 ** SDSS Baryon coustic osciallation likelihood
 ** *********************************************************************************** */
class AngularBAO: public DataSet
{
  public:
    void initialize() {
      McTaskInfo::addEntry(McTaskInfo::LogLike, "angular-bao-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      Cosmos& c = *model.cosmos();

      double ns = c.InitialPower[0];
      double ode = 0;
      QuintCosmos* qc = dynamic_cast<QuintCosmos*>(&c);
      if (qc) {  // if there's quintessence, take early dark energy into account
        ode = qc->omebar();
      }

      VdeCosmos* vde = dynamic_cast<VdeCosmos*>(&c);
      if (vde) {  // if there's vector dark energy, take early vde into account
        ode = vde->omega_vde();
      }

      result("angular-bao-loglike") = -0.5*mAnalyzeThis->Angular_BAP_A_chiSquared(c, ns, ode);
    }
};
REGISTER_DATASET(AngularBAO)


/* *************************************************************************************
 ** Lyman Alpha Pat McDonald
 ************************************************************************************* */
class LyaMcDonald: public DataSet
{
  public:
    void initialize() {
      needTransfer();
      McTaskInfo::addEntry(McTaskInfo::LogLike, "Lya-loglike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      double chi2 = mAnalyzeThis->lymanAlphaPatMcDonaldChi2(*model.cosmos());
      result("Lya-loglike") =  -0.5*chi2;
    }
};
REGISTER_DATASET(LyaMcDonald)

