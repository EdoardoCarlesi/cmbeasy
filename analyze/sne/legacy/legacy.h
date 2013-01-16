#ifndef LEGACY_H
#define LEGACY_H

#include "snedata.h"
#include "basecosmos.h"
#include "controlpanel.h"
#include "data.h"
#include "snedata.h"

#include <vector>
#include <vector>

class LegacySNeData: public SNeData
{
  public:
    LegacySNeData(): mInitialized(false) {}

    void init(bool includeSystematicErrors=true) {
      readData(ControlPanel::cmbeasyDir("/resources/sne1a/legacy/snls_3rdyear_lcparams.txt"), false);
      mInitialized = true;
    }
    double chi2(const baseCosmos& cosmos);

	// Ok this is very unelegant I think... creating a Data object inside a SN class.
	// But I found no other way to store these data - except rewriting the Sn1aCore routine
	// from scratch...
    Data *Legacy;
 
  private:
    void   readData(const string& fileName, bool verbose);
    double  sum;
    typedef vector<double> Vec;
    typedef vector<Vec>    Matrix;
    Matrix  mInverseMatrix;
    bool    mInitialized;
};

#ifdef MONTECARLO
#include "mclikelihoodcalculator.h"

class LegacySN: public DataSet
{
  public:
  void initialize() {
    bool includeSysError=true;
    mSneLegacy.init(includeSysError);
    McTaskInfo::addEntry(McTaskInfo::LogLike,"LegacySN-loglike");
  }

  void computeLogLike(const McModel& model, McTaskInfo& result) {
    result("LegacySN-loglike") = -0.5*mSneLegacy.chi2(*model.cosmos());
	//cout << " Legacy SNeIa 2009 used a sample of 392 supernovae." << endl;
  }

 private:
    LegacySNeData mSneLegacy;
};

REGISTER_DATASET(LegacySN)

#endif // MONTECARLO

#endif // UNION_H
