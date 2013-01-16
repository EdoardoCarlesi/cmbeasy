#ifndef UNION2_1_H
#define UNION2_1_H

#include "snedata.h"
#include "basecosmos.h"
#include "controlpanel.h"

#include <vector>

class Union21SNeData: public SNeData
{
  public:
    Union21SNeData(): mInitialized(false) {}

    void init(bool includeSystematicErrors=true) {
      readData(ControlPanel::cmbeasyDir("cmbeasy/resources/sne1a/union2.1/SCPUnion2.1_mu_vs_z.txt"), false);
      readCovMat(ControlPanel::cmbeasyDir("cmbeasy/resources/sne1a/union2.1/SCPUnion2.1_covmat_sys.txt"),
                 ControlPanel::cmbeasyDir("cmbeasy/resources/sne1a/union2.1/SCPUnion2.1_covmat_nosys.txt"),
                 includeSystematicErrors);
      mInitialized = true;
    }
    double chi2(const baseCosmos& cosmos);

 private:
    void   readData(const string& fileName, bool verbose);
    void   readCovMat(const string& fileNameCovMatSys, const string& fileNameCovMatNoSys, bool includeSysError);

    double  sum;
    typedef vector<double> Vec;
    typedef vector<Vec>    Matrix;
    Matrix  mInverseMatrix;
    bool    mInitialized;
};

#ifdef MONTECARLO
#include "mclikelihoodcalculator.h"

class Union21: public DataSet
{
  public:
  void initialize() {
    bool includeSysError=true;
    mSneUnion21.init(includeSysError);
    McTaskInfo::addEntry(McTaskInfo::LogLike,"Union21SNeData-loglike");
  }

  void computeLogLike(const McModel& model, McTaskInfo& result) {
    result("Union21SNeData-loglike") = -0.5*mSneUnion21.chi2(*model.cosmos());
  }

 private:
    Union21SNeData mSneUnion21;
};

REGISTER_DATASET(Union21)

#endif // MONTECARLO

#endif // UNION_H
