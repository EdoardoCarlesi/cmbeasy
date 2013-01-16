#ifndef UNION2_H
#define UNION2_H

#include "snedata.h"
#include "basecosmos.h"
#include "controlpanel.h"

#include <vector>

class Union2SNeData: public SNeData
{
  public:
    Union2SNeData(): mInitialized(false) {}

    void init(bool includeSystematicErrors=true) {
      readData(ControlPanel::cmbeasyDir("/resources/sne1a/union2/sn_z_mu_dmu_union2.txt"), false);
      readCovMat(ControlPanel::cmbeasyDir("/resources/sne1a/union2/sn_covmat_sys_union2.txt"),
                 ControlPanel::cmbeasyDir("/resources/sne1a/union2/sn_covmat_nosys_union2.txt"),
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

class Union2: public DataSet
{
  public:
  void initialize() {
    bool includeSysError=true;
    mSneUnion2.init(includeSysError);
    McTaskInfo::addEntry(McTaskInfo::LogLike,"Union2SNeData-loglike");
  }

  void computeLogLike(const McModel& model, McTaskInfo& result) {
    result("Union2SNeData-loglike") = -0.5*mSneUnion2.chi2(*model.cosmos());
  }

 private:
    Union2SNeData mSneUnion2;
};

REGISTER_DATASET(Union2)

#endif // MONTECARLO

#endif // UNION_H
