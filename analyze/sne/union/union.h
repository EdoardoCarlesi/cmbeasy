#ifndef UNION_H
#define UNION_H

#include "snedata.h"
#include "basecosmos.h"
#include "controlpanel.h"

#include <vector>

class UnionSNeData: public SNeData
{
  public:
    UnionSNeData(): mInitialized(false) {}

    void init(bool includeSystematicErrors=true) {
      readData(ControlPanel::cmbeasyDir("/resources/sne1a/union/SCPUnion_mu_vs_z.txt"), false);
      readCovMat(ControlPanel::cmbeasyDir("/resources/sne1a/union/SCPUnion_covmat_sys.txt"),
                 ControlPanel::cmbeasyDir("/resources/sne1a/union/SCPUnion_covmat_nosys.txt"),
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

class UnionSNeData08: public DataSet
{
  public:
  void initialize() {
    bool includeSysError=true;
    mSneUnion.init(includeSysError);
    McTaskInfo::addEntry(McTaskInfo::LogLike,"UnionSNeData08-loglike");
  }

  void computeLogLike(const McModel& model, McTaskInfo& result) {
    result("UnionSNeData08-loglike") = -0.5*mSneUnion.chi2(*model.cosmos());
  }

 private:
    UnionSNeData mSneUnion;
};

REGISTER_DATASET(UnionSNeData08)

#endif // MONTECARLO

#endif // UNION_H
