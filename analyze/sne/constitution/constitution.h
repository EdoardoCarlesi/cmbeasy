#ifndef CONSTITUTION_H
#define CONSTITUTION_H

#include "snedata.h"
#include "basecosmos.h"
#include "controlpanel.h"
#include "data.h"
#include "snedata.h"

#include <vector>
#include <vector>

class ConstitutionSNeData: public SNeData
{
  public:
    ConstitutionSNeData(): mInitialized(false) {}

    void init(bool includeSystematicErrors=true) {
      readData(ControlPanel::cmbeasyDir("/resources/sne1a/constitution/constitution.dat"), false);
      mInitialized = true;
    }
    double chi2(const baseCosmos& cosmos);

	// Ok this is very unelegant I think... creating a Data object inside a SNeData class.
	// But I found no other way to store these data - except rewriting the Sn1aCore routine
	// from scratch...
    Data *Constitution_09;
 
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

class Constitution09: public DataSet
{
  public:
  void initialize() {
    bool includeSysError=true;
    mSneConstitution.init(includeSysError);
    McTaskInfo::addEntry(McTaskInfo::LogLike,"ConstitutionSNeData09-loglike");
  }

  void computeLogLike(const McModel& model, McTaskInfo& result) {
    result("ConstitutionSNeData09-loglike") = -0.5*mSneConstitution.chi2(*model.cosmos());
	//cout << " Constitution SNeIa 2009 used a sample of 392 supernovae." << endl;
  }

 private:
    ConstitutionSNeData mSneConstitution;
};

REGISTER_DATASET(Constitution09)

#endif // MONTECARLO

#endif // UNION_H
