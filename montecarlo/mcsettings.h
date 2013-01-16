#ifndef MCSETTINGS_H
#define MCSETTINGS_H

#include "mcerror.h"

#include "controlpanel.h"
#include "analyzethis.h"

#include <string>

class McModel;
class McMaster;
class McRunner;
class DataSet;

class McSettings
{
  public:
    static McSettings* self();

    ControlPanel& controlPanel() { return mControlPanel; }
    AnalyzeThis&  analyzeThis() { return mAnalyzeThis; }

    McModel&      model() const;
    void          setModel(McModel *m) { mModel=m; }

    McMaster&     master() const;
    void          setMaster(McMaster *m) { mMaster=m; }

    McRunner&     slave() const;
    void          setSlave(McRunner *s) { mSlave=s; }

    typedef DataSet*(*DataSetGenerator)();
    typedef std::map<string, DataSetGenerator> GeneratorMap;
    static void registerDataSet(const std::string& name, DataSetGenerator create) {
      availableDataSets()[name] = create;
    }
    const std::string namesOfAvailableDataSets() const {
      std::string s;
      GeneratorMap::const_iterator it, end = availableDataSets().end();
      for ( it = availableDataSets().begin(); it != end; ++it) {
        s += it->first + ", ";
      }
      return s;
    }

    void useDataSet(const std::string& name);
    DataSet* dataSet(const std::string& name) const;

    std::string logDir() const { return mLogDir; }
    void        setLogDir(const std::string& d) { mLogDir = d; }
    std::string errorFileName() const;

    int numberOfChains() const { return mNumberOfChains; }
    void setNumberOfChains(int n) { mNumberOfChains = n; }

    //! minimum number of points for Gelman-Rubin R-statistic calculation
    unsigned int minRStatPoints() const { return mMinRStatPoints; }
    void setMinRStatPoints(unsigned int i) { mMinRStatPoints=i; }

    //! recalculate R-statistics every statisticsUpdateSteps()
    unsigned int statisticsUpdateSteps() const { return mStatisticsUpdateSteps; }
    void setStatisticsUpdateSteps(unsigned int s) { mStatisticsUpdateSteps=s; }

    bool estimateCovariance() const { return mEstimateCovariance; }
    void setEstimateCovariance(bool b) { mEstimateCovariance=b; }

    // settings for the adaptive stepsize algorithm
    //! use the adaptive stepsize algorithm
    bool adaptiveStepSize() const { return mAdaptiveStepSize; }
    void setAdaptiveStepSize(bool b) { mAdaptiveStepSize=b; }

    /*! If you want to stop adaptive stepsize after convergence has been reached;
      freeze-in will be performed once Gelman-Rubin statistics indicate convergence
      sets adaptiveStepSize() to false once converged.*/
    bool freezeIn() const { return mFreezeIn; }
    void setFreezeIn(bool b) { mFreezeIn=b; }

    /*! only if freezeIn() is true: If R-statistic for all parameters less than freezeInR() and
      number of points more than minFreezeInPoints() then freeze-in.
      If you want to be really conservative, you can call setFreezeInR(1.1), default is 1.2 */
    double freezeInR() const { return mFreezeInR; }
    void   setFreezeInR(double r) { mFreezeInR = r; }

    //! minimum size of chains for freeze-in
    unsigned int minFreezeInPoints() const { return mMinFreezeInPoints; }
    void         setMinFreezeInPoints(unsigned int p) { mMinFreezeInPoints=p; }


    /*! For the adaptive stepsize algorithm, variable step length; the defaults
     *  of 1.15/0.9 and 5/3 are rather arbitrary, but seem to work well.
     * Please note, that after freeze-in the covariance matrix without the
     * step factor will be used. Hence, the average stay at a certain parameter point
     * is approx 3
     * */
    double stepIncrease() const { return mStepIncrease; }
    void setStepIncrease(double d) { mStepIncrease = d; }
    double stepDecrease() const { return mStepDecrease; }
    void setStepDecrease(double d) { mStepDecrease = d; }
    int  highStepBound() const { return mHighStepBound; }
    void setHighStepBound(int h) { mHighStepBound=h; }
    int lowStepBound() const { return mLowStepBound; }
    void setLowStepBound(int l) { mLowStepBound=l; }
    //*! For the adaptive stepsize algorithm: when true, ignore out-of-bound steps
    //*  for the step count when determining whether to decrease the stepsize. Default is true.
    //* */
    void setUseReallyInvestigated(bool b) { mUseReallyInvestigated=b; }
    bool useReallyInvestigated() const { return mUseReallyInvestigated; }


    bool useAverageEntireFactor() const { return mUseAverageEntireFactor; }
    void setUseAverageEntireFactor(bool b) { mUseAverageEntireFactor=b; }

    //! For the covariance updater; recompute covariance matrix after updateTime() steps
    int updateTime() const { return mUpdateTime; }
    void setUpdateTime(int t) { mUpdateTime=t; }

    /*! Start using the estimated covariance matrix for steps after beginCovUpdate() points have
     * been computed in a chain
     * */
    unsigned int beginCovUpdate() const { return mBeginCovUpdate; }
    void setBeginCovUpdate(unsigned int i) { mBeginCovUpdate=i; }

   /*! For the adaptive stepsize: Do not compute covariance Matrix
    *  with more than maxPreviousCovPoints() */
    unsigned int maxPreviousCovPoints() const { return mMaxPreviousCovPoints; }
    void setMaxPreviousCovPoints(unsigned int i) { mMaxPreviousCovPoints=i; }

    /*! when using FastSlowStepper(): the number of fast steps to take
     *  per slow step */
    unsigned int fastStepRatio() const { return mFastStepRatio; }
    void setFastStepRatio(unsigned int i) { mFastStepRatio=i; }

  private:
    McSettings();
    McSettings(const McSettings&);

  private:
    ControlPanel mControlPanel;
    AnalyzeThis  mAnalyzeThis;
    std::string  mLogDir;
    int          mNumberOfChains;
    unsigned int mMinRStatPoints;
    unsigned int mStatisticsUpdateSteps;
    bool         mEstimateCovariance;
    bool         mAdaptiveStepSize;
    bool         mFreezeIn;
    double       mFreezeInR;
    unsigned int mMinFreezeInPoints;
    double       mStepIncrease;
    double       mStepDecrease;
    int          mHighStepBound;
    int          mLowStepBound;
    int          mUseReallyInvestigated;
    bool         mUseAverageEntireFactor;
    int          mUpdateTime;
    unsigned int mBeginCovUpdate;
    unsigned int mMaxPreviousCovPoints;
    unsigned int mFastStepRatio;
    McModel*     mModel;
    McMaster*    mMaster;
    McRunner*    mSlave;

    static GeneratorMap& availableDataSets();
    std::map<string, DataSet*> mActiveDataSets;
};


#endif // MCSETTINGS_H
