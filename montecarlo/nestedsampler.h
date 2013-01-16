#ifndef NESTEDSAMPLING_H
#define NESTEDSAMPLING_H

#include "mcrunner.h"
#include "mcmodel.h"

#include <deque>

class NestedSampler: public McMaster
{
  public:
    NestedSampler();
    virtual ~NestedSampler();

    virtual void init(bool restart=false);
    virtual void runMaster();

    enum SamplingMethod { UniformOverPrior,
                          EllipsoidalSampling,
                          MultiNest
    };
    void setSamplingMethod(const SamplingMethod m) { mSamplingMode=m; }
    void setLivePointsCount(unsigned int c) { mLivePointsCount=c; }
    void setMinimumIterations(double i) { mMinIterations=i; }
    void setMaxRemainingLogEvidenceFraction(double f) { mMaxRemainingLogEvidenceFraction=f; }
    void setEllipsoidalEnlargementFactor(double f) { mEllipsoidEnlargementFactor=f; }

    void readLivePoints();
    void produceCurrentPosterior();

  protected:
    void updatePoints(McTaskInfo& result);
    void updateEvidence();
    void completeEvidence();
    McTaskInfo nextParameterSet();
    McTaskInfo nextParameterSetEllipsoidal();
    //McTaskInfo nextParameterSetMultiNest();

    bool terminateCondition() const;
    void writePosteriorSamples();

    double information() const { return mInformation; }
    double logEvidence() const { return mLogEvidence; }

  protected:
    unsigned int mLivePointsCount;
    std::vector<McTaskInfo> mLivePoints;

    typedef struct {
      McTaskInfo point;
      double logWeight;
    } DiscardedPoint;
    friend struct LogWeightComp;
    std::vector<DiscardedPoint> mDiscardedPoints;

    std::vector<std::vector<double> > mInitialPoints;

    double mBestLiveLnLike, mCurrentLogWeight;
    double mCurrentIteration, mCurrentLogX_i;
    double mInformation, mLogEvidence;
    unsigned long int mPointsComputed, mPointsUsed, mErrorCount;
    double mEllipsoidEnlargementFactor;

    SamplingMethod mSamplingMode;
    bool           mWasRestarted;

    double mMinIterations, mMaxRemainingLogEvidenceFraction;
};

#endif // NESTEDSAMPLING_H
