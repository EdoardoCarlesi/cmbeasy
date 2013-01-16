// vim:ts=2:sw=2:et

#include "fastslowstepper.h"
#include "mcmodel.h"

#include "mclikelihoodcalculator.h"
#include "multigaussian.h"

#include <iterator>
#include <iostream>
#include <deque>
#include <algorithm>

using namespace std;



FastSlowStepper::FastSlowStepper() : McStandardSampler()
{
}

FastSlowStepper::~FastSlowStepper()
{
  vector<MultiGaussian*>::iterator it, end;
  it = mFastGauss.begin();
  end = mFastGauss.end();
  while (it != end) {
    delete *it;
    ++it;
  }
}

void FastSlowStepper::setupChains(bool restart)
{
  if (restart) {
    throw McError("FastSlowStepper::setupChains() - restart not implemented");
  }
  McStandardSampler::setupChains(restart);

  vector<double> lowerBounds, upperBounds, initialSigmas;

  int fastParameterCount = 0;
  const McTaskInfo::ParameterInfoMap& infos=McTaskInfo::parameterInfos();
  for (int i=0; i<McTaskInfo::paramCount(); ++i) {
    const McTaskInfo::ParameterInfo& curInfo = infos.find(McTaskInfo::parameterName(i))->second;
    mIsFastParam.push_back(curInfo.isFast);
    if (curInfo.isFast) {
      ++fastParameterCount;
      lowerBounds.push_back(curInfo.lowerBound);
      upperBounds.push_back(curInfo.upperBound);
      initialSigmas.push_back(curInfo.initialSigma);
      cout << "added fast parameter: " << McTaskInfo::parameterName(i) << endl;
    }
  }
  if (fastParameterCount<=0) {
    throw McError("FastSlowStepper::setupChains() - using FastSlowStepper without any fast parameters.");
  }

  mFastCovariance.resize(numberOfChains(), CovMatrix(fastParameterCount,
                                    vector<double>(fastParameterCount, 0.)));

  for (unsigned int i = 0; i < numberOfChains(); i++) {
    McChain& chain = *mChain[i];

    MultiGaussian *mg = new MultiGaussian(fastParameterCount);

    mg->setBounds(lowerBounds, upperBounds);
    mFastGauss.push_back(mg);

    for(unsigned int k=0; k<fastParameterCount; k++) {
      for(unsigned int j=0; j<fastParameterCount; j++) {
        if (j==k)
          mFastCovariance[i][j][j]=initialSigmas[j]*initialSigmas[j];
        else
          mFastCovariance[i][j][k]=0.;
      }
    }
    mFastGauss[i]->generateEigenvectors(mFastCovariance[i], 1.0);
    //mFastGauss[i]->printExtInfo(cout);
  }


  mFastStepsTaken.resize(fastParameterCount, 0);
  mLastProposalWasFast.resize(fastParameterCount, false);
}

void FastSlowStepper::stepTaken(McChain& chain, McTaskInfo& step)
{
  const int id = chain.id;
  if (mLastProposalWasFast[id]) {
    ++mFastStepsTaken[id];
    cout << "took fast step"  << endl;
  } else {
    mFastStepsTaken[id]=0;
    cout << "took slow step"  << endl;
  }

  McStandardSampler::stepTaken(chain, step);
}

static bool throwDice(const vector<double>& oldPoint, vector<double>& newPoint, MultiGaussian *gauss)
{
  bool success = gauss->throwDice(oldPoint);
  newPoint=gauss->getRandomValues();
  return success;
}

void FastSlowStepper::computeNextStepForChain(McChain& chain)
{
  if ((chain.size<=1)  || mConverged
      || (mFastStepsTaken[chain.id]>=cfg().fastStepRatio())) {
    cout << "propose slow step:\n";
    mLastProposalWasFast[chain.id]=false;
    McStandardSampler::computeNextStepForChain(chain);
    return;
  }

  mLastProposalWasFast[chain.id]=true;
  McTaskInfo& lastPoint = chain.points.back();

  cout << "propose fast step:\n";
  //cout << prettyPrint << lastPoint << endl;
  vector<double> allOldParameters,
                 oldPoint, newPoint; // fast parameters only
  allOldParameters = lastPoint.parameterVector();
  for (int i=0; i<allOldParameters.size(); ++i) {
    if (mIsFastParam[i]) {
      oldPoint.push_back(allOldParameters[i]);
    }
  }

  while (!throwDice(oldPoint, newPoint, mFastGauss[chain.id])) {
    lastPoint.Multiplicity++;
  }

  chain.nextPoint = lastPoint;
  int idx=0;
  for (int i=0; i<allOldParameters.size(); ++i) {
    if (mIsFastParam[i]) {
      chain.nextPoint.setParameterValue(i, newPoint[idx++]);
    }
  }
  //cout << "new point: " << prettyPrint << chain.nextPoint << endl;
}


void FastSlowStepper::adjustAdaptiveStepSize(McChain& chain, bool tookStep)
{
  McStandardSampler::adjustAdaptiveStepSize(chain, tookStep);
  if ((tookStep && chain.stepsSinceUpdate!=0) || !tookStep || !cfg().adaptiveStepSize()) {
    // the full covariance matrix was not updated,  so we keep the old one, too
    return;
  }
  const unsigned int paramCount = McTaskInfo::paramCount();
  unsigned int k,l;
  k=l=0;
  for (int i=0; i<paramCount; ++i) {
    if (!mIsFastParam[i])
      continue;
    for (int j=0; j<paramCount; ++j) {
      if (!mIsFastParam[j])
        continue;
      mFastCovariance[chain.id][k][l++] = chain.covMatrix[i][j];
    }
    k++;
  }
  mFastGauss[chain.id]->generateEigenvectors(mFastCovariance[chain.id], chain.EntireFactor);
}

class ModelCache
{
  public:
    struct CacheEntry {
      McTaskInfo params;
      Cosmos* c;
      CL* rawCls;
    };

  private:

    struct SameCosmosPointer {
      Cosmos* mC;
      SameCosmosPointer(Cosmos* c): mC(c) {}
      bool operator()(const CacheEntry& e) { return e.c==mC; }
    };

    struct SameParams {
      const bool mSlowOnly;
      const vector<double> mPV;

      SameParams(const McTaskInfo& p, bool b)
        : mPV(p.parameterVector()), mSlowOnly(b) {}

      bool operator()(const CacheEntry& p) {
        vector<double> pV = p.params.parameterVector();
        if (!mSlowOnly) {
          return mPV==pV;
        }
        const McTaskInfo::ParameterInfoMap& info = McTaskInfo::parameterInfos();
        for (int i=0; i<pV.size(); ++i) {
          const string& name =p.params.parameterName(i);
          if (!info.find(name)->second.isFast && (pV[i]!=mPV[i])) {
            return false;
          }
        }
        return true;
      }
    };

  public:
    ModelCache(unsigned int size): mMaxSize(size), mSize(0) {}

    void insert(const McTaskInfo& params, McModel& m) {
      CacheItem* i = dynamic_cast<CacheItem*>(&m);
      if (!i) {
        throw McError("ModelCache::insert() used with a model that doesn't"
                      "implement the CacheItem interface; do not inherit from McCustomModel<>,"
                      "but from CacheableModel<> to prevent this error.");
      }
      CacheEntry e = { params, i->stealCosmos(), i->stealRawCmbSpectra() };
      mCache.push_back(e);
      mSize++;
      checkSize();
    }

    enum CompareMode { CompareSlowParamsOnly, CompareAllParams };

    CacheEntry* findModel(McTaskInfo& params, CompareMode m = CompareAllParams) {
      deque<CacheEntry>::iterator it;
      it = find_if (mCache.begin(), mCache.end(),
                    SameParams(params, m==CompareSlowParamsOnly));
      if (it==mCache.end()) {
        return 0;
      }
      return &(*it);
    }

    static CacheItem* setDataFromCache(McModel& m, CacheEntry* e) {
      CacheItem* i = dynamic_cast<CacheItem*>(&m);
      if (!i) {
        throw McError("ModelCache::setDataFromCache() used with a model that doesn't"
                      "implement the CacheItem interface; do not inherit from McCustomModel<>,"
                      "but from CacheableModel<> to prevent this error.");
      }
      i->setRawCmbSpectra(e->rawCls->clone());
      i->setCosmos(e->c);
      return i;
    }

  private:
    void checkSize() {
      if (mSize<=mMaxSize)
        return;
      if (count_if(mCache.begin(), mCache.end(), SameCosmosPointer(mCache.front().c))<=1) {
        // the cosmos objects are not copied, so we only delete them if no
        // other CacheEntry points to the same Cosmos
        delete mCache.front().c;
      }
      // but every CacheEntry has its own copy of the Cls, so we delete them.
      delete mCache.front().rawCls;
      mCache.pop_front();
      --mSize;
    }

  private:
    deque<CacheEntry> mCache;
    unsigned int mMaxSize, mSize;
};

void FastSlowSlave::runSlave()
{
  McModel& model = cfg().model();
  McLikelihoodCalculator calcLike;
  bool quitRequested;

  ModelCache cache(1/*size*/);
  ModelCache::CacheEntry* e;
  CacheItem* cachedModel;

  while (true) {
    McTaskInfo params = parametersFromMaster(&quitRequested);
    if (quitRequested) {
      return;
    }
    e = cache.findModel(params, ModelCache::CompareSlowParamsOnly);
    bool cachedResultsAvailable = (e!=0);
    try {
      model.setParameters(params);
      if (cachedResultsAvailable) {
        cachedModel=ModelCache::setDataFromCache(model, e);
        cachedModel->recomputeFast();
        calcLike.computeLogLike(model, params);
      } else {
        model.compute();
        calcLike.computeLogLike(model, params);
        cache.insert(params, model);
      }
    } catch (Bad_Error e) {
      logError(e, &params);
      params.setErrorLogLike();
    } catch (SafeVectorOutOfRange) {
      Bad_Error e;
      e.s = "SafeVector out of Range.";
      logError(e, &params);
      params.setErrorLogLike();
    }

    sendResultToMaster(params);
  }
}
