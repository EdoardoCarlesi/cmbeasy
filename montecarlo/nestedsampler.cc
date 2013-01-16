// vim:ts=2:sw=2:et

#include "nestedsampler.h"

#include "mclikelihoodcalculator.h"
#include "mcutils.h"

#include "multigaussian.h"
//#include "kmlocal-1.7.1/src/KMlocal.h" // k-means algorithms

#include <algorithm>
#include <numeric>
#include <limits>
#include <iostream>
#include <iomanip>
#include <iterator>

using namespace std;

NestedSampler::NestedSampler()
             : McMaster(), mLivePointsCount(100),
               mPointsComputed(0), mPointsUsed(0),
               mEllipsoidEnlargementFactor(1.8),
               mSamplingMode(UniformOverPrior),
               mMinIterations(100),
               mMaxRemainingLogEvidenceFraction(0.1)
{
}

NestedSampler::~NestedSampler()
{
}


void NestedSampler::init(bool restart)
{
  McMaster::init(restart);

  mWasRestarted=restart;
  if (restart) {
    readLivePoints();
    return;
  }

  McTaskInfo::writeParameterNamesFile();

  mInitialPoints.resize(mLivePointsCount);
  generate(mInitialPoints.begin(), mInitialPoints.end(), McUtils::randomParameterPoint);

  mCurrentIteration=0;
  mErrorCount=0;
  mCurrentLogX_i=mCurrentLogWeight=mInformation=0;
  mLogEvidence=mBestLiveLnLike=-numeric_limits<double>::infinity();
}

void NestedSampler::readLivePoints()
{
  ifstream file1("livepointsfile.dat");
  ifstream file2("livepointsfile1.dat");
  unsigned long int c1=count( istream_iterator<char>(file1), istream_iterator <char> (), '\n');
  unsigned long int c2=count( istream_iterator<char>(file2), istream_iterator <char> (), '\n');
  file1.close();
  file2.close();
  ifstream livePointsFile;
  if (c1>c2) {
    livePointsFile.open("livepoints.dat");
  } else {
    livePointsFile.open("livepoints.dat");
  }

  if (!livePointsFile) {
    throw McError("NestedSampler::init() - restarting, but can't open file 'livepoints.dat'");
  }
  livePointsFile >> mLivePointsCount >> mPointsUsed
    >> mErrorCount >> mPointsComputed
    >> mCurrentLogX_i >> mCurrentLogWeight
    >> mInformation >> mLogEvidence
    >> mBestLiveLnLike >> mCurrentIteration;
  McTaskInfo readPoint;
  livePointsFile >> readPoint;
  while (livePointsFile) {
    mLivePoints.push_back(readPoint);
    livePointsFile >> readPoint;
  }
  if (mLivePoints.size() == mLivePointsCount) {
    cout << "restarting with " << mLivePointsCount << " live points." << endl;
  } else {
    cout << "restarting with " << mLivePoints.size() << " live points, expected " << mLivePointsCount
      << "; continuing anyway." << endl;
  }
}

void NestedSampler::runMaster()
{
  //testSampling();
  McTaskInfo params;
  MPI::Status nextMsg;
  while (!terminateCondition()) {
    nextMessageFromSlaves(nextMsg);
    int slaveNo = nextMsg.Get_source();
    if (nextMsg.Get_tag() == TAKERESULT) {
      params = resultFromSlave(nextMsg.Get_source());
      updatePoints(params);
    } else if ( nextMsg.Get_tag() == REQUESTTASK) {
      params = nextParameterSet();
      sendParametersToSlave(params, slaveNo);
    }
  }
  completeEvidence();
  writePosteriorSamples();
  cout << "Sampling efficiency (points used/points computed): "
       << (double(mPointsUsed)/double(mPointsComputed))
       << " (" << mPointsUsed << "/" << mPointsComputed
       << ")"  << endl;
  cout << "ln Z = " << mLogEvidence << ", z = " << exp(mLogEvidence) << endl;
}

McTaskInfo NestedSampler::nextParameterSet()
{
  McTaskInfo newSet;

  if (mInitialPoints.empty() && (mLivePointsCount>mLivePoints.size())) {
    mInitialPoints.resize(mLivePointsCount-mLivePoints.size());
    generate(mInitialPoints.begin(), mInitialPoints.end(), McUtils::randomParameterPoint);
  }

  if(!mInitialPoints.empty()) {
    newSet.setParameterValues(mInitialPoints.back());
    mInitialPoints.pop_back();
    return newSet;
  }

  switch (mSamplingMode) {
    case UniformOverPrior:
      newSet.setParameterValues(McUtils::randomParameterPoint());
      break;
    case EllipsoidalSampling:
      newSet = nextParameterSetEllipsoidal();
      break;
    //case MultiNest:
    //  newSet = nextParameterSetMultiNest();
    //  break;
    default:
      throw McError("NestedSampler::nextParameterSet() - "
                    "unknown sampling method");
  }

  return newSet;
}


McTaskInfo NestedSampler::nextParameterSetEllipsoidal()
{
  using namespace McUtils;

  static unsigned int paramCount=McTaskInfo::paramCount();
  vector<double> mean(paramCount);
  static vector<vector<double> > mLiveCovarianceMatrix(paramCount, vector<double>(paramCount));

  static const int updateInterval=1;
  static int updateCounter=updateInterval;
  if (++updateCounter>=updateInterval) {
    updateCounter=0;
    updateCovarianceMatrix(mLivePoints.begin(), mLivePoints.end(),
                           mLiveCovarianceMatrix, &mean);
  }

  McTaskInfo newSet;

  Matrix eVecs = eigenVectors(mLiveCovarianceMatrix);
  Matrix eVecsT = transposed(eVecs);

  Matrix diag = matrixMatrixProd(mLiveCovarianceMatrix, eVecs);
  diag = matrixMatrixProd(eVecsT, diag);
  /*
     cout << "covMatrix:\n" << covMatrix << endl;
     cout << "evecs:\n" << eVecs << endl;
     cout << "evecsT:\n" << eVecsT << endl;
     cout << "diag:\n" << diag << endl;
  */
  unsigned int dim = paramCount;
  Matrix sqrtDiag(dim, vector<double>(dim));
  for (int i=0; i<diag.size(); ++i) {
    sqrtDiag[i][i] = sqrt(diag[i][i]);
  }
  Matrix invCov = invertMatrix(mLiveCovarianceMatrix);
  double k_max = -std::numeric_limits<double>::infinity();
  vector<McTaskInfo>::iterator pointsIterator=mLivePoints.begin();
  for ( ; pointsIterator!=mLivePoints.end(); ++pointsIterator) {
    vector<double> pV = pointsIterator->parameterVector();
    transform(pV.begin(), pV.end(), mean.begin(), pV.begin(), minus<double>());
    vector<double> v = matrixVectorProd(invCov, pV);
    k_max = max(k_max, inner_product(pV.begin(), pV.end(), v.begin(), 0.));
  }
  Matrix trafo = matrixMatrixProd(eVecs, sqrtDiag);

  vector<double> res, random(dim);
  do {
    double norm=0;
    for (int j=0; j<dim; ++j) {
      random[j]=Miscmath::gaussRnd();
      norm += random[j]*random[j];
    }
    double dist=pow(Miscmath::posRnd(), 1./dim);
    for (int j=0; j<dim; ++j) {
      random[j]*=dist/sqrt(norm);
    }
    res = matrixVectorProd(trafo, random);
    transform(res.begin(), res.end(), res.begin(), bind1st(multiplies<double>(), mEllipsoidEnlargementFactor*sqrt(k_max)));
    transform(res.begin(), res.end(), mean.begin(), res.begin(), plus<double>());
  } while (!McUtils::insideBounds(res));

  newSet.setParameterValues(res);
  return newSet;
}

/*
McTaskInfo NestedSampler::nextParameterSetMultiNest()
{
  throw McError("Multi nest not implemented in this version.");
}
*/

// for inverse sorting w.r.t. likelihood
struct SmallerLikelihood
{
  bool operator()(const McTaskInfo& i1, const McTaskInfo& i2) const {
    return (i1.totalLogLike() > i2.totalLogLike());
  }
};

void NestedSampler::updatePoints(McTaskInfo& result)
{
  double resultLogLike = result.totalLogLike();
  if (isnan(resultLogLike)) { // error occured
    ++mErrorCount;
    return;
  }

  ++mPointsComputed;
  static ofstream investigatedLog("investigated.dat", ios::app);
  investigatedLog << noPrettyPrint << result << endl;

  if (mLivePoints.size()<mLivePointsCount) {
    mLivePoints.push_back(result);
    mBestLiveLnLike=max(mBestLiveLnLike, result.totalLogLike());
    ++mPointsUsed;
    return;
  }

  static bool madeHeapDone=false;
  if (!madeHeapDone) {
    make_heap(mLivePoints.begin(), mLivePoints.end(), SmallerLikelihood());
    madeHeapDone=true;
  }

  double smallestLikeInLivePoints=mLivePoints.front().totalLogLike();
  if (result.totalLogLike()<=smallestLikeInLivePoints) {
    return;
  }

  updateEvidence();

  mLivePoints.push_back(result);
  push_heap(mLivePoints.begin(), mLivePoints.end(), SmallerLikelihood());
  mBestLiveLnLike=max(mBestLiveLnLike, result.totalLogLike());
  ++mPointsUsed;

  static ofstream progressLog("progress.txt", ios::app);
  progressLog << "used/error/computed: " << mPointsUsed << "/" << mErrorCount << "/"
              << mPointsComputed << " || mLogEvidence (max. est. remaining fraction): " << mLogEvidence
              << " (" << (mCurrentLogX_i+mBestLiveLnLike-mLogEvidence) << ")"
              << "   || live LnLike max/min: " << mBestLiveLnLike << "/"
              << mLivePoints.front().totalLogLike()
              << endl;

  static int fileNo=0; // we use two different files to save the live points and state,
                       // so that when we're interrupted while writing to one of them,
                       // we still have the complete one from before
  string liveFileName="livepoints.dat";
  if (++fileNo>1) {
    fileNo=0;
    liveFileName="livepoints1.dat";
  }
  ofstream livePointsFile(liveFileName.c_str());
  std::string sep("    ");
  livePointsFile << scientific << setprecision(10) << noPrettyPrint;
  livePointsFile << mLivePointsCount << sep << mPointsUsed << sep << mErrorCount << sep << mPointsComputed
                 << sep << mCurrentLogX_i << sep << mCurrentLogWeight
                 << sep << mInformation << sep << mLogEvidence
                 << sep << mBestLiveLnLike << sep << mCurrentIteration << "\n";
  copy(mLivePoints.begin(), mLivePoints.end(), ostream_iterator<McTaskInfo>(livePointsFile, "\n"));
  livePointsFile.flush();
}

void NestedSampler::updateEvidence()
{
  double i = ++mCurrentIteration;
  mCurrentLogX_i=-i/mLivePointsCount;
  //log(0.5*(X_(i-1)-X_(i+1)))
  double logWeight=-(i-1.)/mLivePointsCount+log(0.5*(1.-exp(-2./mLivePointsCount)));
  mCurrentLogWeight=logWeight;
  DiscardedPoint p;
  p.point=mLivePoints.front();
  p.logWeight=logWeight;
  mDiscardedPoints.push_back(p);
  pop_heap(mLivePoints.begin(), mLivePoints.end(), SmallerLikelihood());
  mLivePoints.pop_back();

  double discardedPointLnLike=p.point.totalLogLike();
  double newLogEvidence=Miscmath::logAdd(mLogEvidence, logWeight+discardedPointLnLike);
  //mInformation=exp(logWeight-newLogEvidence)*smallestLikeInLivePoints
  //             +exp(mLogEvidence-newLogEvidence)
  //             *(mInformation+mLogEvidence)-newLogEvidence;

  mLogEvidence=newLogEvidence;

  // log discarded points to a file
  double weight=exp(logWeight+discardedPointLnLike);
  static ofstream logfile("discardedpoints.dat", ios::app);

  logfile << noPrettyPrint << p.point
          << "    " << logWeight
          << "    " << mLogEvidence
          << "    " << mInformation
          << endl;
}

void NestedSampler::completeEvidence()
{
  double i = mCurrentIteration+1;
  double logWeight=-i/mLivePointsCount-log(mLivePointsCount);
  mCurrentLogWeight=logWeight;

  while (!mLivePoints.empty()) {
    DiscardedPoint p;
    p.point=mLivePoints.back();
    p.logWeight=logWeight;
    mDiscardedPoints.push_back(p);
    mLivePoints.pop_back();

    double discardedPointLnLike=p.point.totalLogLike();
    mLogEvidence=Miscmath::logAdd(mLogEvidence, logWeight+discardedPointLnLike);

    double weight=exp(logWeight+discardedPointLnLike);

    ofstream logfile("discardedpoints.dat", ios::app);
    logfile << noPrettyPrint << p.point
            << "    " << logWeight
            << "    " << mLogEvidence
            << "    " << mInformation
            << endl;
  }
}

bool NestedSampler::terminateCondition() const
{
  if (mLivePoints.size()<mLivePointsCount) {
    return false;
  }
  if (mDiscardedPoints.size()<mMinIterations) {
    return false;
  }
  return (mCurrentLogX_i+mBestLiveLnLike-mLogEvidence<mMaxRemainingLogEvidenceFraction);
}


void NestedSampler::produceCurrentPosterior()
{
  readLivePoints();
  completeEvidence();
  writePosteriorSamples();
  cout << "Sampling efficiency (points used/points computed): "
       << (double(mPointsUsed)/double(mPointsComputed))
       << " (" << mPointsUsed << "/" << mPointsComputed
       << ")"  << endl;
  cout << "ln Z = " << mLogEvidence << ", z = " << exp(mLogEvidence) << endl;
}


void NestedSampler::writePosteriorSamples()
{
  cout << "writing posterior file." << endl;
  if (mWasRestarted) {
    ifstream discardedPointsFile("discardedpoints.dat");
    double logEv, information;
    unsigned long int count = 0;
    while (discardedPointsFile) {
      DiscardedPoint p;
      discardedPointsFile >> p.point
                          >> p.logWeight
                          >> logEv >> information;
      if (discardedPointsFile) {
        mDiscardedPoints.push_back(p);
        ++count;
      }
    }
    cout << "used " << count << " discarded points from prior runs." << endl;
  }


  ofstream file("posterior.dat");
  for (int i=0; i<=McTaskInfo::entryCount(); ++i) {
    file << setprecision(8)<< -1.0 << "    ";
  }
  file << "    " << 0 << endl;

  unsigned int size=mDiscardedPoints.size();
  double maxLogWeight=-1e100;
  for (unsigned int i=0; i<size; ++i) {
    double logWeight=mDiscardedPoints[i].logWeight+mDiscardedPoints[i].point.totalLogLike()-mLogEvidence;
    maxLogWeight = max(maxLogWeight, logWeight);
  }
  for (unsigned int i=0; i<size; ++i) {
    double logWeight=mDiscardedPoints[i].logWeight+mDiscardedPoints[i].point.totalLogLike()-mLogEvidence;
    unsigned long int multiplicity = (unsigned long int)nearbyint(exp(logWeight)/exp(maxLogWeight)*3e4);
    if (multiplicity>0) {
      file << mDiscardedPoints[i].point << "    " << multiplicity << endl;
    }
  }
  file.close();
  cout << "done." << endl;
}
