// vim:ts=2:sw=2:et

#include "mcrunner.h"

#include "mcmodel.h"
#include "mclikelihoodcalculator.h"
#include "mcerror.h"
#include "mcutils.h"

#include "multigaussian.h"
#include "rollingaverage.h"

#include "controlpanel.h"
#include "cmbcalc.h"
#include "analyzethis.h"

#include <sstream>
#include <limits>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <ostream>

#include <iomanip>

using namespace std;

McChain::McChain()
       : size(0), totalPerformed(0), stepsSinceUpdate(0),
         EntireFactor(1.), lastPointLogged(false), id(-1)
{
  const int paramCount = McTaskInfo::paramCount();
  if (paramCount<=0) {
    throw McError("McChain::McChain() - parameter count is 0. Call McTaskInfo::initializeParameters() first.");
  }

  roll = new RollingAverage(500);

  sigma.resize(paramCount);
  for (int i=0; i < paramCount; ++i) {
    std::vector<double> v;
    v.resize(paramCount);
    covMatrix.push_back(v);
  }
}

McChain::~McChain()
{
  delete roll;
  roll = 0;
}

static const int MAGIC_CHECK=12345678;
void McChain::readFromDisk(const string& filename)
{
  ifstream readcov;
  readcov.open(filename.c_str());
  if (!readcov) {
    throw McError("Covariance file corrupt or non-existent. Cannot restart.");
  }

  McTaskInfo t = McTaskInfo::fromBinaryStream(readcov);
  cout << prettyPrint << " read in last point in chain: " << t << endl;

  points.push_back(t);
  ++size;

  const unsigned int paramCount = McTaskInfo::paramCount();
  for(unsigned int j=0; j<paramCount; ++j) {
    for (unsigned int m=0; m<paramCount; ++m) {
      covMatrix[j][m] = AnalyzeThis::read<double>(readcov);
    }
  }

  int check = AnalyzeThis::read<int>(readcov);
  if (check != MAGIC_CHECK) {
    throw McError("MAGIC CHECK FAILED");
  }

#warning fix restart when using average entirefactor
  double OptimalFactor = 2.4 / sqrt((double)paramCount);
  EntireFactor = OptimalFactor;
}

const int McRunner::REQUESTTASK=1;
const int McRunner::TAKERESULT =2;
const int McRunner::TAKETASK   =3;
const int McRunner::FINISH   =4;

McRunner::McRunner()
          : mRank(-1)
{
}

McRunner::~McRunner()
{
}

void McRunner::globalInit()
{
  McModel::initialize();
  McLikelihoodCalculator::initialize();
  initBesselFunctionFile();
}

void McRunner::start(bool restart)
{
  cout << "starting rank: " << mRank << endl;
  if (mRank<0) {
    throw McError("rank not set or < 0");
    return;
  }
  init(restart);
  run();
}

void McRunner::init(bool restart)
{
}

McSettings& McRunner::cfg() const
{
  return *McSettings::self();
}

int McRunner::taskInfoSize() const
{
  return McTaskInfo::entryCount();
}

McTaskInfo McRunner::parametersFromMaster(bool* quitRequested)
{
  MPI::COMM_WORLD.Send(0, 0, MPI::DOUBLE, 0, REQUESTTASK);

  MPI::Status nextMsg;
  MPI::COMM_WORLD.Probe(0, MPI::ANY_TAG, nextMsg);

  // check whether master wants us to quit
  if (nextMsg.Get_tag() == FINISH) {
    MPI::COMM_WORLD.Recv(0, 0, MPI::DOUBLE, 0, FINISH);
    if (quitRequested) { // i.e. quitRequested pointer is not null, so the caller wants this info
      *quitRequested=true;
    }
    return McTaskInfo();
  }

  int size = taskInfoSize();
  double* array = new double[size];

  MPI::COMM_WORLD.Recv(array, size, MPI::DOUBLE, 0, TAKETASK);

  McTaskInfo params = McTaskInfo::fromArray(array, size);
  delete[] array;
  return params;
}

void McRunner::sendResultToMaster(const McTaskInfo& result)
{
  double *msg = result.array();
  MPI::COMM_WORLD.Send(msg, result.entryCount(), MPI::DOUBLE, 0, TAKERESULT);
  delete[] msg;
}

void McRunner::askSlavesToQuit()
{
  int remainingSlaves=numberOfChains();
  while (remainingSlaves>0) {
    MPI::Status nextMsg;
    nextMessageFromSlaves(nextMsg);
    int slaveNo = nextMsg.Get_source();
    if (nextMsg.Get_tag() == TAKERESULT) {
      McTaskInfo discard = resultFromSlave(slaveNo);
    } else if ( nextMsg.Get_tag() == REQUESTTASK) {
      sendFinishRequestToSlave(slaveNo);
      --remainingSlaves;
    }
  }
}


void McRunner::nextMessageFromSlaves(MPI::Status& status)
{
  MPI::COMM_WORLD.Probe(MPI::ANY_SOURCE, MPI::ANY_TAG, status);
}

McTaskInfo McRunner::resultFromSlave(const int slaveNo)
{
  int size = taskInfoSize();
  double* array = new double[size];

  MPI::Status status;
  MPI::COMM_WORLD.Recv(array, size, MPI::DOUBLE, slaveNo, TAKERESULT, status);

  McTaskInfo result = McTaskInfo::fromArray(array, size);
  delete[] array;
  return result;
}


void McRunner::sendParametersToSlave(const McTaskInfo& params, const int slaveNo)
{
  double *msg = params.array();
  MPI::COMM_WORLD.Recv(0, 0, MPI::DOUBLE, slaveNo, REQUESTTASK); // consume the request message
  MPI::COMM_WORLD.Send(msg, params.entryCount(), MPI::DOUBLE, slaveNo, TAKETASK);
  delete[] msg;
}

void McRunner::sendFinishRequestToSlave(const int slaveNo)
{
  MPI::COMM_WORLD.Recv(0, 0, MPI::DOUBLE, slaveNo, REQUESTTASK); // consume the request message
  MPI::COMM_WORLD.Send(0, 0, MPI::DOUBLE, slaveNo, FINISH);
}

void McRunner::runSlave()
{
  McModel& model = cfg().model();
  McLikelihoodCalculator calcLike;
  bool quitRequested=false;

  while (true) {
    McTaskInfo params = parametersFromMaster(&quitRequested);
    if (quitRequested) {
      return;
    }
    model.setParameters(params);
    try {
      model.compute();
      calcLike.computeLogLike(model, params);
    } catch (Bad_Error e) {
      logError(e, &params);
      params.setErrorLogLike();
    } catch (bad_alloc oomError) {
      Bad_Error e;
      e.s = "Out of memory: " + string(oomError.what());
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

void McRunner::logError(const Bad_Error& e, McTaskInfo* params)
{
  cout << "\n\n******** BAD ERROR OCCURED IN SLAVE " << mRank << "*****\n\n";
  cout << e.s << "\n\n";
  cout << "******************************************\n";
  ofstream errorFile;
  errorFile.open(cfg().errorFileName().c_str(), ios::app);
  errorFile << "\n\n******** BAD ERROR OCCURED IN SLAVE Nr: " << mRank << "  *****\n\n";
  errorFile<< e.s << "\n\n";
  if (params) {
    errorFile << "Parameters have been: " << endl;
    errorFile << prettyPrint << *params << endl;
    errorFile << "******************************************\n";
    // also save the parameter values only for this error, suitable for use w/ -singlemodel
    ofstream errorParamsFile;
    string name="parameters-"+cfg().errorFileName();
    errorParamsFile.open(name.c_str(), ios::app);
    errorParamsFile << noPrettyPrint << *params << endl;
  }
}

void McRunner::initBesselFunctionFile()
{
  // if there is no jlgen.dat in the resource directory, this will automatically generate
  // bessel functions up to l = 5000
  CmbCalc *cmbcalc = new CmbCalc();
  string cmbeasydir = ControlPanel::cmbeasyDir();
  cmbcalc->initjl(cmbeasydir + "/resources/jlgen.dat", 5000);
  delete cmbcalc;
}

void McStandardSampler::runMaster()
{
  McTaskInfo params;
  MPI::Status nextMsg;
  while (true) {
    nextMessageFromSlaves(nextMsg);
    int slaveNo = nextMsg.Get_source();
    int chainNo = slaveNo-1;
    if (nextMsg.Get_tag() == TAKERESULT) {
      params = resultFromSlave(nextMsg.Get_source());
      updateChain(*mChain[chainNo], params);
    } else if ( nextMsg.Get_tag() == REQUESTTASK) {
      computeNextStepForChain(*mChain[chainNo]);
      sendParametersToSlave(mChain[chainNo]->nextPoint, slaveNo);
    }
  }
}

void McStandardSampler::init(bool restart)
{
  McRunner::init(restart);
  initMaster(restart);
}

void McStandardSampler::initMaster(bool restart)
{
  bool appendToLogFiles = false;
  if (restart) {
    appendToLogFiles = true;
  }
  setupLogFiles(appendToLogFiles);

  McTaskInfo::writeParameterNamesFile();

  setupChains(restart);
  if (restart) {
    restartMaster();
  }
}


void McStandardSampler::setupChains(bool restart)
{
  const int parameterCount = McTaskInfo::paramCount();
  const vector<double> lowerBounds = McTaskInfo::lowerParameterBounds();
  const vector<double> upperBounds = McTaskInfo::upperParameterBounds();

  vector<vector<double> > covMatrix;
  if (cfg().estimateCovariance()) {
      McTaskInfo t;
      for(unsigned int j=0; j<parameterCount; j++) {
        t.setParameterValue(j, lowerBounds[j] + Miscmath::posRnd(upperBounds[j]-lowerBounds[j]));
      }
      t = McUtils::findBestLogLike(t);
      cout << "best loglike found at: " << prettyPrint << t << endl;
      covMatrix = McUtils::estimateCovarianceAroundPoint(t, -t.totalLogLike());
  }
  // this sets the starting points and stepsizes for each chain for the first step
  for (unsigned int i = 0; i < numberOfChains(); i++) {

    McChain *c = new McChain();
    c->id = i;
    mChain.push_back(c);

    MultiGaussian *mg = new MultiGaussian(parameterCount);

    mg->setBounds(lowerBounds, upperBounds);
    mGauss.push_back(mg);

    if (!restart) {
      if (mStartingPoints.empty()) {
        //Generate Starting points (inside the prior)
        McTaskInfo t;
        for(unsigned int j=0; j<parameterCount; j++) {
          t.setParameterValue(j, lowerBounds[j] + Miscmath::posRnd(upperBounds[j]-lowerBounds[j]));
        }
        //cout << "chain " << c->id << ": set starting parameter values: " << prettyPrint << t << endl;

        t.setEntryValue("TotalLogLike", -1e100);
        mChain[i]->points.push_back(t);
      } else {
        std::vector<double> pV = mStartingPoints.back().parameterVector();
        for(unsigned int j=0; j<parameterCount; j++) {
          if (isnan(pV[j])) {
            throw McError("McStandardSampler::setupChains() - trying to set the given starting point: forgot to give "
                           + McTaskInfo::parameterName(j) + " an initial value? ");
          }
        }
        mChain[i]->points.push_back(mStartingPoints.back());
        if (mStartingPoints.size()>1)
          mStartingPoints.pop_back();
      }

      mChain[i]->size = 1;
      mChain[i]->totalPerformed=0;
      logProgress(*mChain[i]);

      mChain[i]->covMatrix.resize(parameterCount);
      for(unsigned int k=0; k< parameterCount;k++) {
        mChain[i]->covMatrix[k].resize(parameterCount);
      }

      if (cfg().estimateCovariance()) {
        mChain[i]->covMatrix = covMatrix;
      } else {
        std::vector<double> initialSigma = McTaskInfo::initialParameterSigmas();
        for(unsigned int k=0; k<parameterCount; k++) {
          for(unsigned int j=0; j<parameterCount; j++) {
            if (j==k)
              mChain[i]->covMatrix[j][j]=initialSigma[j]*initialSigma[j];
            else
              mChain[i]->covMatrix[j][k]=0.0;
          }
        }
      }

      mGauss[i]->generateEigenvectors(mChain[i]->covMatrix, 1.0);
      //mGauss[i]->printExtInfo(cout);
    }
  }
}

void McStandardSampler::restartMaster()
{
#warning fix restart when not frozen
  for (unsigned int k=0; k < numberOfChains(); ++k) {
    McChain& chain = *mChain[k];
    chain.readFromDisk(mCovarianceFileNames[k]);
    mGauss[k]->generateEigenvectors(chain.covMatrix, chain.EntireFactor*chain.EntireFactor);
    mGauss[k]->lock();
    cfg().setAdaptiveStepSize(false);
    cfg().setFreezeIn(false);
  }
  mConverged=true;
}

McStandardSampler::McStandardSampler()
          : McMaster(), mLogDir(""), mProgressLog(0), mGelmanRubinLog(0), mCovarianceLog(0),
            mConverged(false)
{
}

McStandardSampler::~McStandardSampler()
{
  cleanupLogFiles();

  vector<MultiGaussian*>::iterator it, end;
  it = mGauss.begin();
  end = mGauss.end();
  while (it != end) {
    delete *it;
    ++it;
  }

  vector<McChain*>::iterator cit, cend = mChain.end();
  cit = mChain.begin();
  while (cit != cend) {
    delete *cit;
    ++cit;
  }
}

void McStandardSampler::setupLogFiles(bool append)
{
  mLogDir = cfg().logDir();

  const string montecarloName = mLogDir+"montecarlo_chain";
  const string headName = mLogDir+"head";
  const string investigatedName = mLogDir+"investigated";
  const string covName = mLogDir+"covfile";

  ios::openmode Mode = ios::out;
  if (append) {
    Mode = ios::app;
  }


  for (unsigned int k=0; k<numberOfChains(); ++k) {
    stringstream numStr;
    numStr << k;
    string chainName = montecarloName+numStr.str()+".dat";
    string headFileName = headName+numStr.str()+".txt";
    string investigatedFileName = investigatedName+numStr.str()+".dat";
    string covFileName = covName+numStr.str()+".dat";

    mChainLog.push_back(new ofstream(chainName.c_str(), Mode));
    mHeadLog.push_back(new ofstream(headFileName.c_str(), Mode));
    mInvestigatedLog.push_back(new ofstream(investigatedFileName.c_str(), Mode));
    mCovarianceFileNames.push_back(covFileName.c_str());


    mChainLog[k]->setf(ios::scientific);
    mHeadLog[k]->setf(ios::scientific);
    mInvestigatedLog[k]->setf(ios::scientific);
    (*mChainLog[k]) << convertNanToZero;
    (*mHeadLog[k]) << convertNanToZero;
    (*mInvestigatedLog[k]) << convertNanToZero;
  }

  mProgressLog = new ofstream((mLogDir+"progress.txt").c_str(), Mode);
  if (append) {    // so we are restarting the chain
    (*mProgressLog) << "\n\n::::::::::::::::::::::::::::::::::::::::::";
    (*mProgressLog) << "::::::::::::::::: RESTARTED FROM HERE ::::::::";
    (*mProgressLog) << ":::::::::::::::::::::::::::::\n\n";
  }

  mGelmanRubinLog = new ofstream((mLogDir+"gelmanRubin.txt").c_str(), Mode);
  mCovarianceLog = new ofstream((mLogDir+"covMatrix.txt").c_str(), Mode);
  mCovarianceLog->setf(ios::scientific);
}

void McStandardSampler::cleanupLogFiles()
{
  for (unsigned int k=0; k<mChainLog.size(); ++k) {
    delete mChainLog[k];
    delete mHeadLog[k];
    delete mInvestigatedLog[k];
  }

  delete mProgressLog;
  delete mGelmanRubinLog;
  delete mCovarianceLog;
}

int McRunner::numberOfChains() const
{
  return cfg().numberOfChains();
}

static bool throwDice(const McTaskInfo& oldPoint, McTaskInfo& newPoint, MultiGaussian *gauss)
{
  bool success = gauss->throwDice(oldPoint.parameterVector());
  newPoint.setParameterValues(gauss->getRandomValues());
  return success;
}

void McStandardSampler::computeNextStepForChain(McChain& chain)
{
  const int chainNo = chain.id;
  McTaskInfo& lastPoint = chain.points.back();

  //cout << prettyPrint << lastPoint << endl;

  while (!throwDice(lastPoint, chain.nextPoint, mGauss[chainNo])) {
    // as long as dice roll gives parameters
    // outside the parameter bounds
    // increase the weight of current point
    lastPoint.Multiplicity++;
  }
  //cout << "new point: " << prettyPrint << chain.nextPoint << endl;
  vector<double> pV = chain.nextPoint.parameterVector();
  for (int i = 0; i < pV.size(); ++i) {
    if (isnan(pV[i]) || isinf(pV[i])) {
        stringstream s;
        s << McTaskInfo::parameterName(i) << " has new value nan or inf: " << pV[i];
        throw McError("McStandardSampler::computeNextStepForChain() -  " + s.str());
    }
  }
}


/*! Dynamically adjust the stepsize
 * In addition and to help the adaptive covariance,
 * we use an adaptive step size multiplicator, called EntireFactor
 * Whenever we take a step, we consider increasing the step size,
 * as frequent approval corresponds to small stepsizes.
 * We also consider the opposite: If we hang around
 * for some time at the same spot, we decrease the step size..
 *
 * \param tookStep true if called after we just accepted a new point
 */
void McStandardSampler::adjustAdaptiveStepSize(McChain& chain, bool tookStep)
{
  if(!cfg().adaptiveStepSize()) {
    return;
  }

  chain.roll->push(chain.EntireFactor);

  // decrease if we didn't take for a long time
  double stepCount;
  if (cfg().useReallyInvestigated()) {
    stepCount = chain.points.back().ReallyInvestigated;
  } else {
    stepCount = chain.points.back().Multiplicity;
  }
  if (stepCount > cfg().highStepBound() /*&& chain.size>=cfg().beginCovUpdate()*/) {
    if (chain.EntireFactor > 1e-1) {
      chain.EntireFactor *= cfg().stepDecrease();
      mGauss[chain.id]->generateEigenvectors(chain.covMatrix, chain.EntireFactor*chain.EntireFactor);
      //cout << chain.id << "[" << chain.points.back().ReallyInvestigated << "]"
      //     <<" decreased to " << chain.EntireFactor << endl;
    }
  }

  // the rest of the adaptive stepsize algorithm is only run once per new point in the chain
  if (!tookStep) {
    return;
  }

  // increase stepsize, if we take steps  too often
  list<McTaskInfo>::iterator prev = chain.points.end();
  prev--; prev--; // prev points to the point before the new one
  McTaskInfo& previousPoint = *prev;
  if (previousPoint.ReallyInvestigated < cfg().lowStepBound() /*&& chain.size>=cfg().beginCovUpdate()*/) {
    if (chain.EntireFactor < 10) {
      chain.EntireFactor *= cfg().stepIncrease();
      mGauss[chain.id]->generateEigenvectors(chain.covMatrix, chain.EntireFactor*chain.EntireFactor);
      //cout << chain.id << "[" << previousPoint.ReallyInvestigated << "]"
      //     <<" increased to " << chain.EntireFactor << endl;
    }
  }

  /* Adaptivly setting stepsize according to covariance matrix
   * Here, we estimate the covariance
   */

  if(chain.stepsSinceUpdate>=cfg().updateTime() && chain.size>=cfg().beginCovUpdate()) {
    const int paramCount = McTaskInfo::paramCount();

    // Here, we determine the number of distinct points in the chain
    // that will be taken into account to estimate the covariance matrix
    unsigned int covSize;
    if (chain.size < 2*cfg().beginCovUpdate()) { // early on
      covSize = chain.size/2;
    } else {
      covSize = chain.size - cfg().beginCovUpdate();  // all but the first few begincovupdate points
    }
    covSize = min(covSize, cfg().maxPreviousCovPoints());

    list<McTaskInfo>::iterator chainIter, startPointIter;
    chainIter = chain.points.begin();

    // ignore the first part of the chain
    for (unsigned int j=0; j<chain.size-covSize; j++) {
      ++chainIter;
    }

    startPointIter=chainIter;

    McUtils::updateCovarianceMatrix(startPointIter, chain.points.end(), chain.covMatrix);

    mGauss[chain.id]->generateEigenvectors(chain.covMatrix, chain.EntireFactor*chain.EntireFactor);
    //(*mProgressLog) << "Extended Information for chain: " << chain.id << ": ";
    //mGauss[chain.id]->printExtInfo(*mProgressLog);


    // output information, in text file as well as binary format for re-starting
    (*mCovarianceLog) << "Chain: " << chain.id << " Step: "
                      << chain.size << " Points used: " << covSize << endl;

    ofstream covfile(mCovarianceFileNames[chain.id].c_str());

    double *taskInfo = chain.points.back().array();
    for (unsigned int j=0; j<McTaskInfo::entryCount(); ++j) {
      AnalyzeThis::write<double>(covfile, taskInfo[j]);
    }

    AnalyzeThis::write<int>(covfile, chain.points.back().Multiplicity);
    AnalyzeThis::write<int>(covfile, chain.points.back().ReallyInvestigated);

    for(unsigned int j=0; j<paramCount; ++j) {
      for (unsigned int k=0; k<paramCount; ++k) {
        AnalyzeThis::write<double>(covfile, chain.covMatrix[j][k]);
        (*mCovarianceLog) << chain.covMatrix[j][k] << "  ";
      }
      (*mCovarianceLog) << endl;
    }
    AnalyzeThis::write<int>(covfile, MAGIC_CHECK);

    (*mCovarianceLog) << "Entire factor: " << chain.EntireFactor << endl;
    (*mCovarianceLog) << "***********************************************************" << endl;

    chain.stepsSinceUpdate=0;
  }
}

void McStandardSampler::logProgress(const McChain& chain)
{
  (*mHeadLog[chain.id]) << chain.nextPoint << endl;

  if (chain.lastPointLogged || chain.size == 0) {
    return;
  }

  ofstream& log = *mChainLog[chain.id];
  log << chain.points.back() << endl;
  chain.lastPointLogged = true;

  vector<McChain*>::const_iterator it, end;
  (*mProgressLog) << "performed/Chainsize: ";

  end = mChain.end();

  unsigned int minSize=chain.size;
  for( it = mChain.begin(); it !=end; ++it) {
    McChain *c = *it;
    (*mProgressLog) << c->totalPerformed << " / " << c->size << "  ";
    minSize = min(minSize, c->size);
  }
  (*mProgressLog) << "min_size: " << minSize   << "  Multiplicity [Really]: ";

  for( it = mChain.begin(); it !=end; ++it) {
    McChain *c = *it;
    (*mProgressLog) << c->points.back().Multiplicity << " [" << c->points.back().ReallyInvestigated << "]  ";
  }
  (*mProgressLog) << endl;
}

void McStandardSampler::updateChain(McChain& chain, McTaskInfo& result)
{
  chain.totalPerformed++;

  bool take = false;
  double loglike1 = chain.points.back().totalLogLike();
  double resultLogLike = result.totalLogLike();
  bool   noError=true;
  if (isnan(resultLogLike)) { // error occured
    resultLogLike = -1e100;
    noError=false;
  }

  if (resultLogLike<=loglike1 && noError) {
    // calculate likeli2 / likeli1. We subtract loglike2 first, doesn't change
    // result, but helps prevent 0/0
    double ratio = 1.0/exp( loglike1 - resultLogLike);
    double x = Miscmath::posRnd(1.0);
    if (x < ratio) {
      take=true;
    }
  } else if (noError) {
    take = true; // always take if the likelihood is larger than the previous
  }

  if (take) {
    logProgress(chain);
    chain.points.push_back(result);
    chain.lastPointLogged = false;
    chain.size++;
    chain.stepsSinceUpdate++;  // a taken step contributes towards stepsSinceUpdate counting
    stepTaken(chain, result);
  } else {
    chain.points.back().Multiplicity++;  // enhance multiplicty of the old one again
    chain.points.back().ReallyInvestigated++; // we did a simulation and we didn't take the step
    stepRejected(chain, result);
  }
}

void McStandardSampler::stepTaken(McChain& chain, McTaskInfo& step)
{
  (*mInvestigatedLog[chain.id]) << step << endl;
  updateGelmanRubinStatistics();
  adjustAdaptiveStepSize(chain, true /* just accepted a new point */);
}

void McStandardSampler::stepRejected(McChain& chain, McTaskInfo& step)
{
  (*mInvestigatedLog[chain.id]) << step << endl;
  adjustAdaptiveStepSize(chain, false /* rejected last proposal */);
}

/*!Gelman and Rubin(1992) statistic
 * compare variances within the chain with the variances between the chains
 * R[k] should be smaller than 1.1 for convergence
 * If you would like to run a synthetic distribution for checks, this statistics
 * will slow you down. You can speed things up if you replace the if
 * statement below by some other test, or even by something like:
 *        if ( (min_size) > RMIN_POINTS && Miscmath::posRnd(1.0) > 0.95) {
 * hence, only about every 20 times, the statistics is re-calculated.
 */
void McStandardSampler::updateGelmanRubinStatistics()
{
  unsigned int minSize=mChain[0]->size;
  for (unsigned int i=0; i<numberOfChains(); ++i) {
    minSize = min(minSize, mChain[i]->size);
  }

  static unsigned int minSizeAtLastUpdate = 0;
  if (minSize-minSizeAtLastUpdate < cfg().statisticsUpdateSteps()) {
    return;
  }

  minSizeAtLastUpdate = minSize;

  const unsigned int paramCount = McTaskInfo::paramCount();
  const unsigned int chainCount = numberOfChains();

  if (minSize < cfg().minRStatPoints()) {
    return;
  }
  unsigned int breakAt = (minSize)/2;
  double N = minSize - breakAt;
  vector<list<McTaskInfo>::iterator> iter(numberOfChains());
  for (unsigned int i = 0; i<numberOfChains(); ++i) {
    iter[i] = mChain[i]->points.begin();
    // ignore the first half of the minimal chain size
    for (unsigned int k =1; k<breakAt; k++, iter[i]++)
      ;
  }


  // now let us calculate the mean of each parameter
  double y[chainCount][paramCount];
  double dist_y[paramCount]; // distribution mean
  double B[paramCount]; // variance between chains
  double W[paramCount]; // variance within chains
  double R[paramCount]; // monitoring parameter
  for (unsigned int k =  0; k<paramCount; ++k) {
    dist_y[k]=0;
    B[k]=0;
    W[k]=0;
    for (unsigned int i=0; i < chainCount; ++i)
      y[i][k]=0;
  }
  double M = chainCount;
  for (unsigned int i=0; i<chainCount; ++i) {
    int TotalMultiplicity=0;
    for (unsigned int n=breakAt; n<minSize; ++n) {  // traverse through the chains
      for (unsigned int k=0; k<paramCount; ++k) { // for all parameters
        // add up all parameters
        vector<double> params = iter[i]->parameterVector();
        y[i][k] += params[k]*iter[i]->Multiplicity;
      }
      TotalMultiplicity += iter[i]->Multiplicity;
      iter[i]++;
    }
    // (*mProgressLog) << "new eval: N: " << N <<  "  check_count: " << check_count << endl;
    N = (double) TotalMultiplicity;
    for (unsigned int k =0; k<paramCount; ++k){
      y[i][k] /= (double) N;
      dist_y[k] += y[i][k]/M;  // mean is sum / chainCount
    }
  }

  // variance between chains
  for (unsigned int k=0; k<paramCount; ++k) {
    for (unsigned int i=0; i<chainCount; ++i) {
      double t1 = y[i][k]-dist_y[k];
      B[k] += t1*t1/(M-1.0);  // as defined in (22) (why 1/(M-1) and not M ???)
    }
  }
  // once again, go back to break_at in the chains
  for (unsigned int i=0; i<chainCount; ++i) {
    iter[i] = mChain[i]->points.begin();
    for (unsigned int k=1; k<breakAt; ++k, iter[i]++)
      ;
  }

  // get W, the variance within chains
  for (unsigned int i=0; i<chainCount; ++i) {  // run over all chains
    for (unsigned int n=breakAt; n<minSize; ++n) { // run over all points
      for (unsigned int k=0; k<paramCount; ++k) {
        vector<double> params = iter[i]->parameterVector();
        double t2 = params[k]-y[i][k];
        W[k] += t2*t2*iter[i]->Multiplicity;
      }
      iter[i]++; //next point
    }
  }
  // normalize W and get R
  for (unsigned int k =0; k<paramCount; ++k) {
    W[k] /= (M*(N-1.0));
    R[k] = (N-1.0)/N*W[k]+B[k]*(1.0+1.0/N);
    R[k] /= W[k];
  }

  vector<double> mean(dist_y, dist_y+paramCount);
  vector<double> vB(B, B+paramCount);
  vector<double> vW(W, W+paramCount);
  vector<double> vR(R, R+paramCount);

  bool toProgressLogOnly = true;
  static unsigned int lastLoggedChainSize=0;
  if (lastLoggedChainSize>minSize) { // after resetting the chains (e.g. after they've converged)
    lastLoggedChainSize=0;
  }
  if (minSize>lastLoggedChainSize) {
    lastLoggedChainSize=minSize;
    toProgressLogOnly=false;
  }
  logGelmanRubinStatistics(mean, vB, vW, vR, toProgressLogOnly, minSize, breakAt);
  checkConvergence(minSize, vR);
}

void McStandardSampler::logGelmanRubinStatistics(const vector<double>& mean, const vector<double>& B,
                                         const vector<double>& W, const vector<double> R,
                                         bool toProgressLogOnly, unsigned int minSize,
                                         unsigned int startPoint)
{
  const unsigned int paramCount = McTaskInfo::paramCount();
  const unsigned int chainCount = numberOfChains();

  (*mProgressLog) << "Statistics: " << endl;

  for (unsigned int k=0; k<paramCount; ++k) {
    (*mProgressLog) << k << " - " << McTaskInfo::parameterName(k) << ": "
                    << mean[k] << "  B: " << B[k]
                    << "  W: " << W[k] << "    R[k]:  " << R[k] << endl;
  }

  (*mProgressLog) <<  "Multiplicity (and " << (!cfg().freezeIn()?"frozen":"") << " EntireFactor): ";
  for (unsigned int i = 0; i<chainCount; i++) {
   (*mProgressLog) << mChain[i]->points.back().Multiplicity << " ("<< mChain[i]->EntireFactor << ")  ";
  }
  (*mProgressLog) << endl;

  if (!toProgressLogOnly) {
    (*mGelmanRubinLog) << startPoint << "  ";
    for (unsigned int n=0 ; n<paramCount; ++n) {
      (*mGelmanRubinLog) << R[n] << "  ";
    }
    (*mGelmanRubinLog) << minSize << endl;
  }
}

void McStandardSampler::checkConvergence(unsigned int minChainSize, const vector<double>& R)
{
  if(!cfg().freezeIn() || minChainSize<cfg().minFreezeInPoints()) {
    return;
  }

  const unsigned int paramCount = McTaskInfo::paramCount();
  const unsigned int chainCount = numberOfChains();

  for(unsigned int m=0; m<paramCount; ++m) {
    if(R[m]>cfg().freezeInR())
      return; // not converged yet
  }

  mConverged = true;

  // ok, let's freeze stepsize, covariance matrices and Entire factor:

  cfg().setAdaptiveStepSize(false); //stop with adaptive stepsize
  cfg().setFreezeIn(false); //never check again (otherwise no real freeze-in)

  ofstream final((mLogDir+"freezeInEigenvectors.txt").c_str());

  for (unsigned int j=0; j<chainCount; j++) {
    McTaskInfo keep=mChain[j]->points.back();

    keep.Multiplicity = 0; // zero Multiplicity signals freeze-in
    double averageEntireFactor = mChain[j]->roll->average();
    std::vector<std::vector<double> > covMatrix = mChain[j]->covMatrix;

    // all models before freeze-in have to be discarded
    delete mChain[j];
    mChain[j] = new McChain();
    McChain& chain = *mChain[j];
    chain.id = j;

    // and we push back the last one before freeze in...
    chain.points.push_back(keep);
    chain.size=1;
    chain.covMatrix=covMatrix;
    logProgress(chain);

    double optimalFactor;
    if (cfg().useAverageEntireFactor()) {
      optimalFactor = averageEntireFactor;
    } else {       // use Dunkley et. al. result for optimal sigma_t: 2.4 / sqrt(Parameters)
      optimalFactor = 2.4 / sqrt((double)paramCount);
    }

    chain.EntireFactor=optimalFactor;
    mGauss[j]->generateEigenvectors(chain.covMatrix, optimalFactor*optimalFactor);
    mGauss[j]->lock();

    final << "::: CHAIN[" << j << "] Eingenvector and value info: " << endl;
    mGauss[j]->printInfo(final);
    final << endl << endl;
  }

  (*mCovarianceLog) << "Stopped adaptive stepsize at "
                    << minChainSize << " for all chains (FREEZE_IN=true)" << endl;
  mCovarianceLog->close();
}


