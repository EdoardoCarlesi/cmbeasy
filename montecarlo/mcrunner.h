#ifndef MCRUNNER_H
#define MCRUNNER_H
// vim:ts=2:sw=2:et

#include "mcsettings.h"
#include "mctaskinfo.h"

#include "mpi.h"

#include <fstream>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <vector>


class MultiGaussian;
class Bad_Error;

class RollingAverage;

class McChain
{
  public:
    McChain();
    ~McChain();

    void readFromDisk(const std::string& filename);

    std::list<McTaskInfo> points;

    unsigned int size;
    unsigned int totalPerformed;

    // for dynamical stepsize
    int stepsSinceUpdate;
    double EntireFactor;
    RollingAverage* roll;   // rolling average of EntireFactor for after burn-in
    std::vector<double> sigma;
    std::vector<std::vector<double> > covMatrix;
    McTaskInfo nextPoint;

    mutable bool lastPointLogged;

    int id;
};


class McRunner
{
  public:

    McRunner();
    virtual ~McRunner();

    void setRank(int r) { mRank = r; }
    void start(bool restart=false);

    static  void globalInit();
    virtual void init(bool restart=false);

  protected:
    /*! run() calls runSlave() in the default implementation.
     *  A master should inherit from McMaster,
     *  and reimplement McMaster::runMaster();
     */
    virtual void run() { runSlave(); }
    virtual void runSlave();

    McTaskInfo parametersFromMaster(bool* quitRequested=0);
    void       sendResultToMaster(const McTaskInfo& result);

    void       nextMessageFromSlaves(MPI::Status& status);
    McTaskInfo resultFromSlave(const int slaveNo);
    void       sendParametersToSlave(const McTaskInfo& params, const int slaveNo);
    void       askSlavesToQuit();
    void       sendFinishRequestToSlave(const int slaveNo);

    int         taskInfoSize() const;
    int         numberOfChains() const;
    McSettings& cfg() const;

    void fatalError(const std::string&);
    void logError(const Bad_Error& e, McTaskInfo* params=0);

    static void initBesselFunctionFile();

    // Some constants needed for MPI communications:
    static const int REQUESTTASK; //!< MPI flag indicating slave's wish to receive a task from master()
    static const int TAKERESULT;  //!< MPI flag indicating slave's request to master() to take result
    static const int TAKETASK; //!< MPI flag indicating master's request to slave to take task.
    static const int FINISH; //!< MPI flag indicating master's request to slave to quit

  private:
    McRunner(const McRunner&);

  private:
    int mRank;
};

class McMaster: public McRunner
{
  public:
    McMaster(): McRunner() {}
    virtual ~McMaster()  {};

    virtual void run() { runMaster(); askSlavesToQuit(); return; }
    virtual void runMaster()=0;
};

class McStandardSampler: public McMaster
{
  public:
    McStandardSampler();
    virtual ~McStandardSampler();
    virtual void init(bool restart=false);
    virtual void initMaster(bool restart=false);
    virtual void restartMaster();

    virtual void runMaster();

    virtual void updateChain(McChain& chain, McTaskInfo& result);
    virtual void stepTaken(McChain& chain, McTaskInfo& step);
    virtual void stepRejected(McChain& chain, McTaskInfo& step);
    virtual void computeNextStepForChain(McChain& chain);

    virtual void setupLogFiles(bool append=false);
    virtual void setupChains(bool restart=false);
    virtual void cleanupLogFiles();

    virtual void adjustAdaptiveStepSize(McChain& chain, bool stepTaken);
    virtual void updateGelmanRubinStatistics();
    virtual void checkConvergence(unsigned int minChainSize, const std::vector<double>& R);
    virtual void logGelmanRubinStatistics(const std::vector<double>& mean, const std::vector<double>& B,
                                  const std::vector<double>& W, const std::vector<double> R,
                                  bool toProgressLogOnly=true, unsigned int minSize = 0,
                                  unsigned int startPoint=0);
    virtual void logProgress(const McChain& chain);

    virtual void addStartingPoint(const McTaskInfo& t) { mStartingPoints.push_back(t); }

  protected:
    std::vector<MultiGaussian*> mGauss;
    std::vector<McChain*> mChain;

    std::string mLogDir;
    std::vector<std::ofstream*> mInvestigatedLog;
    std::vector<std::ofstream*> mHeadLog;
    std::vector<std::ofstream*> mChainLog;
    std::vector<std::string>    mCovarianceFileNames;
    std::ofstream *mProgressLog, *mGelmanRubinLog, *mCovarianceLog;
    bool mConverged;

    std::list<McTaskInfo> mStartingPoints;
};

#endif // MCRUNNER_H
