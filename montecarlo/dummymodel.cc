#include "mcmodel.h"

#include "gauge.h"
#include "lensing.h"
#include "analyzethis.h"
#include "quintcosmos.h"

#include "fastslowstepper.h"
#include "nestedsampler.h"
#include "mclikelihoodcalculator.h"

#include <algorithm>
#include <numeric>
#include <sstream>

using namespace std;

class DummyModel: public McCustomModel<Cosmos>
{
  void compute();
};

class NdGaussian
{
  double gaussExponent(double x, double center, double sigma) {
    double distSq = pow(x-center, 2.);
    double chi2 = distSq/(2.*sigma*sigma);
    return chi2;
  }

  public:
  double logLike(McTaskInfo& result) {
    vector<OneDGauss>::const_iterator it, end = mGauss.end();
    double lnLike = 0;
    for (it=mGauss.begin(); it!=end; ++it) {
      double sigma=it->sigma;
      double x=result(it->name);
      lnLike+=-gaussExponent(x, it->mean, it->sigma);
    }
    return lnLike;
  }

  void addDimension(const std::string& name, double mean, double sigma) {
    mGauss.push_back(OneDGauss(name, mean, sigma));
  }

  struct OneDGauss
  {
    OneDGauss(std::string n, double m, double s): name(n), mean(m), sigma(s) {}
    std::string name;
    double mean, sigma;
  };

  std::vector<OneDGauss> mGauss;
};

class DummyData: public DataSet
{
  public:
    void initialize() {
      McTaskInfo::addEntry(McTaskInfo::LogLike, "DummyLogLike");
    }

    void computeLogLike(const McModel& model, McTaskInfo& result) {
      result("DummyLogLike") = mNdGauss.logLike(result);
    }

    typedef NdGaussian Configurator;
    Configurator& config() { return mNdGauss; }

  private:
    NdGaussian mNdGauss;
};
REGISTER_DATASET(DummyData)

void McModel::initialize()
{
  //NestedSampler* master = new NestedSampler();
  //cfg().setMaster(master);
  //master->setLivePointsCount(300);
  //master->setMinimumIterations(1e3);
  //master->setMaxRemainingLogEvidenceFraction(0.1);
  //master->setSamplingMethod(NestedSampler::MultiNest);
  //master->setEllipsoidalEnlargementFactor(1.1);

  cfg().setModel(new DummyModel());

  cfg().useDataSet("DummyData");
  configureDummyData().addDimension("x", 10, 2);
  configureDummyData().addDimension("y", 3, 2);

  // addParameter(name, lowerBound, upperBound, initialSigma)
  McTaskInfo::addMcParameter("x", -10, 30, 2);
  McTaskInfo::addMcParameter("y", -10, 10, 2);
  McTaskInfo::addMcParameter("zeta", 1, 5, 2);

  // add additional things you would like to keep track of
  //McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "sum");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "x/zeta");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "x/zeta/zeta");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "x/zeta/zeta/zeta");
  McTaskInfo::addEntry(McTaskInfo::ExtraInfo, "x/zeta/zeta/zeta/zeta");
}

void DummyModel::compute()
{
  McTaskInfo& param = *mParams;

  double zeta=param("zeta");
  double x=param("x");
  param("x/zeta") = x/zeta;
  param("x/zeta/zeta") = x/zeta/zeta;
  param("x/zeta/zeta/zeta/zeta") = x/zeta/zeta/zeta;
  //param("sum") = param("param1") + param("param2");
}
