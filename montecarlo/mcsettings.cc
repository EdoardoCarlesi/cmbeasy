// vim: sw=2 ts=2 et

#include "mcsettings.h"

#include "mcerror.h"
#include "mcmodel.h"
#include "mclikelihoodcalculator.h"

static McSettings* mcSettingsInstance = 0;

McSettings* McSettings::self()
{
  if (!mcSettingsInstance) {
    mcSettingsInstance = new McSettings();
  }
  return mcSettingsInstance;
}

// static
McSettings::GeneratorMap& McSettings::availableDataSets()
{
  static GeneratorMap* availableDataSetMap = new GeneratorMap();
  return *availableDataSetMap;
}


McSettings::McSettings()
{
  mLogDir = "";
  mNumberOfChains = 4;
  mMinRStatPoints = 30;
  mStatisticsUpdateSteps = 10;
  mEstimateCovariance=false;
  mAdaptiveStepSize = true;
  mFreezeIn = true;
  mFreezeInR = 1.2;
  mMinFreezeInPoints=250;
  mStepIncrease=1.15;
  mStepDecrease=0.9;
  mHighStepBound=5;
  mLowStepBound=3;
  mUseReallyInvestigated=true;
  mUseAverageEntireFactor=true;
  mUpdateTime=5;
  mBeginCovUpdate=100;
  mMaxPreviousCovPoints=5000;
  mFastStepRatio=1;

  //mControlPanel.cmb = mControlPanel.power_cdm = false;

  mMaster = 0;
  mSlave = 0;
  mModel = 0;
}

string McSettings::errorFileName() const
{
  return logDir()+"errorlog.txt";
}

McModel& McSettings::model() const
{
  if (!mModel) {
    throw McError("McSetting::model() - no model set: McSettings::setModel() should be called in McModel::initialize()");
  }
  return *mModel;
}

McMaster& McSettings::master() const
{
  // if not explicitely set, return the default master with adaptive covariance sampling
  if (!mMaster) {
    const_cast<McSettings*>(this)->mMaster = new McStandardSampler();
  }
  return *mMaster;
}

McRunner& McSettings::slave() const
{
  // if not explicitly set, return the default slave
  // with no bells and whistles, i.e. just use the default
  // McRunner as slaves.
  if (!mSlave) {
    const_cast<McSettings*>(this)->mSlave = new McRunner();
  }
  return *mSlave;
}

void McSettings::useDataSet(const std::string& name)
{
  GeneratorMap::iterator it = availableDataSets().find(name);
  if ( it == availableDataSets().end()) {
    string errorStr = "McLikelihoodCalculator::useDataSet() - DataSet ";
    errorStr += name + " not available. DataSets registerd via REGISTER_DATASET() are: ";
    errorStr += namesOfAvailableDataSets();
    throw McError(errorStr);
  }
  DataSet* set=(*it).second(); // calls DataSetGenerator::create()
  McLikelihoodCalculator::useData(set);
  mActiveDataSets.insert(make_pair(name, set));
}

DataSet* McSettings::dataSet(const std::string& name) const
{
  std::map<string, DataSet*>::const_iterator it;
  it = mActiveDataSets.find(name);
  if (it==mActiveDataSets.end()) {
    throw McError("McSettings::dataSet() - dataset '"+name+"' not in use.");
  }
  return (*it).second;
}
