#ifndef MCLIKELIHOODCALCULATOR_H
#define MCLIKELIHOODCALCULATOR_H
// vim:ts=2:sw=2:et

#include "mcmodel.h"

#include "analyzethis.h"

#include <list>

#define REGISTER_DATASET(x) static DataSetGenerator<x> x##Generator(#x); \
                            static x::Configurator& configure##x() { return dynamic_cast<x*>(McSettings::self()->dataSet(#x))->config(); }

class DataSet
{
  public:
    virtual void initialize() = 0;
    virtual void computeLogLike(const McModel& model, McTaskInfo& result) = 0;

    virtual ~DataSet() {}
    void setAnalyzeThis(AnalyzeThis& ai) { mAnalyzeThis = &ai; }

    virtual void needCmb() const {
      McSettings::self()->controlPanel().cmb=true;
    }

    virtual void needTransfer() const {
      McSettings::self()->controlPanel().power_cdm=true;
    }

    typedef struct ThisDataSetIsNotConfigurable {} Configurator;
    Configurator& config() { return mDummyConfigurator; }

  private:
    Configurator mDummyConfigurator;

  protected:
    AnalyzeThis* mAnalyzeThis;
};

template<class DataS>
struct DataSetGenerator
{
  DataSetGenerator<DataS>(string name) {
    McSettings::registerDataSet(name, DataSetGenerator<DataS>::create);
  }
  static DataSet* create() {
    DataSet* t = new DataS();
    t->setAnalyzeThis(McSettings::self()->analyzeThis());
    return t;
  }
};

class McLikelihoodCalculator
{
  public:
    static void useData(DataSet* data);
    static void initialize();
    double computeLogLike(const McModel& model, McTaskInfo& result);

  private:
    static std::list<DataSet*> mDataSets;
};

#endif // MCLIKELIHOODCALCULATOR_H
