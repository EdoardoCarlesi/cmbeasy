#ifndef MCMODEL_H
#define MCMODEL_H
// vim:ts=2:sw=2:et

#include "mcrunner.h" // for McTaskInfo
#include "mcsettings.h"

#include "controlpanel.h"
#include "cosmos.h"
#include "cmbcalc.h"

class McModel
{
  public:
    McModel();
    ~McModel();

    static void initialize();

    void setParameters(McTaskInfo& t);
    virtual void compute() = 0;
    //void compute();
    //void compute(GridXYZ *g);

    virtual McTaskInfo* currentParameters() const;
    virtual CL*         cmbSpectra() const;
    virtual Cosmos*     cosmos() const = 0;
    virtual void        setCosmos(Cosmos*) = 0;

    virtual void outputDebugInfo();

  protected:
    static McSettings& cfg() { return *McSettings::self(); }

  protected:
    McTaskInfo*   mParams;

    CL*           mCl;
    ControlPanel* mControlPanel;
    CmbCalc*      mCmbCalc;
};

template<class CosmosClass>
class McCustomModel: public McModel
{
  public:
    McCustomModel<CosmosClass>(): mCosmos(new CosmosClass()) {}
    virtual ~McCustomModel<CosmosClass>() { delete mCosmos; }

    void setCosmos(Cosmos* c) {
      CosmosClass *cc = dynamic_cast<CosmosClass*>(c);
      mCosmos=cc;
    }

    CosmosClass* cosmos() const { return mCosmos; }

    typedef CosmosClass CosmosType;

  protected:
    CosmosClass *mCosmos;
};
#endif // MCMODEL_H
