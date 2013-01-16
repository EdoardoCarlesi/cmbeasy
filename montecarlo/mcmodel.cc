// vim:ts=2:sw=2:et
#include "mcmodel.h"
#include "mcrunner.h" // for McTaskInfo

McModel::McModel()
     : mParams(0), mControlPanel(0), mCmbCalc(0)
{
  mControlPanel = &McSettings::self()->controlPanel();

  mCl = new CL;
  mCmbCalc = new CmbCalc();

  string cmbeasydir = ControlPanel::cmbeasyDir();
  string filename = cmbeasydir + "/resources/jlgen.dat";
  mCmbCalc->initjl(filename, 2000 /* max nr of l for cmb*/);
}

McModel::~McModel()
{
  delete mCmbCalc;
  delete mCl;
}

void McModel::setParameters(McTaskInfo& t)
{
  mParams = &t;
}

CL* McModel::cmbSpectra() const
{
  return mCl;
}

McTaskInfo* McModel::currentParameters() const
{
  return mParams;
}

void McModel::outputDebugInfo()
{
  cosmos()->printStatus();
  if (cfg().controlPanel().cmb) {
    mCmbCalc->dumpCl(cosmos()->InitialPower, *mCl, *mControlPanel, "scalarCl.dat" , "tensorCl.dat");
  }
  if (cfg().controlPanel().power_cdm) {
    cosmos()->dumpPower(0, "cdmPower", cosmos()->power_cdm(), cosmos()->z2tau(0));
  }
}

