#ifndef FASTSLOWSTEPPER_H
#define FASTSLOWSTEPPER_H

#include "mcrunner.h"
#include "mcmodel.h"

class FastSlowStepper: public McStandardSampler
{
  public:
    FastSlowStepper();
    virtual ~FastSlowStepper();
    virtual void setupChains(bool restart);
    virtual void computeNextStepForChain(McChain& chain);
    virtual void stepTaken(McChain& chain, McTaskInfo& step);
    virtual void adjustAdaptiveStepSize(McChain& chain, bool tookStep);

  protected:
    vector<int> mFastStepsTaken;
    vector<bool> mLastProposalWasFast, mIsFastParam;
    std::vector<MultiGaussian*> mFastGauss;
    typedef std::vector<std::vector<double> > CovMatrix;
    std::vector<CovMatrix> mFastCovariance;
};

class FastSlowSlave: public McRunner
{
  virtual void runSlave();
};

class CacheItem
{
  public:
    virtual CL*     stealRawCmbSpectra()=0;
    virtual Cosmos* stealCosmos()=0;
    virtual void    setRawCmbSpectra(CL* cl)=0;
    virtual void    setCosmos(Cosmos* c)=0;
    virtual void    recomputeFast()=0;
};

template<class CosmosClass>
class CacheableModel: public McCustomModel<CosmosClass>, public CacheItem
{
  public:
    CacheableModel(): mRawCl(new CL()) {}
    virtual CL*          rawCmbSpectra() const { return mRawCl; }
    virtual CL*          stealRawCmbSpectra() {
      CL* tmp=mRawCl;
      mRawCl=new CL();
      return tmp;
    }
    virtual void         setRawCmbSpectra(CL* cl) { mRawCl=cl; }
    virtual CosmosClass* stealCosmos() {
      CosmosClass* tmp=McCustomModel<CosmosClass>::mCosmos;
      setCosmos(new CosmosClass());
      return tmp;
    }
    virtual void setCosmos (Cosmos* c) {
      McCustomModel<CosmosClass>::setCosmos(c);
    }

  protected:
    CL* mRawCl;
};
#endif // FASTSLOWSTEPPER_H
