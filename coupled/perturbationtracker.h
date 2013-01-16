#ifndef PERTURBATIONTRACKER_H
#define PERTURBATIONTRACKER_H


#include "spline.h"
#include "perturbation.h"

#include <string>

class CmbCalc;

class TrackerHelper;

class PerturbationTracker
{
  public:
    PerturbationTracker(Perturbation* p);
    ~PerturbationTracker();
    static PerturbationTracker* self(Perturbation* p=0);

    void dumpNonlinearPlot(std::string prefix = "");
    void setTrackExtraSplinesEnabled(const bool b) { mTrackExtraSplines = b; }
    void setTrackExtraSplineWebsEnabled(const bool b) { mTrackExtraSplineWebs = b; }
    void dumpExtraSplines(const std::string& prefix);
    void dumpPerturbationsForK(const double k, std::string& prefix);
    void dumpPerturbationsAtTau(const double tau, std::string& prefix);
    Spline* extraSpline(const std::string& name, Anchor* a=0) const;
    std::vector<Spline*> extraSplines();
    SplineWeb* extraSplineWeb(const std::string& name, const vector<double>& taus, Anchor* a=0);
    std::vector<SplineWeb*> extraSplineWebs(const vector<double>& taus, Anchor* a=0);

    void setTausToTrack(const std::vector<double> v) { mTausToTrack = v; }
    std::vector<double> tausTracked() const { return mTausToTrack; }

    void reset();
    void trackPerturbationQuantities(const double tau);
    std::vector<Spline*> perturbationSplines(const double k, const std::string& prefix);

    void trackNonlinear(const double z, const double k, const double h, const double delta_c,
                        const double delta_nu, const double delta_q);

    void track(const std::string& var, const double tau,
               const double value);

    friend class TrackerHelper;

    void disarmWebs();

    SplineWeb* deltaCdmWeb() const { return mDeltaCdm; }
    SplineWeb* vCdmWeb() const { return mVCdm; }
    SplineWeb* deltaBWeb() const { return mDeltaB; }
    SplineWeb* vBWeb() const { return mVB; }
    SplineWeb* phiWeb() const { return mPhi; }
    SplineWeb* phiDotWeb() const { return mPhi; }

    void setPerturbationInstance(Perturbation* p);

  protected:
    void allocateSplineWebs();
    void trackExtraSplines(const std::string& var, const double tau,
                           const double value);
    void trackExtraSplineWebs(const std::string& var, const double tau,
                              const double value);

  private:
    Perturbation* mPerturbations;
    Spline *mDeltaCdmNl, *mDeltaPhiNl, *mDeltaNuNl;

    SplineWeb* mDeltaCdm, *mVCdm;
    SplineWeb* mDeltaB, *mVB;
    SplineWeb* mDeltaPhi, *mVQ;
    SplineWeb* mDeltaNuNR, *mVNuNr;
    SplineWeb  *mPhi, *mPhiDot, *mPsi, *mPsiDot;

    Anchor mSplineAnchor;
    double mLastTau, mLastK;

    std::map<std::string, Spline*> mExtraTrackingSplines;
    std::map<std::string, SplineWeb*> mExtraTrackingSplineWebs;
    // while tracking, we keep pert vars in a map of Splines
    // and create the SplineWebs at the end
  public:
    typedef std::map<std::string /*name of var*/, std::map<double /*k*/, Spline*> > SplineMap;

  private:
    SplineMap mSplineVectors;

    bool mTrackExtraSplines, mTrackExtraSplineWebs;

    // these are the old ones and should be removed
    std::vector<double> mTaus;
    int nextTauIndex;
    bool firstK;

    static PerturbationTracker* mSelf;

    std::vector<double> mTausToTrack;
    unsigned int mNextTauIdx;
    double mCurrentTau;
    double mK;

    friend class Backreaction;
};


#endif // PERTURBATIONTRACKER_H
