#include "perturbationtracker.h"
#include "cninvariant.h"
#include "spline.h"
#include "cmbcalc.h"
#include "speedydeinvariant.h"

#include <algorithm>

using namespace std;

static const int numWebs = 12;

PerturbationTracker* PerturbationTracker::mSelf=NULL;

PerturbationTracker* PerturbationTracker::self(Perturbation* p)
{
  if ( !mSelf )
    mSelf = new PerturbationTracker(p);
  mSelf->setPerturbationInstance(p);
  /*  if (p) {
    if (p->k != mSelf->mLastK) {
      mSelf->mNextTauIdx=0;
      mSelf->mCurrentTau=0;
    }
    mSelf->mLastK = p->k;
  }*/
  return mSelf;
}

PerturbationTracker::PerturbationTracker(Perturbation* p)
           : mLastTau(0), mLastK(0)
{
  mPerturbations = p;
  if (p) {
    mLastK = p->k;
  } else {
    mLastK = 0;
  }
  cout << "Created PerturbationTracker(): " << endl;
 allocateSplineWebs();
  mDeltaCdmNl= new Spline(1000, "DeltaCdmNl", &mSplineAnchor);
  mDeltaNuNl = new Spline(1000, "DeltaNuNl", &mSplineAnchor);
  mDeltaPhiNl = new Spline(1000, "DeltaPhiNl", &mSplineAnchor);

  nextTauIndex = 0;
  firstK = true;

  mTrackExtraSplines = mTrackExtraSplineWebs = false;
  mNextTauIdx = 0;
  mCurrentTau = 0;
}

PerturbationTracker::~PerturbationTracker()
{
}

void PerturbationTracker::setPerturbationInstance(Perturbation* p)
{
  mPerturbations = p;
}

void PerturbationTracker::allocateSplineWebs()
{
	//cout << " Allocating Splines... " << endl;
     	mDeltaCdm = new SplineWeb("DeltaCdm", &mSplineAnchor, 1000, 1000);
  mVCdm = new SplineWeb("VCdm", &mSplineAnchor,  1000, 1000);
  mDeltaB = new SplineWeb("DeltaB", &mSplineAnchor, 1000, 1000);
  mVB = new SplineWeb("VB", &mSplineAnchor,  1000, 1000);
  mDeltaPhi = new SplineWeb("DeltaPhi", &mSplineAnchor,  1000, 1000);
  mVQ = new SplineWeb("Vq", &mSplineAnchor,  1000, 1000);
  mDeltaNuNR = new SplineWeb("DeltaNuNr", &mSplineAnchor,  1000, 1000);
  mVNuNr = new SplineWeb("VnuNr", &mSplineAnchor,  1000, 1000);

  mPhi = new SplineWeb("Phi", &mSplineAnchor,  1000, 1000);
  mPhiDot = new SplineWeb("PhiDot", &mSplineAnchor,  1000, 1000);
  mPsi = new SplineWeb("Psi", &mSplineAnchor,  1000, 1000);
  mPsiDot = new SplineWeb("PsiDot", &mSplineAnchor,  1000, 1000);
  
}

void PerturbationTracker::disarmWebs()
{
  SplineWeb* webs[numWebs] = { mDeltaCdm, mVCdm, mDeltaB, mVB, mPhi, mPhiDot,
                               mPsi, mPsiDot, mDeltaPhi, mVQ,
                               mDeltaNuNR, mVNuNr
                             };
  for_each( webs, webs+numWebs, mem_fun(&SplineWeb::disarm));
}

//###fixme: Factor out duplicate code with DumpForK below...
class CreateForK
{
public:
  CreateForK(const double k, const string& prefix, vector<Spline*>& vec)
        : mk(k), mPrefix(prefix), mVec(vec) {}
  void operator() (SplineWeb* web) {
    if (web->ySize() == 0) { // if we've calculated the perturbations only for one single k
      web->ySplineAt(0)->dump(mPrefix+"-"+web->name);
      return;
    }
    Spline* s = web->createAlongX(mk, mPrefix + "-" + web->name);
    mVec.push_back(s);
  }

private:
  double mk;
  const string& mPrefix;
  vector<Spline*>& mVec; // the Vector we want to append all the splines to
};

vector<Spline*> PerturbationTracker::perturbationSplines(const double k, const string& prefix)
{
  vector<Spline*> pertubationSplines;
  SplineWeb* webs[numWebs] = { mDeltaCdm, mVCdm, mDeltaB, mVB, mPhi, mPhiDot,
    mPsi, mPsiDot, mDeltaPhi, mVQ,
    mDeltaNuNR, mVNuNr
  };
  for_each( webs, webs+numWebs, CreateForK(k, prefix, pertubationSplines));

  return pertubationSplines;
}

class DumpForK
{
public:
  DumpForK(const double k, const string& prefix): mk(k), mPrefix(prefix) {}
  void operator()(SplineWeb* web) {
    if (web->xSize() <= 0 && web->ySize() <= 0) {
      cout << "perturbationtracker.cc: DumpForK::operator() : tried to write empty SplineWeb "
                                                              << web->name << " to disk" << endl;
      return;
    }
    if (web->ySize() == 0) { // if we've calculated the perturbations only for one single k
      web->ySplineAt(0)->dump(mPrefix+"-"+web->name);
      return;
    }
    Spline* s = web->createAlongX(mk, "created by DumpForK from " + web->name);
    s->dump(mPrefix+"-"+web->name);
    delete s;
  }

private:
  double mk;
  const string& mPrefix;
};

class DumpAtTau
{
public:
  DumpAtTau(const double tau, const string& prefix) : mTau(tau), mPrefix(prefix) {}
  void operator() (SplineWeb* web) {
    Spline* s = web->createAlongY(mTau, "created by DumpAtTau from " + web->name);
    s->dump(mPrefix+"-"+web->name);
    delete s;
  }

private:
  double mTau;
  const string& mPrefix;
};


void PerturbationTracker::dumpPerturbationsForK(const double k, string& prefix)
{
  SplineWeb* webs[numWebs] = { mDeltaCdm, mVCdm, mDeltaB, mVB, mPhi, mPhiDot,
    mPsi, mPsiDot, mDeltaPhi, mVQ,
    mDeltaNuNR, mVNuNr
  };
  for_each( webs, webs+numWebs, DumpForK(k, prefix)); // original method

  // modified method from now on - override DumpForK
  /*for(int i=0; i<numWebs; i++) 
{
    if (webs[i]->xSize() <= 0 && webs[i]->ySize() <= 0) { 
	cout << "perturbationtracker.cc: DumpForK::operator() : tried to write empty SplineWeb " << 
							webs[i]->name << " to disk" << endl;
	return;		    }
	cout << " SplineWebs[" << i << "]" << endl;
	webs[i]->ySplineAt(0)->dump(prefix+"-"+webs[i]->name);
}
	DumpForK(k,prefix, webs[i]);
*/
}

void PerturbationTracker::dumpPerturbationsAtTau(const double tau, string& prefix)
{
  SplineWeb* webs[numWebs] = { mDeltaCdm, mVCdm, mDeltaB, mVB, mPhi, mPhiDot,
    mPsi, mPsiDot, mDeltaPhi, mVQ,
    mDeltaNuNR, mVNuNr
  };
#warning fix this finally (createAlongY not working correctly)
  Spline* s = mDeltaCdm->ySplineAt(0);
  double closeTau = s->start();
  int i;
  for (i = s->first(); i < s->last() && s->x(i) < tau; ++i) {
      closeTau = s->x(i);
  }
  if (fabs(tau-closeTau) > (tau-s->x(i)))
      closeTau = s->x(i);
  cout << "PerturbationTracker::dumpPerturbationsAtTau - wanted tau: " << tau
    << " got tau " << closeTau << " rel diff= " << ((tau-closeTau)/tau) << endl;
  for_each( webs, webs+numWebs, DumpAtTau(closeTau, prefix));
}

void PerturbationTracker::trackNonlinear(const double z, const double k, const double h, const double delta_c,
                                                   const double delta_nu, const double delta_q)
{
  static const double amplitude = sqrt(2.3476501859e-9);
  static const double threshold = 1;
  const double koh = k/h;
  if (fabs(delta_c)*amplitude > threshold) {
    if (mDeltaCdmNl->size()==0  || (mDeltaCdmNl->stop()!=koh)) {
//X       cout << "set cdm at z=" << z << endl;
      mDeltaCdmNl->set(koh, z);
    }
  }
  if (fabs(delta_nu)*amplitude > threshold) {
    if (mDeltaNuNl->size()==0  || (mDeltaNuNl->stop()!=koh)) {
//X       cout << "set nu at z=" << z << endl;
      mDeltaNuNl->set(koh, z);
    }
  }
  if (fabs(delta_q)*amplitude > threshold) {
    if (mDeltaPhiNl->size()==0  || (mDeltaPhiNl->stop()!=koh)) {
//X       cout << "set phi at z=" << z << endl;
      mDeltaPhiNl->set(koh, z);
    }
  }
}

void PerturbationTracker::dumpNonlinearPlot(string prefix)
{
  if (mDeltaCdmNl->size() > 2) {
    mDeltaCdmNl->arm();
    mDeltaCdmNl->dump(prefix+"first-nl-cdm", false);
  }
  if (mDeltaNuNl->size() > 2) {
    mDeltaNuNl->arm();
    mDeltaNuNl->dump(prefix+"first-nl-nu", false);
  }
  if (mDeltaPhiNl->size() > 2) {
    mDeltaPhiNl->arm();
    mDeltaPhiNl->dump(prefix+"first-nl-phi", false);
  }
}

void PerturbationTracker::track(const string& var, const double tau, const double value)
{
  if (mTausToTrack.empty())
    return;

  if (tau == mCurrentTau) {
    ; // good, track this variable too
  } else if (tau < mCurrentTau || tau<mTausToTrack[mNextTauIdx]) {
    // don't track anything at this tau, wait for tau to reach the value of
    // the next entry in mTausToTrack
    return;
  } else {
    // we have a new tau at which  we do want to note size of perturbations
    mCurrentTau=tau;
    ++mNextTauIdx;
    //cout << "tracking " << " at tau=" << tau << " (idx=" << mNextTauIdx
    //  << " of " << mTausToTrack.size() << ")" << endl;
  }
  //cout << "tracking " << var << "=" << value << " at tau=" << tau << " (idx=" << mNextTauIdx
  //     << " of " << mTausToTrack.size() << ")" << endl;
  if (mTrackExtraSplines) {
    trackExtraSplines(var, tau, value);
  } else if (mTrackExtraSplineWebs) {
    trackExtraSplineWebs(var, tau, value);
  }
}

void PerturbationTracker::trackExtraSplines(const string& var, const double tau, const double value)
{
  map<string, Spline*>::iterator result;
  result = mExtraTrackingSplines.find(var);
  Spline* s = 0;
  if (result == mExtraTrackingSplines.end()) {
    s = new Spline(mTausToTrack.size(), var, &mSplineAnchor);
    mExtraTrackingSplines.insert(make_pair(var, s));
  }
  if (!s) {
    s = result->second;
  }
  //if (s->size()<1e2 || s->stop() < 0.99*tau)
  if (s->size()==0 or tau > s->stop())
    s->set(tau, value);
}

void PerturbationTracker::trackExtraSplineWebs(const string& var, const double tau, const double value)
{
  if (!mPerturbations) {
    throw Bad_Error("PerturbationTracker::trackExtraSplineWebs() - no Perturbation set, needed to get value of k.");
  }
  map<string, SplineWeb*>::iterator result;
  result = mExtraTrackingSplineWebs.find(var);
  SplineWeb* w = 0;
  if (result == mExtraTrackingSplineWebs.end()) {
    w = (SplineWeb*)0;
    mExtraTrackingSplineWebs.insert(make_pair(var, w));
    //w->set(mPerturbations->k, tau, value);
    //cout << w << " " << w->xSize() << " " << w->ySize() << endl;
    //cout << w << " " << w->igrid << " " << w->jgrid << endl;
    //cout << var << ": new web at k=" << mPerturbations->k << " at tau=" << tau << endl;
    //return;
  } else {
    w = result->second;
  }

  //cout << w << " " << w->xSize() << " " << w->ySize() << endl;
  double lastTau = 0;
  double curSize = 0;
  /*
  map<double, unsigned int>::iterator it;
  it = w->xgrid.find(mPerturbations->k);
  if (it != w->xgrid.end()) {
    //cout << "out" << endl;
    //cout << var << ": error for tau=" << tau << " at k=" << mPerturbations->k << endl;
    //throw Bad_Error("PerturbationTracker::trackExtraSplineWebs() - couldn't find Spline for k");
    Spline *current =  w->xSplineAt((*it).second);
    if (current) {
      lastTau = current->stop(); // the last tau for which we have put a value in the SplineWeb
      curSize = current->size();
    }
  }
  */

  Spline* current = mSplineVectors[var][mPerturbations->k];
  if (!current) {
    current = mSplineVectors[var][mPerturbations->k] = new Spline(500, "created by PerturbationTracker::trackExtraSplineWebs");
  }
  if (current->size()>0) {
    lastTau = current->stop(); // the last tau for which we have put a value in the SplineWeb
    curSize = current->size();
  }
  static bool warned = false;
  if (!warned) {
    cout << "WARNING: consider sampling a little more points" << endl; warned = true;
    cout << "also: switch off FadeOut..." << endl;
  }
  if (true || lastTau<tau && (curSize<1e2 || lastTau < 0.85*tau)) {
    if (current->size()==0 or tau > current->stop())
      current->set(tau, value);
  } else {
    ;
  }
}

class DumpSecond
{
  public:
    DumpSecond(const string& p) :  mPrefix(p) {}
    void operator ()(pair<string, Spline*> p) {
      p.second->dump(mPrefix + p.second->name);
      cout << "wrote to disc: " << (mPrefix + p.second->name) << " (" << p.second->size() << " points)" << endl;
    }
  private:
    const string& mPrefix;
};

void PerturbationTracker::dumpExtraSplines(const string& prefix)
{
  using namespace std;
  if (!mTrackExtraSplines) {
    cout << "Not tracking extra perturbation quantities - no output generated." << endl;
    return;
  }
  for_each(mExtraTrackingSplines.begin(), mExtraTrackingSplines.end(),
           DumpSecond(prefix));
}

class AppendSpline
{
    vector<Spline*>& mVec;
  public:
    AppendSpline(vector<Spline*>& v) : mVec(v) {}
    void operator ()(pair<string, Spline*> p) {
      mVec.push_back(p.second);
    }
};

vector<Spline*> PerturbationTracker::extraSplines()
{
  using namespace std;
  vector<Spline*> vec;
  for_each(mExtraTrackingSplines.begin(), mExtraTrackingSplines.end(),
           AppendSpline(vec));
  return vec;
}

Spline* PerturbationTracker::extraSpline(const string& name, Anchor* a) const
{
  map<string, Spline*>::const_iterator it;
  it = mExtraTrackingSplines.find(name);
  if (it == mExtraTrackingSplines.end()) {
    throw Bad_Error("PerturbationTracker::extraSpline(string): no spline named " + name + " tracked. ");
  }
  Spline *s = new Spline(*(*it).second, (*it).second->name, a);
  return s;
}

class AppendSplineWeb
{
  vector<SplineWeb*>& mVec;
  public:
  AppendSplineWeb(vector<SplineWeb*>& v) : mVec(v) {}
  void operator ()(pair<string, SplineWeb*> p) {
    mVec.push_back(p.second);
  }
};

static bool convertedToSplineWeb = false;

static SplineWeb* createWeb(const string& name,  const vector<double>& taus, Anchor* a,
                            const PerturbationTracker::SplineMap& theSplines,
                            map<string, SplineWeb*>&        theWebs)
{
  map<string, SplineWeb*>::const_iterator it = theWebs.find(name);
  if (it == theWebs.end()) {
    throw Bad_Error("perturbationtracker.cc: createWeb()- no splineWeb named " + name + " tracked. ");
  }
  SplineWeb *w = (*it).second;
  if (w) {
    if (w->MyAnchor != a) {
      throw Bad_Error(" perturbationtracker.cc: createWeb() - called twice, but with different Anchors");
    }
    return w;
  }

  typedef map<double, Spline*> SplinePMap;
  const SplinePMap& splines = (*(theSplines.find(name))).second;
  SplinePMap::const_iterator it1, end1 = splines.end();
  it1 = splines.begin();

  unsigned int kCount = splines.size();
  unsigned int tauCount = taus.size();

  w = new SplineWeb(name, a, kCount, tauCount);

  for (; it1 != end1; ++it1) {
    double k = it1->first;
    Spline* s = it1->second;
    if (!s->armed) {
      s->arm();
    }
    for (int i=0; i<taus.size(); ++i) {
      double tau = taus[i];
      //if (tau < s->start())
      //  continue;
#warning DONT use safeY here!!
      double value = s->safeY(tau);
      w->set(k, tau, value);
    }
    delete s;
  }
  theWebs[name] = w;
  return w;
}


vector<SplineWeb*> PerturbationTracker::extraSplineWebs(const vector<double>& taus, Anchor* a)
{
  vector<SplineWeb*> v;
  SplineMap::iterator it, end = mSplineVectors.end();
  for ( it = mSplineVectors.begin(); it != end; ++it) {
     v.push_back(createWeb(it->first, taus, a, mSplineVectors, mExtraTrackingSplineWebs));
  }
  return v;
}

SplineWeb* PerturbationTracker::extraSplineWeb(const string& name, const vector<double>& taus, Anchor* a)
{
  SplineWeb* w = createWeb(name, taus, a, mSplineVectors, mExtraTrackingSplineWebs);
  return w;
}

void PerturbationTracker::reset()
{
  mSplineAnchor.kill();
  mExtraTrackingSplines.clear();
  allocateSplineWebs();
}

struct TrackerHelper
{
    TrackerHelper(const PerturbationTracker* t, const double tau): mTracker(t), mTau(tau) {}

    void operator()(const Perturbation* const pert) {
      if (!pert)
        return;

      const double h = pert->cosmos->h();
      const double k = pert->k;
      const double kh = k/h;
      const SpeedyInvariant* const sInv = dynamic_cast<const SpeedyInvariant* const>(pert);
      if (sInv) {
        ;
      }
      /*
      const CnInvariant* const cnInv = dynamic_cast<const CnInvariant* const>(pert);
      if (cnInv) {
        mTracker->mDeltaCdm->set(mTau, kh, cnInv->delta_c_longit());
        mTracker->mVCdm->set(mTau, kh, cnInv->v_c_longit());
        mTracker->mDeltaB->set(mTau, kh, cnInv->delta_b_longit());
        mTracker->mVB->set(mTau, kh, cnInv->v_b_longit());
        mTracker->mPhi->set(mTau, kh, cnInv->phi());
        mTracker->mPhiDot->set(mTau, kh, cnInv->phiDot());
        mTracker->mPsi->set(mTau, kh, cnInv->psi());
        mTracker->mPsiDot->set(mTau, kh, cnInv->psiDot());
        mTracker->mDeltaPhi->set(mTau, kh, cnInv->delta_q_longit());
        mTracker->mVQ->set(mTau, kh, cnInv->v_q_longit());
        mTracker->mDeltaNuNR->set(mTau, kh, cnInv->delta_nuNr_longit());
        mTracker->mVNuNr->set(mTau, kh, cnInv->v_nu_nr_longit());
        //cout << "\n set (" << cnInv->k << ", " << mTau << ")\n";
      }
      */
    }

  private:
    const PerturbationTracker* mTracker;
    const double mTau;
};

#warning remove
void PerturbationTracker::trackPerturbationQuantities(const double tau)
{
//X #warning disabled tracking
//X   return;
//X #warning fixme
//X   if (tau <= 7e-7)
//X     return;
//X #warning fixme, otherwise this class is useless
//X   double z = pert->tau2z(tau);
//X   static double lastz = 100;
//X   ofstream *ofs;
//X   ofs = 0;
//X   if (lastz>50 && z<=50)
//X     ofs = z50f;
//X   else if (lastz>5 && z<=5)
//X     ofs = z5f;
//X   else if (lastz>0.5 && z<=0.5)
//X     ofs = z05f;
//X   else if (lastz>0.01 && z<=0.01)
//X     ofs = z0f;
//X   if (ofs) {
//X     cout << ofs << endl;
//X     (*ofs) << pert->k/pert->h() << "      " << pert->delta_c_longit()
//X       << "      " << pert->phi()
//X       << "      " <<  pert->delta_nu_nr_longit()
//X         << endl;
//X    lastz = pert->tau2z(tau);
//X   }

  const double numPoints = mDeltaCdm->xSize();
  static int lastTauIndex = 0;
  if (!mPerturbations) {
    throw Bad_Error("PerturbationTracker::trackPerturbationQuantities() - no perturbation pointer set.");
  }
  if (mLastK != mPerturbations->k ) {
    cout << "used " << nextTauIndex << " entries for " << mLastK << endl;
    mLastK = mPerturbations->k;
    mLastTau = 0;
//X     lastz=100;
    nextTauIndex=0;
    if (mTaus.size()>0)
      firstK = false;
    cout << firstK << " switched k, new one is " << mPerturbations->k << " at tau= " << tau << endl;
  }
  if ( firstK /*&&  numPoints > 10*/ &&  tau < 1.05*mLastTau ) {
    return;
  }
//X   static double lastStatusOutputTau = 0;
//X   if( lastStatusOutputTau == 0 or tau>100.+lastStatusOutputTau) {
//X     lastStatusOutputTau = tau;
//X     //const double  z = pert->tau2z(tau);
//X      cout << "tracked to tau=: " <<  tau << "     \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << flush;
//X   }

#warning check accuracy of this, i.e. that we track things close to their real tau...
  double nextTau;
  if (!firstK) {
    nextTau = mTaus[nextTauIndex];
    if (tau>=nextTau && (nextTauIndex+1<mTaus.size())) {
      ++nextTauIndex;
    } else {
        return;
    }
  } else {
    if (mTaus.size() > 0 && mTaus.back() >= tau)
      return;
    mTaus.push_back(tau);
//X     cout << "index: " << nextTauIndex << " for tau= " << tau << endl;
    ++nextTauIndex;
  }



#warning fix this: first k always has one more tau than the others...
//X     if (nextTauIndex>=219) return;  // for Christof's model
//X     if (nextTauIndex>=119) return;  // for lcdm reference
    if (nextTauIndex>=186) return;  // for lcdm reference
    if (nextTauIndex>=1142) return;  // for Christof's model

  mLastTau = tau;

  TrackerHelper th(this, (firstK?tau:nextTau));
  th(mPerturbations);

  //cout << "tracked at: " << (firstK?tau:nextTau) << endl;


  if(false && tau > 300)  {//false && numPoints>2000) {
    string prefix = "debug-pert-";
    dumpPerturbationsForK(mPerturbations->k, prefix);
    throw Bad_Error("stop here");
  }
}

