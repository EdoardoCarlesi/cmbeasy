// vim: sw=2 ts=2 et

#include "cl.h"

#include "spline.h"

#include <algorithm>

using namespace std;

struct CloneSpline {
  Spline* operator()(const Spline* s) {
    return new Spline(*s, s->name, s->MyAnchor);
  }
};

struct CloneSplineWithNewMother {
  Spline* newMother;
  CloneSplineWithNewMother(Spline* m): newMother(m) {}
  Spline* operator()(Spline* s) {
    Spline *newSpline = new Spline(newMother, s->name, s->MyAnchor);
    for (int i=s->first(); i<=s->last(); ++i) {
      if (s->x(i)!=newMother->x(i)) {
        throw Bad_Error("cl.cc: CloneSplineWithNewMother(): old and new mother have differenct");
      }
      newSpline->setY(i, s->y(i));
    }
    return newSpline;
  }
};

CL* CL::clone()
{
  CL* clonedCL = new CL;

  SplineCleanVector* spectra[] = {/*no &ts, is special case*/ &(clonedCL->tt), &(clonedCL->es), &(clonedCL->et),
                                   &(clonedCL->bt), &(clonedCL->cs), &(clonedCL->ct), &(clonedCL->kk), &(clonedCL->tk),
                                   &(clonedCL->bs)};
  const SplineCleanVector* mySpectra[] = {/*no &ts, is special case*/ &tt, &es, &et,
                                       &bt, &cs, &ct, &kk, &tk, &bs};
  const int count = sizeof(spectra)/sizeof(SplineCleanVector*);

  transform(ts.begin(), ts.end(), back_inserter(clonedCL->ts), CloneSpline());
  for (int i=0; i<count; ++i) {
    transform(mySpectra[i]->begin(), mySpectra[i]->end(), back_inserter(*spectra[i]),
              CloneSplineWithNewMother(clonedCL->ts[0]));
  }
  return clonedCL;
}

void CL::resize(unsigned int n) {
  clear();

  ts.resize(n);
  tt.resize(n);
  es.resize(n);
  et.resize(n);
  bt.resize(n);
  cs.resize(n);
  ct.resize(n);
  kk.resize(n);
  tk.resize(n);
  bs.resize(n);
}

void CL::clear() {
  ts.clear(); tt.clear(); es.clear(); et.clear();
  bt.clear(); cs.clear(); ct.clear(); kk.clear(); tk.clear();
  bs.clear();
}

/*!
  Multiply all splines for spectral index number i with
  factor x
  */
void CL::mul(const double x, const int i) {
  *ts[i] *= x;
  *tt[i] *= x;
  *es[i] *= x;
  *et[i] *= x;
  *bt[i] *= x;
  *cs[i] *= x;
  *ct[i] *= x;
  *kk[i] *= x;
  *tk[i] *= x;
  *bs[i] *= x;
}


Spline* CL::createTotalTTSpline(const int i, const std::string name, Anchor* a) const
{
  return createTotalXXSpline(i, ts, tt, name, a);
}

Spline* CL::createTotalEESpline(const int i, const std::string name, Anchor* a) const
{
  return createTotalXXSpline(i, es, et, name, a);
}

Spline* CL::createTotalTESpline(const int i, const std::string name, Anchor* a) const
{
  return createTotalXXSpline(i, cs, ct, name, a);
}

Spline* CL::createTotalBBSpline(const int i, const std::string name, Anchor* a) const
{
  return createTotalXXSpline(i, bs, bt, name, a);
}

Spline* CL::createTotalXXSpline(const int i, const SplineCleanVector& scalar, const SplineCleanVector& tensor,
                                const std::string name, Anchor* a) const
{
  Spline* totalCls = new Spline(*scalar[i], name, a);
  Spline* tensorCls = tensor[i];
  if (tensorCls->last()>totalCls->last()) {
    throw Bad_Error("CL::createTotalXXSpline() - tensor spline has larger l_max than scalar spline.");
  }

  for (int i=tensorCls->first(); i<=tensorCls->last(); ++i) {
    if (totalCls->x(i) != tensorCls->x(i)) {
      throw Bad_Error("CL::createTotalXXSpline() - index mismatch between tensor and scalar splines.");
    }
    totalCls->add(i, tensorCls->y(i));
  }
  return totalCls;
}
