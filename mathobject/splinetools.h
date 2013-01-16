#ifndef SPLINETOOLS_H
#define SPLINETOOLS_H

#include "spline.h"
#include "cosmos.h"
#include "quintcosmos.h"

//#include <QFileInfo>

#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include <typeinfo>
#include <limits>
#include <functional>

//#include <QFileInfo>

class NewLogStepper
{
  public:
    NewLogStepper( const double start, const double end)
      : mLogStart(log(start)), mLogEnd(log(end)),
      mEnd(end), mLogX(log(start)), mDone(false) {
        if (start==0) {
          throw Bad_Error("LogSteps can't start at 0");
        }
      }
    NewLogStepper& setSteps(const int steps) {
      mLogStep = (mLogEnd-mLogStart)/((double)steps-1.);
      return *this;
    }

    struct DenseSamplingRegion {
      double xMin, xMax, stepSize;
      DenseSamplingRegion(double min, double max, double s): xMin(min), xMax(max), stepSize(s) {}
    };
    struct SampleContains {
      SampleContains(const double x): mX(x) {}
      bool operator()(const DenseSamplingRegion& r) const { return ((r.xMin<=mX) && (mX<=r.xMax)); }
      const double mX;
    };
    typedef std::vector<DenseSamplingRegion>::const_iterator RegIt;

    double next() {
      const double arg = mLogX;

      RegIt curRegion = find_if (mDenseRegions.begin(), mDenseRegions.end(),
                                 SampleContains(arg));
      if (curRegion != mDenseRegions.end()) {
        mLogX += curRegion->stepSize;
      } else {
        mLogX += mLogStep;
      }
      if (arg >= mLogEnd) {
        mDone = true;
        return mEnd;
      }
      return exp(arg);
    }

    NewLogStepper& denseSamplingBetween(const double x1, const double x2, unsigned int steps) {
      double logSteps = (log(x2)-log(x1))/(double(steps)-1.);
      mDenseRegions.push_back(DenseSamplingRegion(log(x1), log(x2), logSteps));
      return *this;
    }
    bool hasNext() { return !mDone; }

  private:
    double mLogStart, mLogEnd, mEnd, mLogX, mLogStep;
    bool mDone;
    std::vector<DenseSamplingRegion> mDenseRegions;
};

namespace SplineCompare
{
  static bool allComparisonsVerbose = false;
  static bool fuzzyShapeCheck = false;
  static std::vector<std::string> dumpList;
  //! returns true if s1 and s2 have the same x-values


  void checkDumpList(const Spline& s1, const Spline& s2, bool verbose=false);
  bool haveSameShape(const Spline& s1, const Spline& s2, bool verbose=false);
  //! returns the maximum of fabs(s1.y(i)-s2.y(i))
  double greatestAbsoluteDifference(const Spline& s1, const Spline& s2, bool verbose=false,
                                    int lastx = std::numeric_limits<int>::max());

  //! same as the above, but s2 is loaded from disc; s2 is assumed to have been
  //! written to disc by Spline::save, not dumped to disc by Spline::dump()
  double greatestAbsoluteDifference(const Spline& s1, const std::string& s2Name, bool verbose=false,
                                    int lastx = std::numeric_limits<int>::max());

  //! returns the maximum of fabs(s1.y(i)-s2.y(i))/fabs(s1.y(i))
  double greatestRelativeDifference(const Spline& s1, const Spline& s2, bool verbose=false,
                                    int lastx = std::numeric_limits<int>::max());

  //! same as the above, but s2 is loaded from disc; s2 is assumed to have been
  //! written to disc by Spline::save, not dumped to disc by Spline::dump()
  double greatestRelativeDifference(const Spline& s1, const std::string& s2Name, bool verbose=false,
                                    int lastx = std::numeric_limits<int>::max());
}

namespace SplineTools
{

  struct ToSplineSaver
  {
    public:
      ToSplineSaver() {}
      void saveToSpline(const double x, const double y, const std::string& name);
      void operator()(const double x, const double y, const std::string& name) {
        saveToSpline(x, y, name);
      }
      void reset();
      ~ToSplineSaver();

    private:
        ToSplineSaver(const ToSplineSaver&) {}

        std::map<std::string, Spline*> mSplineMap;
        // no Anchor, since the staticTracer lives until application exit, and we
        // want to avoid the otherwise necessary c++ trickery for deleting global static non-POD things:
        // Anchor mAnchor;
  };

  enum  YAccessMethod { FastY, SafeY };  //! whether to use Spline::fastY() or Spline::safeY for accessing the spline read fro disc
  double fromSplineOnDisk(const double x, const std::string& filename, YAccessMethod m = FastY);
  static ToSplineSaver staticTracer;
}

Spline* tau2a(Cosmos& c, const Spline& s);
Spline* tau2lna(Cosmos& c, const Spline& s);
Spline* tau2z(Cosmos& c, const Spline& s);
Spline* tau2zp1(Cosmos& c, const Spline& s);
Spline* tau2zp1(Cosmos& c, const Spline& s);

void convertFromTau2zp1(Cosmos& c, const std::vector<Spline*>& v);

void dumpOmegas(Cosmos& c, std::string prefix="");
void dumpDensities(Cosmos& c, std::string prefix = "");

#endif // SPLINETOOLS_H
