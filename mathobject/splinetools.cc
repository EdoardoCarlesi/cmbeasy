#include "splinetools.h"

// vim: sw=2 ts=2 et

//X static bool SplineCompare::allComparisonsVerbose = false;
//X static bool SplineCompare::fuzzyShapeCheck = false;
//X static std::vector<std::string> SplineCompare::dumpList;

namespace SplineCompare
{
  //! returns true if s1 and s2 have the same x-values


#warning disabled checkDumpList due to Qt dependency
  void checkDumpList(const Spline& s1, const Spline& s2, bool verbose)
  {
    cout << "warning: disabled checkDumpList due to Qt dependency" << endl;
//X     std::string all("all");
//X     if ((    std::find(dumpList.begin(), dumpList.end(), s1.name) == dumpList.end())
//X          && (std::find(dumpList.begin(), dumpList.end(), s2.name) == dumpList.end())
//X          && (std::find(dumpList.begin(), dumpList.end(), all) == dumpList.end()))
//X         return;
//X 
//X     QFileInfo info1(QString::fromStdString(s1.name));
//X     QFileInfo info2(QString::fromStdString(s2.name));
//X     std::string name1 = info1.fileName().toStdString();
//X     std::string name2 = info2.fileName().toStdString();
//X     if (verbose)
//X       cout << "Dumping " << name1 << " and " << name2 << endl;
//X     const_cast<Spline&>(s1).dump(name1);
//X     const_cast<Spline&>(s2).dump(name2);
  }

  bool haveSameShape(const Spline& s1, const Spline& s2, bool verbose)
  {
    if (&s1 == &s2) return true;

    std::string errorString = s1.name + " and " + s2.name + " don't have the same shape - ";
    if (s1.size() != s2.size())
    {
      if (verbose)
       cout << errorString << "sizeError" << endl;
      return false;
    }
    if (s1.start() != s2.start())
    {
     if (verbose)
      cout << errorString << "startError" << endl;
       return false;
    }

    if ( fuzzyShapeCheck
         && (fabs(s1.stop()-s2.stop()) < (30. *std::numeric_limits<double>::epsilon())) )
      return true;
    else if ((s1.stop() != s2.stop()))
    {
      cout << errorString << "stopError" << endl;
      return false;
    }

    const int last = s1.last();
    for (int i = s1.first(); i <= last; ++i)
    {
      if (s1.x(i) != s2.x(i))
      {
        if (verbose)
          cout << errorString << ": different at position: "<< i << endl;
        return false;
      }
    }

    return true;
  }

  //! returns the maximum of fabs(s1.y(i)-s2.y(i))
  double greatestAbsoluteDifference(const Spline& s1, const Spline& s2, bool verbose, int lastx)
  {
    checkDumpList(s1, s2, verbose || allComparisonsVerbose);
    if (!haveSameShape(s1, s2, verbose || allComparisonsVerbose))
      return std::numeric_limits<double>::infinity();

    double maxDiff = 0;
    double maxDiffx = std::numeric_limits<double>::infinity();
    int maxDiffi = -1;
    const int last = s1.last();
    for (int i = s1.first(); i <= last && s1.x(i) <= lastx; ++i) {
      double diff = fabs(s1.y(i)-s2.y(i));
      maxDiff = max(maxDiff, diff);
      if (maxDiff == diff)
      {
        maxDiffx = s1.x(i);
        maxDiffi = i;
      }
    }
    if (verbose || allComparisonsVerbose)
    {
      if (maxDiffi != -1)
      {
        double s1y = s1.y(maxDiffi);
        double s2y = s2.y(maxDiffi);
        cout << "greatest absolute difference of " << maxDiff << " between " << endl
             << "    " << s1.name << " and " << s2.name << endl
             << "    " << "at " << maxDiffx << ", index " << maxDiffi << " of " << s1.last() << endl
             << "    " << "s1.y: " << s1y << ", s2.y: " << s2y
             << endl;
      }
      else
      {
        cout << "greatest absolute difference is zero between " << endl
             << "    " << s1.name << " and " << s2.name << endl;
      }
    }

    return maxDiff;
  }

  //! same as the above, but s2 is loaded from disc; s2 is assumed to have been
  //! written to disc by Spline::save, not dumped to disc by Spline::dump()
  double greatestAbsoluteDifference(const Spline& s1, const std::string& s2Name, bool verbose, int lastx)
  {
    ifstream in(s2Name.c_str());
    Spline s2(in, s2Name /* name of the spline*/);
    return greatestAbsoluteDifference(s1, s2, verbose, lastx);
  }

  //! returns the maximum of fabs(s1.y(i)-s2.y(i))/fabs(s1.y(i))
  double greatestRelativeDifference(const Spline& s1, const Spline& s2, bool verbose, int lastx)
  {
    checkDumpList(s1, s2, verbose || allComparisonsVerbose);
    if (!haveSameShape(s1, s2, verbose || allComparisonsVerbose))
      return std::numeric_limits<double>::infinity();

    double maxDiff = 0;
    int maxDiffi = -1;
    double maxDiffx = std::numeric_limits<double>::infinity();
    const int last = s1.last();
    for (int i = s1.first(); i <= last  && s1.x(i) <= lastx; ++i) {
      double diff = 0;
      if (s1.y(i)!=0) diff = fabs(s1.y(i)-s2.y(i))/fabs(s1.y(i));
      else if (s2.y(i)!=0) diff = fabs(s1.y(i)-s2.y(i))/fabs(s2.y(i));

      maxDiff = max(maxDiff, diff);
      if (maxDiff == diff)
      {
        maxDiffx = s1.x(i);
        maxDiffi = i;
      }
    }

    if (verbose || allComparisonsVerbose)
    {
      if (maxDiffi != -1)
      {
        double s1y = s1.y(maxDiffi);
        double s2y = s2.y(maxDiffi);
        cout << "greatest relative difference of " << maxDiff << " between " << endl
             << "    " << s1.name << " and " << s2.name << endl
             << "    " << "at " << maxDiffx << ", index " << maxDiffi << " of " << s1.last() << endl
             << "    " << "s1.y: " << s1y << ", s2.y: " << s2y
             << endl;
      }
      else
      {
        cout << "greatest relative difference is zero between " << endl
             << "    " << s1.name << " and " << s2.name << endl;
      }
    }

    return maxDiff;
  }

  //! same as the above, but s2 is loaded from disc; s2 is assumed to have been
  //! written to disc by Spline::save, not dumped to disc by Spline::dump()
  double greatestRelativeDifference(const Spline& s1, const std::string& s2Name, bool verbose, int lastx)
  {
    ifstream in(s2Name.c_str());
    Spline s2(in, s2Name /* name of the spline*/);
    return greatestRelativeDifference(s1, s2, verbose, lastx);
  }
}

namespace SplineTools
{


  struct SaveSecond
  {
    void operator()(const std::pair<std::string, Spline*>& p) {
      p.second->save();
      p.second->dump();
    }
  };


  ToSplineSaver::~ToSplineSaver()
  {
    using namespace std;
    for_each(mSplineMap.begin(), mSplineMap.end(), SaveSecond());
  }

  void ToSplineSaver::saveToSpline(const double x, const double y, const std::string& name)
  {
    if (mSplineMap.find(name) == mSplineMap.end()) {
      mSplineMap[name] = new Spline(100, name, 0 /* no anchor, on purpose*/);
    }
    mSplineMap[name]->set(x, y);
  }

  struct DeleteSecond
  {
    void operator()(const std::pair<std::string, Spline*>& p) {
      delete p.second;
    }
  };

  void ToSplineSaver::reset()
  {
    using namespace std;
    for_each(mSplineMap.begin(), mSplineMap.end(), DeleteSecond());
    mSplineMap.clear();
  }



  double fromSplineOnDisk(const double x, const std::string& filename, YAccessMethod m)
  {
    static std::map<std::string, Spline*> splines;
    if (splines.find(filename) == splines.end()) {
      ifstream in(filename.c_str());
      if (!in)
        throw Bad_Error("Splinetools::fromSplineOnDisk() - cannot open " + filename);
      splines[filename] = new Spline(in, filename,
                                      0 /* no anchor, to avoid order-of-deletion problems on exit &fromSplineOnDiskAnchor*/);
      splines[filename]->arm();
    }
  switch (m) {
      case FastY: return splines[filename]->fastY(x);
      case SafeY: return splines[filename]->safeY(x);
      default: throw Bad_Error("SplineTools::fromSplineOnDisk - unknown access method");
    }
  }
}


Spline* tau2a(Cosmos& c, const Spline& s)
{
  Spline* ret = new Spline(s.size(), "createdBytau2a from " + s.name);
  for (int i = s.first(); i <= s.last(); ++i)
  {
    const double tau = s.x(i);
    const double a = c.tau2a(tau);

    ret->set(a, s.y(i));
  }
  ret->arm();
  return ret;
}

Spline* tau2lna(Cosmos& c, const Spline& s)
{
  Spline* ret = new Spline(s.size(), "createdBytau2lna from " + s.name);
  for (int i = s.first(); i <= s.last(); ++i)
  {
    const double tau = s.x(i);
    const double a = c.tau2a(tau);

    ret->set(log(a), s.y(i));
  }
  ret->arm();
  return ret;
}


Spline* tau2z(Cosmos& c, const Spline& s)
{
  Spline* ret = new Spline(s.size(), "createdBytau2z from " + s.name);
  for (int i = s.first(); i <= s.last(); ++i)
  {
    const double tau = s.x(i);
    const double z = c.tau2z(tau);
    
    ret->set(z, s.y(i));
  }
  ret->flip();
  ret->arm();
  return ret;
}

Spline* tau2zp1(Cosmos& c, const Spline& s)
{
  Spline* ret = new Spline(s.size(), "createdBytau2zp1 from " + s.name);
  for (int i = s.first(); i <= s.last(); ++i)
  {
    const double tau = s.x(i);
    const double z = c.tau2z(tau);
    
    ret->set(z+1, s.y(i));
  }
  ret->flip();
  ret->arm();
  return ret;
}

void convertFromTau2zp1(Cosmos& c, const std::vector<Spline*>& v)
{
  std::vector<Spline*>::const_iterator it, end;
  it = v.begin();
  end = v.end();
  for ( ; it != end; ++it) {
    Spline* s = *it;
    for (int i = s->first(); i <= s->last(); ++i)
    {
      const double tau = s->x(i);
      const double z = c.tau2z(tau);
      s->setX(i, z+1);
    }
    s->flip();
  }
}

void dumpQuintOmegas(QuintCosmos& c, std::string prefix="")
{
  int points = 500;
  Spline O_q( points, "Omega_q(a)" );
  for (double lna=0.; lna > log(c.a_min()) ; lna -= (1./points)) {
    double a = exp(lna);
    double tau = c.a2tau(a);
    O_q.set(a, c.Tau2Omega_q->fastY(tau));
  }

  O_q.dump(prefix+"Omega_phi");
}

void dumpOmegas(Cosmos& c, std::string prefix)
{
  try {
    QuintCosmos& qc = dynamic_cast<QuintCosmos& >(c);
    dumpQuintOmegas(qc, prefix);
  } catch(std::bad_cast & b) {
    ;
  }
  int points = 500;
  Spline O_m( points, "Omega_m(a)" );
  Spline O_cdm( points, "Omega_cdm(a)" );
  Spline O_b( points, "Omega_b(a)" );
  Spline O_v( points, "Omega_v(a)" );
  Spline O_r( points, "Omega_g(a)" );
  Spline O_n( points, "Omega_nu_massless(a)" );
  Spline O_nm( points, "Omega_nu_massive(a)" );
  Spline H_a(points, "H(a)/H_0");

  for (double lna=0.; lna > log(c.a_min()) ; lna -= (1./points)) {
    double a = exp(lna);
    double tau = c.a2tau(a);
    O_m.set(a, c.tau2rho_m(tau)/c.tau2rho(tau));
    O_b.set(a, c.tau2rho_b(tau)/c.tau2rho(tau));
    O_cdm.set(a, c.tau2rho_cdm(tau)/c.tau2rho(tau));
    O_r.set(a, c.tau2rho_g(tau)/c.tau2rho(tau));
    O_n.set(a, c.tau2rho_nu(tau)/c.tau2rho(tau));
    O_nm.set(a, c.tau2rho_nuNR(tau)/c.tau2rho(tau));
    O_v.set(a, c.tau2rho_v()/c.tau2rho(tau));
    H_a.set(a, );
  }

  O_m.dump("Omega_m");
  O_b.dump("Omega_b");
  O_cdm.dump("Omega_cdm");
  O_v.dump("Omega_v");
  O_r.dump("Omega_gamma");
  O_nm.dump("Omega_massive_nu");
  O_n.dump("Omega_nu");
}


static void dumpRhoSpline(Cosmos& c, Spline* spl, std::string name)
{
  if (spl->size() <= 1 )
    return;
  Spline s;
  for (int i = spl->first(); i <= spl->last(); ++i) {
    s.set(spl->x(i), spl->y(i));
  }
  s *= 1./c.M_p(2);
  tau2zp1(c, s)->dump(name);
}

static void dumpQSpline(QuintCosmos& c, Spline* spl, std::string name)
{
  if (spl->size() <= 1 )
    return;
  Spline s;
  for (int i = spl->first(); i <= spl->last(); ++i) {
    s.set(spl->x(i), spl->y(i));
  }
  s *= 1./c.M_p(1);
  tau2zp1(c, s)->dump(name);
}

static void dumpQuintDensities(QuintCosmos& c, std::string prefix = "" )
{
  dumpRhoSpline(c, c.Tau2Rho_qpot, prefix+"zp1-rho_qpot");
  dumpRhoSpline(c, c.Tau2Rho_qkin, prefix+"zp1-rho_qkin");
  dumpRhoSpline(c, c.Tau2Rho_q, prefix+"zp1-rho_q");
  dumpRhoSpline(c, c.Tau2P_q, prefix+"zp1-pres_q");
  dumpRhoSpline(c, c.Tau2Vprime, prefix+"zp1-Vprime");
  tau2zp1(c, *(c.Tau2W_q))->dump(prefix+"zp1-wq");
}

void dumpDensities(Cosmos& c, std::string prefix)
{
  try {
    QuintCosmos& qc = dynamic_cast<QuintCosmos& >(c);
    dumpQuintDensities(qc, prefix);
  } catch(std::bad_cast & b) {
    ;
  }
  dumpRhoSpline(c, c.Tau2Rho, prefix+"zp1-rho");
  dumpRhoSpline(c, c.Tau2Rho_g, prefix+"zp1-rho_g");
  dumpRhoSpline(c, c.Tau2Rho_b, prefix+"zp1-rho_b");
  dumpRhoSpline(c, c.Tau2Rho_nu, prefix+"zp1-rho_nu");
  dumpRhoSpline(c, c.Tau2Rho_nuNR, prefix+"zp1-rho_nuNR");
  dumpRhoSpline(c, c.Tau2P_nuNR, prefix+"zp1-p_nuNR");
  dumpRhoSpline(c, c.Tau2Rho_cdm, prefix+"zp1-rho_cdm");
  dumpRhoSpline(c, c.Tau2Rho_m, prefix+"zp1-rho_m");
  //tau2zp1(c, *(c.Tau2w_nuNR))->dump(prefix+"zp1-w_nuNR");
  dumpRhoSpline(c, c.Tau2ADotToA, prefix+"Hubble_a");
  dumpRhoSpline(c, c.Tau2P,prefix+"zp1-pres");
}

