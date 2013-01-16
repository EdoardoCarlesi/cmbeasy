#include "exponentialcoupling.h"

#include "coupledleaping.h"
#include "spline.h"

#include <iostream>

ExponentialCoupling::ExponentialCoupling()
  : mb(0), mbetaSpline(0), mbetaPrimeSpline(0), mbetaDoublePrimeSpline(0)
{
}

ExponentialCoupling::~ExponentialCoupling()
{
  delete mbetaSpline;
  mbetaSpline = 0;

  delete mbetaPrimeSpline;
  mbetaPrimeSpline = 0;

  delete mbetaDoublePrimeSpline;
  mbetaDoublePrimeSpline=0;
}

void ExponentialCoupling::setCouplingParameters(const vector<double>&param)
{
  if (param.size() == 1)
    mb=param[0];
  else
    throw Bad_Error("Error in ExponentialCoupling::setCouplingParameters: \
                     Expected parameter vector of size 1");
}

void ExponentialCoupling::setCouplingParameters(double beta){
	mb=beta;
	cout<< " Exp.Coupling set to:  " << beta;
}

vector<double> ExponentialCoupling::getParameters() const
{
  vector<double> vec(1);
  vec[0]=mb;
  return vec;
}

void ExponentialCoupling::postPrepare()
{
  delete mbetaSpline;
  mbetaSpline = 0;

  delete mbetaPrimeSpline;
  mbetaPrimeSpline = 0;

  delete mbetaDoublePrimeSpline;
  mbetaDoublePrimeSpline=0;

/*
  // Special case Leaping Kinetic Term quintessence:
  // if our quintessence model is LKT, we have to set the effective
  // coupling beta after the potential has been rewritten to the standard form.

  CoupledLeaping* cl = dynamic_cast<CoupledLeaping*>(mQuint);

  if (cl)
  {
    mbetaSpline = new Spline(cl->vSpline()->size(), "ExponentialCoupling::betaSpline");
    int last = cl->vSpline()->last();
    for (int i = cl->vSpline()->first(); i <= last; ++i)
    {
      const double phi = cl->vSpline()->x(i)*baseCosmos::M_p();
      mbetaSpline->set(phi, phi2beta(phi));
    }

    mbetaSpline->arm();
    mbetaPrimeSpline = new Spline(mbetaSpline, "ExponentialCoupling::betaPrimeSpline");
    mbetaDoublePrimeSpline = new Spline(mbetaSpline, "ExponentialCoupling::betaPrimeDoubleSpline");

    mbetaSpline->derive(*mbetaPrimeSpline);
    mbetaPrimeSpline->arm();
    mbetaPrimeSpline->derive(*mbetaDoublePrimeSpline);
    mbetaDoublePrimeSpline->arm();

    //mbetaSpline->dump("testbeta");
    //mbetaPrimeSpline->dump("testbetaprime");
    //mbetaDoublePrimeSpline->dump("testbetadoubleprime");
  }
*/
}

