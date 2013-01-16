#include "arbitraryomega.h"
#include "quintcosmos.h"
#include "spline.h"

void ArbitraryOmega::prepare()
{
  anchor.kill();
  W = new Spline(1000, "ArbitraryOmega_w", &anchor);
  DW = new Spline(W, "ArbitraryOmega_dw", &anchor);
  DDW = new Spline(W, "ArbitraryOmega_ddw", &anchor);

  double y = -0.25, y_step = 1e-3;
  do {
    double z = exp(y) - 1;
    double a = 1.0/(z+1.0);
    double Oq = omega(z); 
    double h = max(1e-5, z*1e-5);
    double dOq_dz = (omega(z+h) - omega(z-h) ) / (2*h);
    double dOq_dlnz = dOq_dz * (1+z);

    double rho_m = cosmos.initialRho_b(a)+cosmos.initialRho_cdm(a);
    double rho_rel = cosmos.initialRho_g(a)+cosmos.initialRho_nu(a);
    double rho_nuNR, p_nuNR;
    rho_nuNR = cosmos.initialRho_nuNR(a, &p_nuNR);

    double w = dOq_dlnz / (3. * Oq*(1.-Oq) ) +  (1./3.)*(rho_rel+3.*p_nuNR)/(rho_m+rho_rel+rho_nuNR);

    W->set(y,w);
    y += y_step;
 } while (y < 50);

  W->arm();
  W->derive(*DW,Spline::yes);
  DW->derive(*DDW,Spline::yes);

  Arbitrary::prepare();
}
