#include "omegastep.h"
#include "quintcosmos.h"
#include "spline.h"

void OmegaStep::prepare()
{
  anchor.kill();
  W = new Spline(1000,"Omegastep_w",&anchor);
  DW = new Spline(W,"Omegastep_dw",&anchor);
  DDW = new Spline(W,"Omegastep_ddw",&anchor);

  Omega = new Spline(1000,"Omegastep_omega",&anchor);
  DOmega = new Spline(Omega,"Omegastep_domega",&anchor);


  //
  //  init model parameters
  //

  double w0 = param[0];
  double Oe = param[1];
  double gamma = param[2];
  double Om = cosmos.omega_m()+cosmos.omega_nuNR();
  double Oq = cosmos.omega_q();

 double y = -0.3;
 double y_step = 0.001;

 // Fill the Omega spline
 do {
   double fade = (1-exp(3*w0*y));
   if (fade<=0) fade = 0;
   fade = pow(fade,gamma);

   double omega = ( Oq - fade * Oe ) / (Oq + Om*exp(-3*w0*y)) + Oe*fade; 

   Omega->set(y,omega);
   if (fade > 1e-2 && fade < 0.98) y+= y_step; else y+= 10*y_step;  // make stepping small in "knee"
 } while (y < 50);
 Omega->arm();
 Omega->derive(*DOmega,Spline::yes);

 y = -0.25;

 // compute w(z) from Omega_de(z)
 do {

   double a = exp(-y);
   double rho_m = cosmos.initialRho_b(a)+cosmos.initialRho_cdm(a);
   double rho_rel = cosmos.initialRho_g(a)+cosmos.initialRho_nu(a);
   double rho_nuNR, p_nuNR;
   rho_nuNR = cosmos.initialRho_nuNR(a, &p_nuNR);

   double Oq = Omega->fastY(y);
   double w = DOmega->fastY(y) / (3. * Oq*(1.-Oq) ) +  (1./3.)*(rho_rel+3.*p_nuNR)/(rho_m+rho_rel+rho_nuNR);

   if (y < 0) w = w0;

   W->set(y,w);
   y += y_step;

 } while (y < 50);

 W->arm();
 W->derive(*DW,Spline::yes);
 DW->derive(*DDW,Spline::yes);

 Arbitrary::prepare();  // Arbitrary, not ArbitraryOmega (needs to be cleaned up)
}
