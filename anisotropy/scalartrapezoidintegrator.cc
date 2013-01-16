#include "scalartrapezoidintegrator.h"
#include "spline.h"
#include "cosmos.h"
#include "controlpanel.h"

ScalarTrapezoidIntegrator::ScalarTrapezoidIntegrator(Cosmos *c, double mi, double ma, SplineWeb* d, SplineWeb *dp) : ScalarIntegrator(c,mi,ma,d,dp) {}



ScalarTrapezoidIntegrator::~ScalarTrapezoidIntegrator() {}

/*!
  The inner most core of cmbeasy.
  Performs the convolution of bessel function with source and the
  final k - integration. Use ofstream monitor to take a look at the integrand,
  if you like. 

  Please note, that we stop integration, if after 200 steps, the integral
  doesn't increase significantly. In order to increase precision, you may
  either want to play around with when this happens, or disable it
  altogether. The code for it is in the last few lines.

  One last word on precision: You may think that the precision in the 
  convolution translates directly over to the overall precision. Yet this is not
  true. Just take a look at monitor (the file, as mentioned above): You may introduce
  some "noise", yet for as long as it is not too large, its effects will largely cancel,
  as one has the noise for several hundred k-values.

*/
ClReturn ScalarTrapezoidIntegrator::integrate(int l, Spline *bessel, const ControlPanel &control) {
  //  ofstream monitor("monitor.dat");

  //cout << "Integrating l = " << l << endl;

  ClReturn clReturn;
  clReturn.null(cosmos->InitialPower.size());
  double tau0 = cosmos->tau_0();
 
  vector<double> accu(nko+1); // here we note the integration contribution to  clts to determine end of k integration
  int accuCount =-1;   // counts the accu-vector entry number
  unsigned int mp = (cosmos->InitialPower.size()-1) >> 1;    // middle index of all power-spectras
  int c=0; 
  // small loop to start at high enough k, where bessel function and
  // sources have at least some overlap, otherwise we would catch a lot of noConvolutionOverlaps
  while (K[c] * tau0 < bessel->start() && c <= nko) c++;

  //  int cnt=0,cnt2=0;

  while (c < nko) { 
    double k = K[c];
    double d=0, dp=0;  // the convolution results

    double tau2 = source[c]->start();   // we start from here. note that tau1 = tau2 = tau_start at the beginning  in the loop 
    double tau1=0;  // init
    double xf2 =  k*(tau0 - tau2);  //argument of bessel function 
    double xf1=0; //init 
    double b2 = 0,b1=0; // Bessel function value at pos "1" and "2". Strictly speaking not correct for init,but we take it 
    int sz = source[c]->size();
    double *tt = source[c]->ptrY(0);
    double *pp = sourceP[c]->ptrY(0);
    double *timesteps = source[0]->xdat;
    int j = 0;
    for (int i = 1; i < sz; i++,j++) {
      // first, we init the left point variables "1" with the old right variables 
      tau1 = tau2; 
      xf1 = xf2;
      b1 = b2;	 

      //      if (tau1 > 1000 && l > 400) break;  speeds it up even more
      // now get the updated values for the right "2" variables
      tau2 = timesteps[i];

      xf2 = k*(tau0 - tau2);
      //
      // Catch l=2 case and use explicit formula for bessel function, as tabulated values
      //  are bad behaved
      //
      if (l == 2) {
          if (xf2 < 1e-1) {  // for such small argument, use taylor expansion, better behaved
              b2 = xf2*xf2*(1.0/15.0 + xf2*xf2*(-1.0/210.0 + xf2*xf2/7560.0));
          } else b2 = (3.0*(sin(xf2)/xf2-cos(xf2))/xf2- sin(xf2))/xf2;
      } else {  // l not equal to 2 
          if (xf2 < bessel->start()) b2=0;  // bessel function is 0 in that region
          else b2 = bessel->fastY(xf2);  // ask the bessel spline 
      }
      if (b1 || b2) {  // if either of the two bessel functions is != 0
          double t = tt[i]*b2+tt[j]*b1;
          double p = pp[i]*b2+pp[j]*b1;	
          double tau_step = tau2-tau1;
          d += t*tau_step;
          dp += p*tau_step;
      }
    }
    d *= 0.5;  // the integration was a trapezoidal rule and we have to correct for the factor 2 
    dp *= 0.5;

    for (unsigned int q = 0; q < cosmos->InitialPower.size(); q++) {
      double power = cosmos->scalarSpectrum(k, q)/k;
      clReturn.clts[q] += dK[c] * d*d*power; 
      clReturn.cles[q] += dK[c] * dp*dp*power;
      clReturn.clcs[q] += dK[c] * d*dp*power;
      if (q == mp) accu[++accuCount] = dK[c] * d*d*power;   // note current integration contribution
      //      monitor << k << "     " << dK[c] * d*d*power << "   " <<  clReturn.clts[mp] << "  "  << d <<  endl;
    }

    // if we have at least more than 200 integration points, we check,
    // if the last ten points together make up less than 1e-4 of the total
    // result up to now. If so, we stop at this k value and return;
    double sum=0;
    if (accuCount > 200) {
      for (int i = accuCount; i > accuCount-10; i--) sum+= accu[i];
      if ((sum  / clReturn.clts[mp]) < 1e-4) {  // 1e-4
          // cout << "break at: " << rc << "  k: " << k << "   after # steps: " << accuCount << endl;
          break;
      }
    }
    c++;
  }
  //  cout << "cnt: " << cnt << " cnt2: " << cnt2 << endl;
  //  if (l ==1000) throw Bad_Error("stop");
  return clReturn;
}
