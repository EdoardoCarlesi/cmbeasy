#include "scalarintegrator.h"
#include "spline.h"
#include "cosmos.h"
#include "controlpanel.h"

#ifdef _OPENMP
#include "omp.h"
#endif

ScalarIntegrator::ScalarIntegrator(Cosmos *c, double mi, double ma, SplineWeb* d, SplineWeb *dp) : Integrator(c,mi,ma), D(d), Dp(dp) {

  source.resize(nko+1);
  sourceP.resize(nko+1);

  source[0] = D->createAlongY(K[0], "scalarIntegrator", &anchor);
  sourceP[0] = Dp->createAlongY(K[0], "scalarIntegrator_p", &anchor,source[0]);
  //double time_before = omp_get_wtime();
#pragma omp parallel for
  for (int c=1; c<nko; c++) {
    source[c] = D->createAlongY(K[c], "scalarIntegrator", &anchor,source[0]);
    sourceP[c] = Dp->createAlongY(K[c], "scalarIntegrator_p", &anchor,source[0]);
  }
  //double time_after = omp_get_wtime();
  //cout << "slicing sources " << (time_after-time_before) << endl;

  //time_before = omp_get_wtime();
  source[0]->arm(Spline::all);
  //time_after = omp_get_wtime();
  //cout << "arming sources " << (time_after-time_before) << endl;
}



ScalarIntegrator::~ScalarIntegrator() {}

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
ClReturn ScalarIntegrator::integrate(int l, Spline *bessel, const ControlPanel &control) {
  // ofstream monitor("monitor.dat");

  ClReturn clReturn;
  clReturn.null(cosmos->InitialPower.size());
  double tau0 = cosmos->tau_0();

  vector<double> accu(nko+1); // here we note the integration contribution to  clts to determine end of k integration
  int accuCount =-1;   // counts the accu-vector entry number
  unsigned int mp = (cosmos->InitialPower.size()-1) >> 1;    // middle index of all power-spectras
  int c=0; 
  // small loop to start at high enaugh k, where at least, bessel function and
  // sources have some overlap, otherwise we would catch a lot of noConvolutionOverlaps
  while (K[c] * tau0 < bessel->start() && c <= nko) c++;  

  while (c < nko) { 
    double k = K[c];
    double d=0, dp=0;  // the convolution results
  
    if (Spline::convolutionVektor.size() != 2) Spline::convolutionVektor.resize(2);
    Spline::convolutionVektor[0] = source[c];
    Spline::convolutionVektor[1] = sourceP[c];
    
    // these are the "speedy" lines of code :-)
    // if you want to be on the real good side, keep precision = 1e-4 or smaller
    // and never set it to higher values 
    // i.e. comment the if (l > 100) precision = ... statetements
    double precision = 1e-4;
    if (l > 100) precision = 5e-4;
    if (l > 400) precision = 1e-3;
    if (l > 1000) precision = 1e-3;

    if (l != 2) { // all except l = 2 are calculated using pre-calculated bessel functions
      try {
	double stopTau=tau0;
	Spline::vectorConvolution2(bessel, 0,stopTau, precision,-k,k*tau0); // this is the inner most loop
	d = source[c]->convolutionResult();
	dp = sourceP[c]->convolutionResult();    
      } catch (noConvolutionOverlap) {}   // quite rarely, cause we start from high enough k  
    } else {
      // ======================
      // CATCH l=2 and use explicit formula
      // in this case. Bessel function for l=2 not sooo good from
      // table on disk
      // ======================
      double tau_min=source[c]->start();
      double tau_stop=source[c]->stop();
      double tau_step = (tau_stop-tau_min)*1e-4;
      
      for (double tau=tau_min; tau <= tau_stop; tau+=tau_step) {
	double xf = k*(tau0 - tau);
	double b = 0;
	if (xf < 1e-1) {  // for such small argument, use taylor expansion, better behaved
	  b = xf*xf*(1.0/15.0 + xf*xf*(-1.0/210.0 + xf*xf/7560.0));
	} else { 
	  b = (3.0*(sin(xf)/xf-cos(xf))/xf- sin(xf))/xf;
        }
        d += source[c]->fastY(tau)*b * tau_step;
        dp += sourceP[c]->fastY(tau)*b * tau_step;
      } 
    }

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
  //  if (l ==2) throw Bad_Error("stop");
  return clReturn;
}
