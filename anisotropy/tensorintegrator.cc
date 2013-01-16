#include "tensorintegrator.h"
#include "cosmos.h"
#include "controlpanel.h"

/*!
  Create TensorIntegrator. 
  This sets basic variables, like the source webs Dt etc. 

  Then, it calculates the k values at which the integrand of the 
  final C_l calculation will be evaluated. For each k value, the sources
  are interpolated from the webs Dt etc. The sources are functions
  of conformal time tau.
*/

TensorIntegrator::TensorIntegrator(Cosmos *c,double mi, double ma,double stop,SplineWeb *dt, SplineWeb*  dte, SplineWeb* dtb) :  Integrator(c,mi,ma),  Dt(dt), Dte(dte), Dtb(dtb) , stpt(stop) { 
  source.resize(nko);
  sourceE.resize(nko); 
  sourceB.resize(nko); 
    
  for (int c =0; c < nko; c++) {
    if (c==0) source[c] = Dt->createAlongY(K[c], "tensorIntegrator_s", &anchor); 
    else source[c] = Dt->createAlongY(K[c], "tensorIntegrator_s", &anchor,source[0]); 
    sourceE[c] = Dte->createAlongY(K[c], "tensorIntegrator_s", &anchor,source[0]);
    sourceB[c] = Dtb->createAlongY(K[c], "tensorIntegrator_s", &anchor,source[0]); 
  }
  source[0]->arm(Spline::all);
}

/*!
  Perform an integration for one given value of l and the bessel function j_l.

  Basically, it loops though all k-values determined by the constructor and
  for each of these convolutes the sources with bessel functions j_l.

  These results are integrated over and returned.

  Use ofstream monitor below, to explicitly look at the integrand
  of the k-integral, if you like.

*/
ClReturn TensorIntegrator::integrate(int l, Spline *bessel, const ControlPanel &control) {
  //ofstream monitor("monitor.dat");
  //cout << "+"; cout.flush();

  ClReturn clReturn;
  clReturn.null(cosmos->InitialPower.size());

  double tau0 = cosmos->tau_0();
  vector<double> accu(nko+1); // here we note the integration contribution to  clts to determine end of k integration
  int accuCount =-1;   // counts the accu-vector entry number
  unsigned int mp = (cosmos->InitialTensorPower.size()-1) >> 1;    // middle index of all power-spectras
  int c=0; // c four count :-)
  // small loop to start at high enaugh k, where at least, bessel function and
  // sources have some overlap, otherwise we would catch a lot of noConvolutionOverlaps
  while (K[c] * tau0 < bessel->start() && c <= nko) c++;  

  double lastk=0;
  while (c < nko) { 
    double k = K[c];  // for convenience, small k is this K[c] mode
    // cout << "actual k is : " << k << endl;
   
    // now, we convolute with the bessel function. We ask Spline::vectorConvolution to
    // do this for us
    double d=0, de=0,db=0;
   	  
    if (Spline::convolutionVektor.size() != 3) Spline::convolutionVektor.resize(3);
    Spline::convolutionVektor[0] = source[c];
    Spline::convolutionVektor[1] = sourceE[c];
    Spline::convolutionVektor[2] = sourceB[c];
    
    double precision =2e-4; // convolution precision
    
    try {
      double stopTau = stpt/k;   // we can stop cause from k*tau > stpt on, the sources are 0
      stopTau = min(tau0,stopTau);
      Spline::vectorConvolution2(bessel, 0,stopTau, precision,-k,k*tau0); 
      d = source[c]->convolutionResult();
      de = sourceE[c]->convolutionResult(); 
      db = sourceB[c]->convolutionResult();
    } catch (noConvolutionOverlap) {}


    if (l > 1600 && lastk / k < 0.8) {
      cout << "dump at l : " << l << " k: " << k << endl;
      source[c]->dump("source",false);
      int nxt;
      cin >> nxt;
      lastk=k;
    }
   

    // ======================
    // CATCH l=2 and use explicit formula
    // ======================

    if (l==2) {
      d = 0;
      double tau_min=source[c]->start();
      double tau_stop=source[c]->stop();
      double tau_step = (tau_stop-tau_min)*1e-3;

      for (double tau=tau_min; tau <= tau_stop; tau+=tau_step) {
	double xf = k*(tau0 - tau);
	double b = 0;
	if (xf < 1e-1) {  // for such small argument, use taylor expansion, better behaved
	  b = xf*xf*(1.0/15.0 + xf*xf*(-1.0/210.0 + xf*xf/7560.0));
	} else { 
	  b = (3.0*(sin(xf)/xf-cos(xf))/xf- sin(xf))/xf;
	}
	d += source[c]->fastY(tau)*b * tau_step;
      } 
    }
    // ======================
   

    for (unsigned int q = 0; q < cosmos->InitialPower.size(); q++) {
      double power = cosmos->tensorSpectrum(k, q, control)/k; 
      if (power == 0) throw Bad_Error("tensorSpectrum() returned 0.\nMost propably, this is due to an scalar spectral index of 1, leading to n_t =0.\nIt may also be that you specified explicitly a tensor spectral index of 0");
      clReturn.cltt[q] += dK[c] * d*d*power;
      clReturn.clet[q] += dK[c] * de*de*power;
      clReturn.clbt[q] += dK[c] * db*db*power;
      clReturn.clct[q] += dK[c] * d*de*power; 
      //monitor << k << "  " << power * d*d * dK[c] <<  "  " << clReturn.cltt[q] << "  " << d << " "  << power << endl;
      if (q == mp) accu[++accuCount] =  dK[c] * d*d*power;
    }
    
    // if we have at least more than 200 integration points, we check,
    // if the last ten points together make up less than 1e-4 of the total
    // result up to now. If so, we stop at this k value and return;
    double sum=0;
    if (accuCount > 200) {
      for (int i = accuCount; i > accuCount-10; i--) sum+= accu[i];
      if ((sum  / clReturn.cltt[mp]) < 1e-4) { // 1e-4
	//cout << "break at: " << rc << "  k: " << k << "   after # steps: " << accuCount << endl;
	break;
      }	
    }

    //monitor.flush();
    c++;
  }
  return clReturn;
}
