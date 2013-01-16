#include "tensortrapezoidintegrator.h"
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

TensorTrapezoidIntegrator::TensorTrapezoidIntegrator(Cosmos *c,double mi, double ma,double stop,SplineWeb *dt, SplineWeb*  dte, SplineWeb* dtb) :  TensorIntegrator(c,mi,ma,stop,dt,dte,dtb) {}

/*!
  Perform an integration for one given value of l and the bessel function j_l.

  Basically, it loops though all k-values determined by the constructor and
  for each of these convolutes the sources with bessel functions j_l.

  These results are integrated over and returned.

  Use ofstream monitor below, to explicitly look at the integrand
  of the k-integral, if you like.

*/
ClReturn TensorTrapezoidIntegrator::integrate(int l, Spline *bessel, const ControlPanel &control) {
  //ofstream monitor("tensmonitor.dat");
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
 
  while (c < nko) {  
    double k = K[c];  // for convenience, small k is this K[c] mode
    double d=0, de=0,db=0;

   double tau2 = source[c]->start();   // we start from here. note that tau1 = tau2 = tau_start at the beginning  in the loop 
    double tau1=0;  // init
    double xf2 =  k*(tau0 - tau2);  //argument of bessel function 
    double xf1=0; //init 
    double b2 = 0,b1=0; // Bessel function value at pos "1" and "2". Strictly speaking not correct for init,but we take it 
    int sz = source[c]->size();
    double *tt = source[c]->ptrY(0);
    double *ee = sourceE[c]->ptrY(0);
    double *bb = sourceB[c]->ptrY(0);
    double *timesteps = source[0]->xdat;
    int j = 0;
    for (int i = 1; i < sz; i++,j++) {
      // first, we init the left point variables "1" with the old right variables 
      tau1 = tau2; 
      xf1 = xf2;
      b1 = b2;	 

      // now get the updated values for the right "2" variables
      tau2 = timesteps[i];      
      if (tau2 >= stpt/k) { 
	//cout << "stoppijng for k: " << k << " at tau: " << tau2 << endl;
	break; // stop at the stop tensor, afterwards sources would be 0 
      }
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
	double e = ee[i]*b2+ee[j]*b1;	
	double b = bb[i]*b2 +bb[j]*b1;
	double tau_step = tau2-tau1;

	d += t*tau_step;
	de += e*tau_step;
	db += b*tau_step;
      }
    } 
    d *= 0.5;  // the integration was a trapezoidal rule and we have to correct for the factor 2 
    de *= 0.5;
    db *= 0.5;
   	  
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
    c++;
  }
  //monitor.flush();
  //if (l > 1600) throw Bad_Error("stop");
  return clReturn;
}
