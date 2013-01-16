#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "safevector.h"
#include "anchor.h"

class Cosmos;
class Spline;
class SplineWeb;
class ControlPanel;

/*!
  Small struct to pass the integration
  results, i.e. the multipole values for
  one multipole from the integrators
  back to CmbCalc.
  
  The values come as vectors, because
  there is a result for each spectral
  index.
*/
struct ClReturn {
  SafeVector<double> clts, cles, clcs;
  SafeVector<double> cltt,clet,clbt,clct;
  
  /*! resize all vectors */
  void resize(unsigned int n) {
    clts.resize(n); cles.resize(n); clcs.resize(n);
    cltt.resize(n); clet.resize(n); clct.resize(n);
    clbt.resize(n);
  }
  /*! resize to n and fill all vectors with zeros */
  void null(unsigned int n) {
    resize(n);
    for (unsigned int q = 0; q < n; q++) { 
      clts[q] = 0; cles[q] = 0; clcs[q] = 0;
      cltt[q] = 0; clet[q] = 0; clct[q] = 0;
      clbt[q] = 0;
    }
  }
  
};


/*!
  Base class for all Integrators. Folds the sources with the bessel functions.

  It also provides the values of K and the interval dK at which the sources
  should be sampled. These k-values are used by ScalarIntegrator and
  TensorIntegrator. 
*/
class Integrator {
 protected:
  Anchor anchor;
  Cosmos *cosmos;
  double akmin,akmax;
  int nko;
  vector<double> K, dK;  //!< k values and interval for integration

 public:

  virtual void clean_up() {}; //!< should clean up be necessary, we define it here. Up to now, it isn't. Just for future use.
  Integrator(Cosmos *c, double mi, double ma);

  virtual ClReturn integrate(int l, Spline *bessel, const ControlPanel&)=0; //!< integrate for multipole value l and corresponding bessel function the perturbations. Results are stored in ClReturn    
  virtual ~Integrator();
};

#endif
