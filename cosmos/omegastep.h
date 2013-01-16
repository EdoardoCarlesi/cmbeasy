#ifndef OMEGASTEP_H
#define OMEGASTEP_H

#include "arbitraryomega.h"

/*!
  Omegastep class, the model of Eq. 6 of astro-ph/0601544

  Omega_e and w_0 are the parameters for the standard form,
  param[2] is gamma from the appendix

*/

class OmegaStep : public ArbitraryOmega {
 public:
  OmegaStep(QuintCosmos &c) : ArbitraryOmega(c) {
    param.resize(3);
    param[0] = -0.9; // w0
    param[1] = 0.1; // Omega_e
    param[2] = 1; // gamma
  }

  virtual Type type() { return omegastep;}
  virtual void prepare();
};


#endif
