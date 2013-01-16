#ifndef WMAP_5YR_TEEEBB_PIXLIKE_H
#define WMAP_5YR_TEEEBB_PIXLIKE_H
/****************************************************************
 *   This file is part of the C++ port of the Fortran likelihood
 *   code for the WMAP 7yr release provided by the WMAP team
 *   at http://lambda.gsfc.nasa.gov/ .
 *
 *   The code was ported by Georg Robbers for easier
 *   interfacing with cmbeasy (http://www.cmbeasy.org).
 *   Bugs in this port should be reported to the
 *   cmbeasy authors (bugs@cmbeasy.org).
 ****************************************************************/

#include "port_helper.h"

#include "gsl/gsl_matrix.h"

namespace  wmap_teeebb_lowl
{
  void cholesky_invert (gsl_matrix* A);
  void  teeebb_lowl_like_setup();
  void  teeebb_lowl_likelihood(int nlmaxin, real8_1d& Clttin, real8_1d& Cltein, real8_1d& Cleein, real8_1d& Clbbin, REAL_8& chisq_r3, REAL_8& lndet);
  unsigned int teeebb_pixlike_dof();
}

#endif // wmap_5YR_TEEEBB_PIXLIKE_H

