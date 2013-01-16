#ifndef WMAP_7YR_TT_PIXLIKE_H
#define WMAP_7YR_TT_PIXLIKE_H
/****************************************************************
 *   This file is part of the C++ port of the Fortran likelihood
 *   code for the WMAP 5yr release provided by the WMAP team
 *   at http://lambda.gsfc.nasa.gov/ .
 *
 *   The code was ported by Georg Robbers for easier
 *   interfacing with cmbeasy (http://www.cmbeasy.org).
 *   Bugs in this port should be reported to the
 *   cmbeasy authors (bugs@cmbeasy.org).
 ****************************************************************/

#include "port_helper.h"

namespace  wmap_tlike
{
  void  setup_for_tt_exact(int nlhigh, int* tt_ngood = 0);
  void  compute_tt_pixlike(int nlhigh,double_1d& cl_in, double& chisq, double& lndet);
  unsigned int tt_pixlike_dof();
}

#endif // WMAP_7YR_TT_PIXLIKE_H

