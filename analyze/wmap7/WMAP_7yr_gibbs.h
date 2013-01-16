#ifndef WMAP_7YR_GIBBS_H
#define WMAP_7YR_GIBBS_H
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

namespace wmap_gibbs
{
  void  setup_for_tt_gibbs();
  void  compute_tt_gibbslike(const real8_1d& cl_in, double& like);
}

#endif // WMAP_7YR_GIBBS_H
