#ifndef WMAP_7YR_LIKELIHOOD_H
#define WMAP_7YR_LIKELIHOOD_H
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

// This code is the central likelihood routine from which the other subroutines
// are called.
//
// Parameters are defined in associated WMAP_options module
//
// This code was a collaborative effort by the following WMAP Team Members:
// R. Bean
// O. Dore
// J. Dunkley
// E. Komatsu
// D. Larson
// M. Nolta
// H. Peiris
// L. Verde
//
// ===========================================================================

#include "port_helper.h"

#include "WMAP_7yr_util.h"
#include "WMAP_7yr_tt_beam_and_ptsrc_chisq.h"
#include "WMAP_7yr_tt_pixlike.h"

#ifdef USE_LOWELL_TBEB
#include "WMAP_7yr_tetbeebbeb_pixlike.h"
#else
#include "WMAP_7yr_teeebb_pixlike.h"
#endif


namespace wmap_likelihood_7yr
{
  void wmap_likelihood_init();
  void wmap_likelihood_dof(int& tt_npix, int& teeebb_npix);
#ifdef USE_LOWELL_TBEB
void wmap_likelihood_compute(real8_1d& cltt, real8_1d& clte, real8_1d& cltb, real8_1d&clee,
                             real8_1d& cleb, real8_1d& clbb, Array1D<REAL_8>& like);
#elif USE_HIGHELL_TB
void wmap_likelihood_compute(real8_1d& cltt, real8_1d& clte, real8_1d& cltb, real8_1d& clee,
                             real8_1d& cleb, real8_1d& clbb, Array1D<REAL_8>& like);
#else
void wmap_likelihood_compute(real8_1d& cltt, real8_1d& clte, real8_1d& clee, real8_1d& clbb,
                             Array1D<REAL_8>& like);
#endif
}

#endif // WMAP_7YR_LIKELIHOOD_H
