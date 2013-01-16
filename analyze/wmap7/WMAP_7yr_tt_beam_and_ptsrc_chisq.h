#ifndef WMAP_7YR_TT_BEAM_AND_PTSRC_CORR_H
#define WMAP_7YR_TT_BEAM_AND_PTSRC_CORR_H
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

namespace  wmap_tt_beam_ptsrc_chisq
{
  void  init_tt_beam_and_ptsrc_chisq(int lmin, int lmax);
  void  quit_tt_beam_and_ptsrc_chisq( );
  REAL_8 compute_tt_beam_and_ptsrc_chisq( int lmin, int lmax, Array1D<REAL_8> &cltt,
                                        Array1D<REAL_8> &cltt_dat, Array1D<REAL_8> &neff,
                                        Matrix<REAL_8>& fisher, Array1D<REAL_8> &z,
                                        Array1D<REAL_8> &zbar );
}

#endif // WMAP_7YR_TT_BEAM_AND_PTSRC_CORR_H

