#ifndef WMAP_7YR_UTIL_H
#define WMAP_7YR_UTIL_H
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

#include <string>

namespace  wmap_util
{
  bool  wmap_likelihood_ok();
  void  wmap_likelihood_error_init( );
  void  wmap_likelihood_error( std::string msg, int code );
  void  wmap_likelihood_warning( std::string msg, int code );
  void  wmap_likelihood_error_report();
  void  wmap_timing_start( std::string msg );
  void  wmap_timing_checkpoint( std::string msg );
  void  wmap_timing_end();
}

#endif // WMAP_7YR_UTIL_H

