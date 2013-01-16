#ifndef BR_MOD_DIST_H
#define BR_MOD_DIST_H
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
  // *********************************************************************
  // *      br_mod -- An F90 module for computing the Blackwell-Rao      *
  // *                estimator given signal samples from the posterior  *
  // *                                                                   *
  // *                 Written by Hans Kristian Eriksen                  *
  // *                                                                   *
  // *                Copyright 2006, all rights reserved                *
  // *                                                                   *
  // *                                                                   *
  // *   NB! The code is provided as is, and *no* guarantees are given   *
  // *       as far as either accuracy or correctness goes.              *
  // *                                                                   *
  // *  If used for published results, please cite these papers:         *
  // *                                                                   *
  // *      - Eriksen et al. 2006, ApJ, submitted, astro-ph/0606088      *
  // *      - Chu et al. 2005, Phys. Rev. D, 71, 103002                  *
  // *                                                                   *
  // *********************************************************************

#include "port_helper.h"

namespace  br_mod_dist {

  void  initialize_br_mod(const int lmin_in, const double_3d& sigmas_in);
  void  clean_up_br_mod();
  void  compute_br_estimator(const double_1d& cls, double& lnL);
  void  read_gibbs_chain(const std::string& filename, long int& lmax, long int& numchains,
                                                      long int& numsamples, FloatToDoubleArray3D& data);
  void  compute_largest_term(const double_1d& cls);
}

#endif // BR_MOD_DIST_H
