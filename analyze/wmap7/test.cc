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

//===========================================================
//program test_likelihood;

// RB December '05
//===========================================================

#include "WMAP_7yr_likelihood.h"
#include "WMAP_7yr_options.h"

#include "port_helper.h"

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>


int main()
{
  using namespace std;

  using namespace  wmap_likelihood_7yr;
  using namespace  wmap_util;

  string WMAP_data_dir = WMAP_OPTIONS::self()->WMAP_data_dir();
  WMAP_OPTIONS* options = WMAP_OPTIONS::self();

  string filename;
  REAL_8 like_tot, expected_like_tot;
  Array1D<REAL_8> like( Range( 1, num_WMAP ));

  int  tt_npix, teeebb_npix;

  //---------------------------------------------------

  cout << "" << endl;
  cout << "WMAP 7-year likelihood test program" << endl;
  cout << "===================================" << endl;
  cout << "" << endl;
  cout << "NOTE: This code using namespace s a different CMB spectrum than previous versions." << endl;
  cout << "      The new spectrum (data/test_cls_v4.dat) is a better fit to the data" << endl;
  cout << "      than the old one (data/test_cls_v3.dat)." << endl;
  cout << "" << endl;
  cout << "      As before, a constant offset is now being subtracted from ln(L)." << endl;
  cout << "      The value of the offset is the sum of the determinant" << endl;
  cout << "      contributions to ln(L) computed for the CMB spectrum in" << endl;
  cout << "      data/test_cls_v4.dat, ln_det_C_f:" << endl;
  cout << "" << endl;
  cout << "        -2ln(L) = chi^2 + ln_det_C - ln_det_C_f" << endl;
  cout << "" << endl;


  //---------------------------------------------------

  const unsigned int ttmax = WMAP_OPTIONS::self()->ttmax();
  const unsigned int ttmin = WMAP_OPTIONS::self()->ttmin();
  const unsigned int lowl_max = WMAP_OPTIONS::self()->lowl_max();
  const unsigned int temax = WMAP_OPTIONS::self()->temax();
  real8_1d cl_tt(Range(ttmin,ttmax));
  real8_1d cl_te(Range(ttmin,ttmax));
  real8_1d cl_ee(Range(ttmin,ttmax));
  real8_1d cl_bb(Range(ttmin,ttmax));

  cl_bb = 0.0;

#ifdef USE_LOWELL_TBEB
  real8_1d cl_tb(Range(ttmin,ttmax));
  real8_1d cl_eb(Range(ttmin,ttmax));
  cl_tb = 0.;
  cl_eb = 0.;
#endif

  //---------------------------------------------------
  // read in Cls
  //---------------------------------------------------
  filename = WMAP_data_dir+"test_cls_v4.dat";

  cout << "Reading in Cls from:  " << filename << endl;

  ifstream lun( filename.c_str() );

  int i,l;
  for (unsigned int l=2; l <= ttmax; ++l)
  {
    lun >> i >> cl_tt(l) >> cl_ee(l) >> cl_bb(l) >>cl_te(l);
  }

  lun.close();


  //---------------------------------------------------
  // put in likelihood options
  // see PASS2_options module for the options below
  //---------------------------------------------------

  options->setUsing_TT(true);
  options->setUsing_TE(true);
  options->setUsing_lowl_TT(true);
  options->setUsing_lowl_pol(true);

  //---------------------------------------------------
  // get likelihoods
  //---------------------------------------------------
  like = 0.;
  wmap_likelihood_init();
  wmap_likelihood_dof( tt_npix, teeebb_npix );
#ifdef USE_LOWELL_TBEB
  wmap_likelihood_compute(cl_tt,cl_te,cl_tb,cl_ee,cl_eb,cl_bb,like);
#else
  wmap_likelihood_compute(cl_tt,cl_te,cl_ee,cl_bb,like);
#endif
  wmap_likelihood_error_report();

  like_tot = like.sum();
  //---------------------------------------------------
  // write outputs
  //---------------------------------------------------
  cout.precision( 16 );
  cout << "------------------------------------------------------------------" << endl;
  cout << "Breakdown of -2ln(L)" << endl;
  cout << "------------------------------------------------------------------" << endl;
  cout << "MASTER TTTT             " << (2*like(ttlike)) << " for " << (ttmax-lowl_max) << endl;
  cout << "Beam/ptsrc TT correction" << (2*like(beamlike)) << endl;
  if (options->using_gibbs()) {
    cout << "low-l TTTT gibbs        " << (2*like(ttlowllike)) << endl;
  } else {
    cout << "low-l TTTT chi2         " << (2*like(ttlowllike)) << " for " << tt_npix << endl;
    cout << "low-l TTTT det          " << (2*like(ttlowldet))   << endl;
  }
  cout << "MASTER TETE chi2        " << (2*like(telike)) << " for " << (temax-23.) << endl;;
  cout << "MASTER TETE det         " << (2*like(tedet)) << endl;
  cout << "TT/TE/EE/BB low-l chi2  " << (2*like(lowllike)) << " for " << teeebb_npix << endl;
  cout << "TT/TE/EE/BB low-l det   " << (2*like(lowldet)) << endl;
  cout << "------------------------------------------------------------------" << endl;
  cout <<  "TOTAL -2ln(L)           " << (2*like_tot) << endl;
  cout << "------------------------------------------------------------------" << endl;

  if (options->using_gibbs()) {
    expected_like_tot = 7481.131790;
  } else {
#ifdef FASTERTT
    expected_like_tot = 8257.119952;
#else
    expected_like_tot = 11068.324977;
#endif
  }

  cout << "Expected -2ln(L)         = " << expected_like_tot << endl;
  cout << "Difference         = " << (2*like_tot-expected_like_tot) << endl;
  cout << "------------------------------------------------------------------" << endl;
  cout <<  "" << endl ;
  cout <<  "Differences on the order of O(0.001) are normal between platforms." << endl ;
  cout <<  "" << endl ;

} // end program test_likelihood;
