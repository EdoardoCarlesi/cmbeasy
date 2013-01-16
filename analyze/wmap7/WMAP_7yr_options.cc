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

#include "WMAP_7yr_options.h"

#include <iostream>

WMAP_OPTIONS* WMAP_OPTIONS::mSelf = NULL;

WMAP_OPTIONS::WMAP_OPTIONS()
   :
   mWMAP_data_dir   ("./data/"),
   mTTmax           (1200),              // must be l.le.1200
   mTTmin           (2),                 // must be l.ge.2
   mTEmax           (800),               // must be l.le.800
   mTEmin           (2),                 // must be l.ge.2
   mUsing_lowl_TT   (true),              // include TT pixel likelihood, for l<=lowl_max
   mUsing_lowl_pol  (true),              // include TE,EE,BB pixel likelihood for l<24
   mUsing_TT        (true),              // include MASTER TT in likelihood
   mUsing_TT_beam_ptsrc(true),           // include beam/ptsrc errors
   mUsing_TE        (true),              // include MASTER TE in likelihood
   mUsing_gibbs_pol_cleaning(false),
   mUsing_gibbs     (true),
//X    mUsing_gibbs     (false),
   mLowl_tt_res      (4),                // TT map resolution
   mLowl_max        (32),                // use low l tt code 2<l<lowl_max
//X    mLowl_max        (30),                // use low l tt code 2<l<lowl_max
   mGibbs_sigma_filename (mWMAP_data_dir+"lowlT/gibbs/sigmaEllsHkeChu_run4_kq85.fits"),
   mGibbs_first_iteration  (10),
   mGibbs_last_iteration   (120000),
   mGibbs_skip             (2),
   mGibbs_ell_max          (32),
#ifdef FASTERTT
   mtt_pixlike_lndet_offset4(5024.741512),
#else
   mtt_pixlike_lndet_offset4(-29677.056620),
#endif
   mteeebb_pixlike_lndet_offset( 16078.083180),
   mte_lndet_offset(3584.277805),
   mtb_lndet_offset(3598.152208)
{}

WMAP_OPTIONS* WMAP_OPTIONS::self()
{
  if ( !mSelf )
    mSelf = new WMAP_OPTIONS();
  return mSelf;
}

void WMAP_OPTIONS::wmap_print_options()
{
  using namespace std;
  cout <<  "-----------------------------------------------------" << endl ;
  cout <<  "WMAP_data_dir = " << WMAP_data_dir()  << endl ;
  cout <<  endl ;
  cout <<  "ttmax = " << ttmax()  << endl ;
  cout <<  "ttmin = " << ttmin()  << endl ;
  cout <<  "temax = " << temax() << endl ;
  cout <<  "temin = " << temin() << endl ;
  cout <<  "" << endl ;
  cout <<  "using_lowl_TT =       " << using_lowl_TT() << endl ;
  cout <<  "using_lowl_pol =      " << using_lowl_pol()<< endl ;
  cout <<  "using_TT =            " << using_TT() << endl ;
  cout <<  "using_TT_beam_ptsrc = " << using_TT_beam_ptsrc() << endl ;
  cout <<  "using_TE =            " << using_TE() << endl ;
  cout <<  "" << endl ;
  cout <<  "lowl_tt_res = " << lowl_tt_res() << endl ;
  cout <<  "lowl_max =    " << lowl_max() << endl ;
  cout <<  "" << endl ;
  cout <<  "tt_pixlike_lndet_offset =     " << tt_pixlike_lndet_offset() << endl ;
  cout <<  "teeebb_pixlike_lndet_offset = " << teeebb_pixlike_lndet_offset() << endl ;
  cout <<  "te_lndet_offset =             " << te_lndet_offset() << endl ;
#ifdef use_HIGHELL_TB
  cout <<  "tb_lndet_offset =             " << tb_lndet_offset() << endl ;
#endif
  cout <<  "" << endl ;
  cout <<  "using_gibbs = " <<  using_gibbs() << endl ;
  cout <<  "gibbs_sigma_filename = "  <<  gibbs_sigma_filename() << endl ;
  cout <<  "gibbs_first_iteration = " <<  gibbs_first_iteration() << endl ;
  cout <<  "gibbs_last_iteration =  " <<  gibbs_last_iteration() << endl ;
  cout <<  "gibbs_skip =            " <<  gibbs_skip() << endl ;
  cout <<  "gibbs_ell_max =         " <<  gibbs_ell_max() << endl ;
  cout <<  "-----------------------------------------------------" << endl ;
}
