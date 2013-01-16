#ifndef WMAP_7YR_OPTIONS_H
#define WMAP_7YR_OPTIONS_H

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

#include <string>
#include <iostream>

#ifdef USE_HIGHELL_TB
#define  num_WMAP   10      // number of individual chi2 terms in likelihood
#else
#define  num_WMAP   8      // number of individual chi2 terms in likelihood
#endif

//---------------------------------------------------
// likelihood terms from WMAP
//---------------------------------------------------
#define ttlike     1 // master tttt chisq flag
#define ttlowllike 2 // low tttt chisq flag
#define ttlowldet  3 // low tttt determinant flag
#define beamlike   4 // beam/pt source correction to tttt chisq flag
#define telike     5 // master tete chisq flag
#define tedet      6 // master tete determinant flag
#define lowllike   7 // TE/EE/BB lowl chisq flag
#define lowldet    8 // TE/EE/BB lowl determinant flag
#define tblike     9 // master tbtb chisq flag
#define tbdet     10 // master tbtb determinant flag


// ===========================================================================
//
class  WMAP_OPTIONS
{
// This module contains the options in the likelihood code
//
// ===========================================================================

 public:
    static WMAP_OPTIONS* self();

    void  setWMAP_data_dir(std::string str) { mWMAP_data_dir = str; }
    std::string  WMAP_data_dir() { return mWMAP_data_dir; }
    void  setTTmax(int ttmax) { mTTmax = ttmax; }                       // must be l.le.1000
    int  ttmax() { return mTTmax; }                       // must be l.le.1000
    void  setTTmin(int ttmin) { mTTmin = ttmin; }                       // must be l.ge.2
    int  ttmin() { return mTTmin; }                       // must be l.ge.2
    void  setTEmax(int temax) { mTEmax = temax; }                       // must be l.le.450
    int  temax() { return mTEmax; }                       // must be l.le.450
    void  setTEmin(int temin) { mTEmin=temin; }                       // must be l.ge.2
    int  temin() { return mTEmin; }                       // must be l.ge.2
    void  setLowl_max(int lowl_max) { mLowl_max=lowl_max; }                 // using low l TT code 2<l<lowl_max (max value lowl_max=12)
    int  lowl_max() { return mLowl_max; }                 // using low l TT code 2<l<lowl_max (max value lowl_max=12)
    void  set_lowl_tt_res(int lowl_tt_res) { mLowl_tt_res=lowl_tt_res; }  // must be either 3 or 4
    int   lowl_tt_res() { return mLowl_tt_res; }
    void  setUsing_lowl_TT(bool b) { mUsing_lowl_TT=b; }      // include TT pixel likelihood, for l<=lowl_max
    bool  using_lowl_TT() { return mUsing_lowl_TT; }      // include TT pixel likelihood, for l<=lowl_max
    void  setUsing_lowl_pol(bool b) { mUsing_lowl_pol=b; }    // include TE,EE,BB pixel likelihood for l<24
    bool  using_lowl_pol() { return mUsing_lowl_pol; }    // include TE,EE,BB pixel likelihood for l<24
    void  setUsing_TT(bool b) { mUsing_TT=b; }                // include MASTER TT in likelihood
    bool  using_TT() { return mUsing_TT; }                // include MASTER TT in likelihood
    void  setUsing_TT_beam_ptsrc(bool b) { mUsing_TT_beam_ptsrc=b; }                // include beam/ptsrc error
    bool  using_TT_beam_ptsrc() { return mUsing_TT; }                // include MASTER TT in likelihood
    void  setUsing_TE(bool b) { mUsing_TE=b; }                // include MASTER TE in likelihood
    bool  using_TE() { return mUsing_TE; }                // include MASTER TE in likelihood

    void setUsing_gibbs_pol_cleaning(bool b) { mUsing_gibbs_pol_cleaning = b; }
    bool using_gibbs_pol_cleaning() { return mUsing_gibbs_pol_cleaning; }

    void setUsing_gibbs(bool b) { mUsing_gibbs = b; }
    bool using_gibbs() { return mUsing_gibbs; }

    void setGibbs_sigma_filename(const std::string& s) { mGibbs_sigma_filename = s; }
    std::string gibbs_sigma_filename() const { return mGibbs_sigma_filename; }

    void setGibbs_first_iteration(long unsigned int i) { mGibbs_first_iteration=i; }
    long unsigned int gibbs_first_iteration() { return mGibbs_first_iteration; }

    void setGibbs_last_iteration(long unsigned int i) { mGibbs_last_iteration=i; }
    long unsigned int gibbs_last_iteration() { return mGibbs_last_iteration; }

    void setGibbs_skip (long unsigned int i) { mGibbs_skip=i; }
    long unsigned int gibbs_skip() { return mGibbs_skip; }

    void setGibbs_ell_max (int l) { mGibbs_ell_max = l; }
    int gibbs_ell_max () { return mGibbs_ell_max; }

    double tt_pixlike_lndet_offset() { return mtt_pixlike_lndet_offset4; }
    double teeebb_pixlike_lndet_offset() { return mteeebb_pixlike_lndet_offset; }
    double te_lndet_offset() { return mte_lndet_offset; }
    double tb_lndet_offset() { return mtb_lndet_offset; }

    void wmap_print_options();

    ~WMAP_OPTIONS() {}
 protected:
    WMAP_OPTIONS();
    static WMAP_OPTIONS* mSelf;
 private:
//---------------------------------------------------
// location of input data
// ---------------------------------------------------
    std::string  mWMAP_data_dir;

//---------------------------------------------------
// l range to be used in the likelihood code
// change these to consider a more limited l range in TTTT and TETE
//---------------------------------------------------
    int  mTTmax                ; //!< must be l.le.1000
    int  mTTmin                ; //!< must be l.ge.2
    int  mTEmax                ; //!< must be l.le.450
    int  mTEmin                ; //!< must be l.ge.2

//---------------------------------------------------
// various likelihood options
// change these to include/ exclude various likelihood aspects
//---------------------------------------------------
    bool  mUsing_lowl_TT        ; //!< include TT pixel likelihood, for l<=lowl_max
    bool  mUsing_lowl_pol       ; //!< include TE,EE,BB pixel likelihood for l<24
    bool  mUsing_TT             ; //!< include MASTER TT in likelihood
    bool  mUsing_TT_beam_ptsrc  ; //!< include beam/ptsrc errors
    bool  mUsing_TE             ; //!< include MASTER TE in likelihood
    bool  mUsing_gibbs_pol_cleaning;

//---------------------------------------------------
// *** AN IMPORTANT CHANGE WITH REGARD TO THE TT LIKELIHOOD ***
//---------------------------------------------------
// There are two options to choose from for evaluating the low-l temperature
// likelihood. Both options produce the same results.
//
// (1) The direct evaluation of likelihood in pixel space using a resolution 4 temperature map.
// (2) The Gibbs sampling.
//
// The option (2) is much faster to evaluate than the option (1).
//
// To use (1), set "use_gibbs = .false." and "lowl_max = 30".
// To use (2), set "use_gibbs = .true." and "lowl_max = 32".
//
// Note that the resolution 3 option for (1) has been disabled.
//

  bool  mUsing_gibbs; // =  true by default

//---------------------------------------------------
// (1) Pixel likelihood
//---------------------------------------------------


  int  mLowl_tt_res            ; // = 4      // TT map resolution
  int  mLowl_max               ; // = 32     // use low l TT code 2<l<lowl_max

//---------------------------------------------------
// (2) Gibbs likelihood
//---------------------------------------------------
// For using different sections of the sigmaElls file,
// adjust gibbs_first_iteration, gibbs_last_iteration,
// and gibbs_skip.
//
// For a 50,000 Gibbs sample file, it may be useful to set
// gibbs_first_iteration = 100
// gibbs_last_iteration = 25000
// gibbs_skip = 3
// for one parameter run (to use every third value from the first half
// (approximately) of the file), and
// gibbs_first_iteration = 25100
// gibbs_last_iteration = 50000
// gibbs_skip = 3
// for another parameter run, to use the second half of the file (every third value).
//
// To get really fast (possibly inaccurate) likelihoods,
// set gibbs_skip to be ~ 0.01 * (gibbs_last_iteration - gibbs_first_iteration)
//
// gibbs_first_iteration must be >= 1

  std::string mGibbs_sigma_filename;
  long unsigned int  mGibbs_first_iteration; // = 10;
  long unsigned int  mGibbs_last_iteration; // = 30000;
  long unsigned int  mGibbs_skip; // = 2;

// The sum in the BR estimator goes up to this value:
  int  mGibbs_ell_max; // = 32 ;



//---------------------------------------------------
// ln(det) offsets
//---------------------------------------------------
// The value of ln(L) returned by the likelihood code is renormalized
// by subtracting off a constant offset:
//
//   -2ln(L) = chi^2 + ln_det_C - ln_det_C_f
//
// The value of the offset, ln_det_C_f,  is the sum of the determinant
// contributions to -2ln(L) computed for the CMB spectrum in
// data/test_cls_v3.dat:
//
//   ln_det_C_f = tt_pixlike_lndet_offset(lowl_tt_res)
//              + teeebb_pixlike_lndet_offset
//              + te_lndet_offset
//
  const double mtt_pixlike_lndet_offset4;
  const double mteeebb_pixlike_lndet_offset;
  const double mte_lndet_offset;
  const double mtb_lndet_offset;
//
// If you'd like to compare the values of ln(L) returned by this version
// of the code with previous versions, use the following:
//
//!  double precision, parameter  tt_pixlike_lndet_offset(3:4) &
//!		= (/5476.3672001139175d0,29467.570953238155d0/)
//!  double precision, parameter  teeebb_pixlike_lndet_offset = 0d0
//!  double precision, parameter  te_lndet_offset = 0d0
//
//  c++ port note: you can do this by changing the constructor in WMAP_7yr_options.cc
};
#endif // WMAP_7YR_OPTIONS_H
