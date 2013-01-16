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

#include "WMAP_7yr_gibbs.h"

#include "WMAP_7yr_options.h"
#include "WMAP_7yr_util.h"
#include "br_mod_dist.h"

#include <string>
#include <fstream>

namespace wmap_gibbs
{
  using namespace br_mod_dist;
  using namespace wmap_util;


  namespace
  {
    real8_1d  cl_tt_fiducial, cl_tt_dummy;
    static bool gibbs_verbose = false;
//X     static bool gibbs_verbose = true;
  }

void  setup_for_tt_gibbs()
{
    using namespace std; // cout, ifstream

    long int  ell_max, num_chains, num_samples, ell_min, ell, i;
    ifstream lun;
    FloatToDoubleArray3D  sigmas;
    double  ln_likelihood, dummy_te, dummy_ee, dummy_bb;
    std::string filename;

    WMAP_OPTIONS* options = WMAP_OPTIONS::self();

    if ( !options->using_lowl_TT())  {
       cout <<  "Error: using_lowl_TT ==  false  in setup_for_tt_gibbs()" << endl ;
       cout <<  "Set using_lowl_TT =  true " << endl ;
       exit(-1);
    }

    if ( !options->using_gibbs())  {
       cout <<  "Error: use_gibbs ==  false  in setup_for_tt_gibbs()" << endl ;
       cout <<  "Set use_gibbs =  true " << endl ;
       exit(-1);
    }

    if (options->lowl_max()  <  2)  {
       cout <<  "Error: lowl_max is less than 2." << endl ;
       cout <<  "Raise it to 2 or greater for use with Gibbs sampling." << endl ;
       exit(-1);
    }

    //---------------------------------------------------
    // read in sigmas and initialize br_mod
    //---------------------------------------------------

    string gibbs_sigma_filename=options->gibbs_sigma_filename();
    if (gibbs_verbose)  {
      cout <<  "------------------------------------------------------------" << endl ;
      cout <<  "Reading from "+gibbs_sigma_filename << endl;
    }

    read_gibbs_chain(gibbs_sigma_filename, ell_max,
                     num_chains, num_samples, sigmas);

    if ( gibbs_verbose )  {
      cout <<  "minval sigmas = " << *std::min_element(sigmas.data()+2, sigmas.data()+sigmas.size()) << endl ;
      cout <<  "maxval sigmas = " << *std::max_element(sigmas.data()+2, sigmas.data()+sigmas.size()) << endl ;

      cout <<  "Read in sigmas for Gibbs likelihood:" << endl ;
      cout <<  "Maximum ell = " << ell_max << endl ;
      cout <<  "Number of chains = " << num_chains << endl ;
      cout <<  "Number of Gibbs samples = " << num_samples << endl ;
    }

    if (options->lowl_max() > ell_max)  {
       cout <<  "Error: lowl_max is too large for using namespace with "+gibbs_sigma_filename << endl ;
       cout <<  "lowl_max: " << options->lowl_max() << endl ;
       cout <<  "maximum allowed value: " << ell_max << endl ;
       exit(-1);
    }

    long unsigned int gibbs_first_iteration = options->gibbs_first_iteration();
    long unsigned int gibbs_last_iteration = options->gibbs_last_iteration();
    long unsigned int gibbs_skip = options->gibbs_skip();
    int gibbs_ell_max = options->gibbs_ell_max();

    if (gibbs_ell_max > ell_max)  {
       cout <<  "Error: gibbs_ell_max is too large for using namespace with " << gibbs_sigma_filename << endl ;
       cout <<  "gibbs_ell_max: " << gibbs_ell_max << endl ;
       cout <<  "maximum allowed value: " << ell_max << endl ;
       exit(-1);
    }

    if ( gibbs_verbose )  {
      cout <<  "Using values of ell from 2 to " << gibbs_ell_max << endl ;
      cout <<  "Using Gibbs samples from " << gibbs_first_iteration << " to " << gibbs_last_iteration << endl ;
      cout <<  "in steps of " << gibbs_skip << endl ;
      cout <<  "------------------------------------------------------------" << endl ;
    }

    ell_min = 2;
    int num_samples_tmp = (int)(gibbs_last_iteration-gibbs_first_iteration)/gibbs_skip;
    double_3d sigmas_tmp;
    sigmas_tmp.resize(Range(2,gibbs_ell_max), Range(1,num_chains), Range(1,num_samples_tmp+1));

    for (int i1 = 2; i1 <= gibbs_ell_max; ++i1) {
      for (int i2 = 1; i2 <= num_chains; ++i2) {
        int i3_new = 1;
        for (unsigned int i3 = gibbs_first_iteration; i3 <= gibbs_last_iteration; i3 += gibbs_skip) {
          sigmas_tmp(i1, i2, i3_new++) = sigmas(i1, i2, i3);
        }
      }
    }

    initialize_br_mod(ell_min, sigmas_tmp);

    //---------------------------------------------------
    // read in Cls
    //---------------------------------------------------
    cl_tt_fiducial.resize(Range(options->ttmin(),options->ttmax()));
    cl_tt_dummy.resize(Range(options->ttmin(),options->gibbs_ell_max()));

    filename = options->WMAP_data_dir()+"lowlT/gibbs/test_cls.dat";

    if ( gibbs_verbose )  {
      cout << "Reading in Cls from: " << filename << endl;
    }

    lun.open(filename.c_str());
    for (ell = options->ttmin(); ell <= options->ttmax(); ++ell ) {
       lun >> i >> cl_tt_fiducial(ell) >> dummy_te >>  dummy_ee >> dummy_bb;
       if (i  !=  ell)  {
          cout <<  "Error with file format of: " << filename << endl ;
          exit(-1);
       }
    }
    lun.close();

    //---------------------------------------------------
    // Initialize br_mod with default spectrum
    //---------------------------------------------------
    for (int l=options->ttmin(); l<=gibbs_ell_max; ++l)
      cl_tt_dummy(l) = cl_tt_fiducial(l);
    compute_br_estimator(cl_tt_dummy, ln_likelihood);
    if ( gibbs_verbose )  {
      cout <<  "Initialized log likelihood = " << ln_likelihood << endl ;
    }

  } //  setup_for_tt_gibbs



  // compute_br_estimator() returns the natural logarithm of the 
  // likelihood (plus some arbitrary constant).  The code expects
  // the negative of this quantity (= chisquared/2)
  void  compute_tt_gibbslike(const real8_1d& cl_in, double& like)
  {
#ifdef TIMING
    wmap_timing_start("compute_tt_gibbslike");
#endif
    WMAP_OPTIONS *options = WMAP_OPTIONS::self();

    if (options->gibbs_ell_max() > options->lowl_max())  {
      for (int l=options->lowl_max()+1; l <= options->gibbs_ell_max(); ++l)
       cl_tt_dummy(l) = cl_tt_fiducial(l);
    }
    for (int l=options->ttmin(); l <= options->lowl_max(); ++l)
      cl_tt_dummy(l) = cl_in(l);
    compute_br_estimator(cl_tt_dummy, like);
    like = -like;
#ifdef TIMING
    wmap_timing_end();
#endif
  } //  compute_tt_gibbslike


} // end namespace  wmap_gibbs
