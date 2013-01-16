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

// ===========================================================================

#include "WMAP_7yr_likelihood.h"
#include "WMAP_7yr_options.h"
#include "WMAP_7yr_gibbs.h"

#ifndef WMAP_STANDALONE
#include "global.h"  // for Bad_Error
#endif

#include "gsl/gsl_errno.h"

#include <string>
#include <fstream>

namespace wmap_likelihood_7yr
{

// This code is the central likelihood routine from which the other subroutines
// are called.
//
// Parameters are defined in associated WMAP_options module
//
// This code was a collaborative effort by the following WMAP Team Members:
// R. Bean
// J. Dunkley
// E. Komatsu
// D. Larson
// M. Nolta
// H. Peiris
// L. Verde
// D. Spergel
//
// ===========================================================================

using namespace  wmap_tt_beam_ptsrc_chisq;
#ifdef USE_LOWELL_TBEB
  using namespace  wmap_tetbeebbeb_lowl;
#else
  using namespace  wmap_teeebb_lowl;
#endif
  using namespace  wmap_tlike;
  using namespace  wmap_gibbs;

namespace
{
static const std::string wmap_likelihood_version="v4_c++";

static const int  tt_lmax = 1200;
static const int  te_lmax =  800;
  
static bool  initialise_pass2 = true;

real8_2d R_off_tete(Range( 2,te_lmax ), Range( 2,te_lmax) );
real8_2d R_off_tttt(Range(2,tt_lmax ),Range( 2,tt_lmax));
real8_2d off_log_curv(Range(2,tt_lmax), Range(2,tt_lmax));
real8_2d epsilon(Range(2,tt_lmax),Range(2,tt_lmax));

real8_1d cltt_dat(Range(2,tt_lmax));
real8_1d clte_dat(Range(2,te_lmax));
real8_1d cltb_dat(Range(2,te_lmax));
real8_1d ntt(Range(2,tt_lmax));
real8_1d ntt_te(Range(2,te_lmax));
real8_1d nee_te(Range(2,te_lmax));
real8_1d fskytt(Range(2,tt_lmax));
real8_1d fskyte(Range(2,te_lmax));
real8_1d ntt_tb(Range(2,te_lmax));
real8_1d nbb_tb(Range(2,te_lmax));
real8_1d fskytb(Range(2,te_lmax));
}


void cmbeasy_gsl_error_handler( const char * reason, const char * file, int line, int gsl_errno );


// ===========================================================================
void  wmap_likelihood_init()
{

  // ===========================================================================

  using namespace wmap_util;
  using namespace std;

  WMAP_OPTIONS* options = WMAP_OPTIONS::self();

  int   l,ll;
  int   i,j;
  REAL_8 dummy;
  std::string ttfilename,tefilename,tbfilename,ttofffilename,teofffilename;

#ifdef TIMING
  wmap_timing_start( "wmap_likelihood_init" );
#endif

  std::cout <<  "Initializing WMAP likelihood, version " << wmap_likelihood_version << std::endl;


  //-----------------------------------------------
  // initialise beam uncertainties
  //-----------------------------------------------

  if ( options->using_TT_beam_ptsrc() ) init_tt_beam_and_ptsrc_chisq( 2, tt_lmax );

  //-----------------------------------------------
  // initialise low l codes
  //-----------------------------------------------
  if (options->using_lowl_pol()) {
        // cout << options->using_lowl_pol << " using low ell" << endl;
#ifdef USE_LOWELL_TBEB
        tetbeebbeb_lowl_like_setup();
#else
        teeebb_lowl_like_setup();
#endif
  }

  if (options->using_lowl_TT()) {
    if (options->using_gibbs()) {
      setup_for_tt_gibbs();
    } else {
      setup_for_tt_exact(options->lowl_max());
    }
  }


  //-----------------------------------------------
  // set file names
  //-----------------------------------------------

  std::string WMAP_data_dir = options->WMAP_data_dir();
  ttfilename = WMAP_data_dir+"highl/wmap_likelihood_inputs_tt.p4v4.dat";
  tefilename = WMAP_data_dir+"highl/wmap_likelihood_inputs_te_recalibrated.p4v4.dat";
  tbfilename = WMAP_data_dir+"highl/wmap_likelihood_inputs_tb_recalibrated.p4v4.dat";
  ttofffilename = WMAP_data_dir+"highl/wmap_likelihood_inputs_tt_offdiag.p4v4.dat";
  teofffilename = WMAP_data_dir+"highl/wmap_likelihood_inputs_te_offdiag.p4v4.dat";

  //-----------------------------------------------
  // get TT diag
  //-----------------------------------------------

  ifstream lun ( ttfilename.c_str() );
  if ( !lun )
  {
    cout << "cant find" << ttfilename << endl;
    exit(-1);
  }

  for (l=2; l<= tt_lmax; ++l)
  {
    lun >> dummy >> cltt_dat(l) >> ntt(l) >> fskytt(l);
  }
  lun.close();


  //-----------------------------------------------
  // get TE diag
  //-----------------------------------------------

  lun.open( tefilename.c_str() );
  if ( !lun )
  {
    cout <<  "cant find " << tefilename << endl;
    exit(-1);
  }

  for (l=2; l<= te_lmax; ++l)
  {
    lun >> ll >> clte_dat(l) >>  dummy >>  ntt_te(l) >> nee_te(l) >> fskyte(l);
  }
  lun.close();

  //-----------------------------------------------
  // get TB diag
  //-----------------------------------------------

#ifdef USE_HIGHELL_TB
  lun.open( tbfilename.c_str() );
  if ( !lun )
  {
    cout <<  "cant find " << tbfilename << endl;
    exit(-1);
  }

  for (l=2; l<= te_lmax; ++l)
  {
    lun >> ll >> cltb_dat(l) >>  dummy >>  ntt_tb(l) >> nbb_tb(l) >> fskytb(l);
  }
  lun.close();
#endif

  //-----------------------------------------------
  // get TT off diag
  //-----------------------------------------------

  R_off_tttt = 0.;

  lun.open( ttofffilename.c_str() );

  if ( !lun)
  {
    cout << "cant find " << ttofffilename << endl;
    exit(-1);
  }

  for (i=2; i<= tt_lmax; ++i)
  {
    for (j=i+1; j<= tt_lmax; ++j)
    {

      lun >> l >> ll >> epsilon(i,j) >> R_off_tttt(i,j);
      R_off_tttt(j,i)=R_off_tttt(i,j);
      epsilon(j,i)=epsilon(i,j);

      if (i != l || j != ll)
      {
        cout << "tt off file misread " << i << " " << j << " " << l << " " << ll << endl;
        exit(-1);
      }

    }
  }

  lun.close();

  //-----------------------------------------------
  // get TE off diag
  //----------------------------------------------

  lun.open( teofffilename.c_str() );

  if ( !lun )
  {
    cout <<  "cant find " << teofffilename << endl;
    exit(-1);
  }

  for (i=2; i<= te_lmax-1; ++i)
  {
    for (j=i+1; j<= te_lmax; ++j)
    {

      lun >> l >> ll >> dummy;

      if (l <= te_lmax && ll <= te_lmax)
      {
        R_off_tete(i,j)=dummy;
        R_off_tete(j,i)=R_off_tete(i,j);
      }

      if (l != i || ll != j)
      {
        cout << "TE off diag misread i,j,l,ll " << i << " " << j << " " << l << " " << ll << endl;
        exit(-1);
      }

    }
  }
  lun.close();

  //X set our own GSL error handler
  gsl_set_error_handler ( &wmap_likelihood_7yr::cmbeasy_gsl_error_handler );

  initialise_pass2 =  false;

#ifdef TIMING
  wmap_timing_end();
#endif

  return;
} //  wmap_likelihood_init

void cmbeasy_gsl_error_handler( const char * reason, const char * file, int line, int gsl_errno )
{
  std::string errorString;
  errorString += "An error occured while computing the WMAP3 likelihood at: ";
  errorString += file;
  errorString += ", line ";
  errorString += line;
  errorString += ".\nThe reason was: ";
  errorString += reason;
  errorString += "\nGSL error no";
  errorString += gsl_errno;
  errorString += '\n';
  throw Bad_Error( errorString );
}

void  wmap_likelihood_dof( int& tt_npix, int& teeebb_npix )
{
    tt_npix = tt_pixlike_dof();
#ifdef USE_LOWELL_TBEB
    teeebb_npix = tetbeebbeb_pixlike_dof();
#else
    teeebb_npix = teeebb_pixlike_dof();
#endif
} // wmap_likelihood_dof



// ===========================================================================
#ifdef USE_LOWELL_TBEB
void wmap_likelihood_compute(real8_1d& cltt, real8_1d& clte, real8_1d& cltb, real8_1d&clee,
                             real8_1d& cleb, real8_1d& clbb, Array1D<REAL_8>& like)
#elif USE_HIGHELL_TB
void wmap_likelihood_compute(real8_1d& cltt, real8_1d& clte, real8_1d& cltb, real8_1d& clee,
                             real8_1d& cleb, real8_1d& clbb, Array1D<REAL_8>& like)
#else
void wmap_likelihood_compute(real8_1d& cltt, real8_1d& clte, real8_1d& clee, real8_1d& clbb,
                             Array1D<REAL_8>& like)
#endif
{

  // ===========================================================================

  using namespace  wmap_util;
  WMAP_OPTIONS* options = WMAP_OPTIONS::self();

#ifdef USE_LOWELL_TBEB
  REAL_8 ee2, bb2, eb2;
#elif USE_HIGHELL_TB
  REAL_8 ee2, bb2, eb2;
#endif

  int  il, ill;
  REAL_8    dlnlike_tot, dlnlike, dlnlike_beam, ln_det_TETE, ln_det_TBTB;
  real8_2d  fisher(Range( 2,tt_lmax ), Range( 2,tt_lmax) );
  real8_1d  tttt(Range(2,tt_lmax)), tete(Range(2,te_lmax)), ttte(Range(2,te_lmax)), tbtb(Range(2,te_lmax));
  real8_1d  z(Range(2,tt_lmax )), zbar(Range(2,tt_lmax) );
  real8_1d  cltt_temp(Range(2,tt_lmax) );
  double_1d  lowl_cl;
  int  tt_hi_l_start, te_hi_l_start;
  REAL_8   correlation_coefficient_cl;
  static const REAL_8 tol=1e-10;

#ifdef TIMING
  wmap_timing_start( "wmap_likelihood_compute" );
#endif

  wmap_likelihood_error_init();
  dlnlike_tot = 0.;

  if (initialise_pass2)
  {
    wmap_likelihood_init();
  }

  like = 0.0;

  //--------------------------------------------------------
  // Are cltt, clte, and clee consistent?
  //--------------------------------------------------------
  for (il=2; il<= options->ttmax(); ++il)
  {
    if ( fabs(clte(il)) > 0.0 )  {
      correlation_coefficient_cl = fabs(clte(il))/sqrt(cltt(il)*clee(il));
      if (correlation_coefficient_cl-1.0 > tol)
      {
        wmap_likelihood_error("unphysical input: TE/sqrt(TT*EE) > 1 at l=", il);
        return;
        //      exit(-1);
      }
    }
  }

#ifdef USE_HIGHELL_TB
  for (il=2; il<= options->ttmax(); ++il)
  {
    if ( fabs(cltb(il)) > 0.0 )  {
      correlation_coefficient_cl = fabs(cltb(il))/sqrt(cltt(il)*clbb(il));
      if (correlation_coefficient_cl-1.0 > tol)
      {
        wmap_likelihood_error("unphysical input: TB/sqrt(TT*BB) > 1 at l=", il);
        return;
      }
    }
  }
#endif // USE_HIGHELL_TB

#ifdef USE_LOWELL_TBEB
  for (il=2; il<= options->ttmax(); ++il)
  {
    if ( fabs(cltb(il)) > 0.0 )  {
      correlation_coefficient_cl = fabs(cltb(il))/sqrt(cltt(il)*clbb(il));
      if (correlation_coefficient_cl-1.0 > tol)
      {
        wmap_likelihood_error("unphysical input: TB/sqrt(TT*BB) > 1 at l=", il);
        return;
      }
    }
  }

  for (il=2; il<= options->ttmax(); ++il)
  {
    if ( fabs(cleb(il)) > 0.0 )  {
      correlation_coefficient_cl = fabs(cleb(il))/sqrt(clee(il)*clbb(il));
      if (correlation_coefficient_cl-1.0 > tol)
      {
        using namespace std;
        cout << correlation_coefficient_cl << endl;
        cout << "unphysical input: EB/sqrt(EE*BB) > 1 at l=" << il << endl;
        wmap_likelihood_error("WMAP5 Likelihood: unphysical input: EB/sqrt(EE*BB) at l=", il);
        return;
      }
    }

    if (fabs(cleb(il)) > 0.0) {
      ee2 = clee(il) - clte(il)*clte(il)/cltt(il);
      bb2 = clbb(il) - cltb(il)*cltb(il)/cltt(il);
      eb2 = cleb(il) - clte(il)*cltb(il)/cltt(il);
      if (fabs(eb2)/sqrt(ee2*bb2)-1.0 > tol )  {
        wmap_likelihood_error("unphysical input: EB/sqrt(EE*BB) > 1 at l=", il);
        return;
      }
    }
  }

#endif // USE_LOWELL_TBEB



  cltt_temp = cltt(Range(2, options->ttmax()));

  //---------------------------------------------------------------------------
  // low l TT likelihood
  //---------------------------------------------------------------------------
  if (options->using_lowl_TT())
  {
    if (options->using_gibbs())  {
      compute_tt_gibbslike(cltt_temp(Range(2,options->lowl_max())), like(ttlowllike));
      like(ttlowldet) = 0.0;
    } else {
      lowl_cl.resize(Range(2, options->lowl_max()));
      lowl_cl = cltt_temp(Range(2, options->lowl_max()));
      compute_tt_pixlike(options->lowl_max(),lowl_cl,like(ttlowllike),like(ttlowldet));
      lowl_cl.resize(0);
    }
    tt_hi_l_start = options->lowl_max()+1;
  } else {
    tt_hi_l_start = options->ttmin();
  }

  //---------------------------------------------------------------------------
  // low l TE/EE/BB likelihood
  //---------------------------------------------------------------------------
  if (options->using_lowl_pol())
  {
#ifdef USE_LOWELL_TBEB
    tetbeebbeb_lowl_likelihood(23, cltt_temp,clte,cltb,clee,clbb,cleb,like(lowllike),like(lowldet));
#else
    teeebb_lowl_likelihood(23, cltt_temp,clte,clee,clbb,like(lowllike),like(lowldet));
#endif
    te_hi_l_start = 24;
  }
  else
  {
    te_hi_l_start = options->temin();
  }

  //---------------------------------------------------------------------------
  // TT and TE covariance terms
  //---------------------------------------------------------------------------
  tttt=0.;
  tete=0.;
  ttte=0.;

  if (options->using_TT() || options->using_TE())
  {

    for (il=options->ttmin(); il<= options->ttmax(); ++il)
    {
      tttt(il)= 2.*pow( cltt_temp(il)+ ntt(il),2 ) /((2.0*double(il)+1.0)*pow( fskytt(il),2 ));
    }

    for (il=options->temin(); il<= options->temax(); ++il)
    {
      tete(il)= ((cltt_temp(il)+ntt_te(il))*(clee(il)+nee_te(il))+pow( clte(il),2 ))/
                ((2.0*double(il)+1.0)*fskyte(il)*fskyte(il));
      ttte(il)= 2.*((cltt_temp(il)+ntt(il))*clte(il))/
                ((2.0*double(il)+1.0)*fskyte(il)*fskytt(il));
    }

  }


  //---------------------------------------------------------------------------
  // TB covariance terms (diagonal only)
  //---------------------------------------------------------------------------
  tbtb=0.0;

#ifdef USE_HIGHELL_TB
     for (il=options->temin(); il<= options->temax(); ++il) {
        tbtb(il)= ((cltt_temp(il)+ntt_tb(il))*(clbb(il)+nbb_tb(il))+cltb(il)*cltb(il))/
             ((2.0*double(il)+1.0)*fskytb(il)*fskytb(il));
     }
#endif


  //---------------------------------------------------------------------------
  //TTTT MASTER likelihood
  //---------------------------------------------------------------------------
  if (options->using_TT())
  {

    fisher=0.;
    off_log_curv=0.;

    for (il=options->ttmin(); il<= options->ttmax(); ++il)
    {
      z(il)=log(cltt_dat(il)+ntt(il));
      zbar(il)=log(cltt_temp(il)+ntt(il));
    }

    for (il=options->ttmin(); il<= options->ttmax(); ++il)
    {
      for (ill=il; ill<= options->ttmax(); ++ill)
      {
        if (il == ill)
        {
          if (il <= te_lmax)
          {
            fisher(il,ill) = tete(il)/(tttt(il)*tete(il)-ttte(il)*ttte(il));
          }
          else
          {
            fisher(il,ill) = 1./tttt(il);
          }
        }
        else
        {
          fisher(il,ill) = -R_off_tttt(il,ill)/sqrt(tttt(il)*tttt(ill))
                           +epsilon(il,ill)/(tttt(il)*tttt(ill));
        }
        off_log_curv(il,ill)=(cltt_temp(il)+ntt(il))*fisher(il, ill)*
                             (cltt_temp(ill)+ntt(ill));
      }
    }


    for (il=options->ttmin(); il<= options->ttmax(); ++il)
    {
      for (ill=il; ill<= options->ttmax(); ++ill)
      {
        dlnlike = 2./3.*(z(il)-zbar(il))*off_log_curv(il,ill)*(z(ill)-zbar(ill))+
                  1./3.*(cltt_temp(il)-cltt_dat(il))*
                  fisher(il,ill)*(cltt_temp(ill)-cltt_dat(ill));

        if (il >= tt_hi_l_start || ill >= tt_hi_l_start)
        {
          if (il == ill)
          {
            dlnlike_tot = dlnlike_tot + dlnlike;
          }
          else
          {
            dlnlike_tot = dlnlike_tot + dlnlike*2;
          }
        }

      }
    }

    like(ttlike) = dlnlike_tot/2.;

  }
  //---------------------------------------------------------------------------
  //TTTT Beam and point source correction
  //---------------------------------------------------------------------------
  if (options->using_TT() && options->using_TT_beam_ptsrc())
  {
    dlnlike_beam = compute_tt_beam_and_ptsrc_chisq( options->ttmin(), options->ttmax(),
                   cltt_temp, cltt_dat, ntt, fisher, z, zbar );
    if ( fabs(dlnlike_beam) >= dlnlike_tot/4. )
    {
      wmap_likelihood_warning( "beam correction invalid", 0 );
      dlnlike_tot = 0.;
    }
    like(beamlike) = dlnlike_beam/2.;
  }

  //---------------------------------------------------------------------------
  //TETE MASTER likelihood
  //---------------------------------------------------------------------------
  if (options->using_TE())
  {
    ln_det_TETE=0.0;
    dlnlike_tot=0.;
    fisher=0.     ;

    for (il=options->temin(); il<= options->temax(); ++il)
    {
      for (ill=il; ill<= options->temax(); ++ill)
      {
        if (il == ill)
        {
          if (il <= te_lmax)
          {
            fisher(il,ill) = tttt(il)/(tttt(il)*tete(il)-ttte(il)*ttte(il));
          }
          else
          {
            fisher(il,ill) = 1./tete(il);
          }
        }
        else
        {
          fisher(il,ill) = -R_off_tete(il,ill)/sqrt(tete(il)*tete(ill));
        }
      }
    }

    for (il=options->temin(); il<= options->temax(); ++il)
    {
      for (ill=il; ill<= options->temax(); ++ill)
      {

        dlnlike = (clte(il)-clte_dat(il))*fisher(il,ill)* (clte(ill)-clte_dat(ill));

        if (il >= te_hi_l_start || ill >= te_hi_l_start)
        {
          if (il == ill)
          {
            dlnlike_tot=dlnlike_tot+dlnlike;
            ln_det_TETE=ln_det_TETE-log(fisher(il,ill));
          }
          else
          {
            dlnlike_tot=dlnlike_tot+dlnlike*2.;
          }
        }

      }
    }
    like(telike) = dlnlike_tot/2.;

    like(tedet) = (ln_det_TETE - WMAP_OPTIONS::self()->te_lndet_offset())/2.;
  }

  //---------------------------------------------------------------------------
  //TBTB MASTER likelihood
  //---------------------------------------------------------------------------
#ifdef USE_HIGHELL_TB
     ln_det_TBTB=0.0;
     dlnlike_tot=0.0;
     fisher=0.0;

     for (il=options->temin(); il<= options->temax(); ++il) {
       fisher(il,il) = 1.0/tbtb(il);
     }

     for (il=options->temin(); il<= options->temax(); ++il) {
       dlnlike = pow((cltb(il)-cltb_dat(il)),2.)*fisher(il,il);

       if (il >= te_hi_l_start)  {
         dlnlike_tot=dlnlike_tot+dlnlike;
         ln_det_TBTB=ln_det_TBTB-log(fisher(il,il));
       }

     }
     like(tblike) = dlnlike_tot/2.0;

     like(tbdet) = (ln_det_TBTB - options->tb_lndet_offset())/2.0;
//X      cout << ln_det_TBTB << endl;
#endif

  //X 10 continue;
  //
#ifdef TIMING
	wmap_timing_end();
#endif

} // wmap_likelihood_compute


} // end namespace wmap_likelihood_7yr

