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
//=========================================
// Module to calculated the correction to the TTTT chisq
// from beam and point source corrections
//
// Written by Mike Nolta
//=========================================

#include "WMAP_7yr_tt_beam_and_ptsrc_chisq.h"
#include "WMAP_7yr_util.h"
#include "WMAP_7yr_options.h"

#include "gsl/gsl_linalg.h"

#include <string>
#include <fstream>

namespace  wmap_tt_beam_ptsrc_chisq
{

using namespace std;
namespace
{    //X simulate 'private'

static const int n_beam_modes = 9;
int  nmodes, ptsrc_mode_index;
static const string ifn_ptsrc_mode="highl/clps.v4a4.dat",
                    ifn_beam_modes="highl/top_ten_modes.beam_covariance_VW_combined_7yr.dat",
                    ifn_fiducial_cltt="test_cls_v4.dat";
static const REAL_8 ptsrc_err = 0.1; // 10%
int lmin0, lmax0;
Matrix<REAL_8> mode, F_mode, beam_mode;
Array1D<REAL_8> fiducial_cltt, ptsrc_mode;
Array1D<REAL_8> a, a2, c;
Matrix<REAL_8> b;
Array1D<int> iwork;

// beam options See Appendix of Hinshaw, et.al. (2006) for a description.
static bool  beam_diagonal_sigma =  true ;
static bool  beam_gaussian_likelihood =  true ;
static bool  beam_fixed_fiducial_spectrum =  false ;
static bool  beam_include_beam_modes =  true ;
static bool  beam_include_ptsrc_mode =  true ;
}

void  init_tt_beam_and_ptsrc_chisq(int lmin, int lmax)
{


  using namespace  wmap_util;
  std::string WMAP_data_dir = WMAP_OPTIONS::self()->WMAP_data_dir();


  int   i, l;
  REAL_8  x;
  string  ifn;

  lmin0 = lmin;
  lmax0 = lmax;

  nmodes = 0;
  if ( beam_include_beam_modes )
  {
    nmodes = nmodes + n_beam_modes;
  }
  if ( beam_include_ptsrc_mode )
  {
    nmodes = nmodes + 1;
    ptsrc_mode_index = nmodes;
  }

  if ( nmodes == 0 )
  {
    return;
  }

  a.resize( Range( 1, nmodes ) );
  b.resize( Range( 1, nmodes ),Range( 1,nmodes ));
  c.resize(Range( 1,nmodes));
  a2.resize(Range( 1,nmodes));
  iwork.resize(Range( 1,nmodes));
  fiducial_cltt.resize(Range( lmin,lmax) );
  ptsrc_mode.resize(Range( lmin,lmax) );
  beam_mode.resize( Range(lmin,lmax ),Range( 1, n_beam_modes) );
  mode.resize(Range(lmin,lmax ),Range( 1,nmodes));
  F_mode.resize(Range( lmin,lmax ),Range( 1,nmodes) );

  if ( beam_include_beam_modes )
  {
    ifn = WMAP_data_dir+ifn_beam_modes;
    ifstream lun( ifn.c_str() );
    while ( lun )
    {
      lun >> i >> l >> x;
      if ( i <= n_beam_modes )
      {
        beam_mode(l,i) = x;
      }
    }
    lun.close();
  }

  if ( beam_fixed_fiducial_spectrum )
  {
    ifn = WMAP_data_dir+ifn_fiducial_cltt;
    ifstream lun( ifn.c_str() );
    while ( lun )
    {
      lun >> l >> x;
      if ( l >= lmin  &&  l <= lmax )
      {
        fiducial_cltt(l) = x;
      }
    }
    lun.close();
  }

  if ( beam_include_ptsrc_mode )
  {
    ifn = WMAP_data_dir+ifn_ptsrc_mode;
    ifstream lun( ifn.c_str() );
    while ( lun )
    {
      lun >> l >> x;
      if ( l >= lmin  &&  l <= lmax )
      {
        ptsrc_mode(l) = x;
      }
    }
    lun.close();
  }

}


void  quit_tt_beam_and_ptsrc_chisq( )
{
  if ( nmodes > 0 )
  {
    beam_mode.free();
    mode.free();
    F_mode.free();
    ptsrc_mode.free();
    fiducial_cltt.free();
    a.free();
    b.free();
    c.free();
    a2.free();
    iwork.free();
  }
}


REAL_8 compute_tt_beam_and_ptsrc_chisq( int lmin, int lmax, Array1D<REAL_8> &cltt,
                                        Array1D<REAL_8> &cltt_dat, Array1D<REAL_8> &neff,
                                        Matrix<REAL_8>& fisher, Array1D<REAL_8> &z,
                                        Array1D<REAL_8> &zbar )
{
  using namespace  wmap_util;

  REAL_8  compute_tt_beam_and_ptsrc_chisq_value;

  REAL_8  dgauss, dlnnorm, dlndet;
  int  i, j, l, l1, l2, stat;

  if ( nmodes == 0 )
  {
    return 0.0;
  }

  mode = 0.0;

  //! beam modes
  if ( beam_include_beam_modes )
  {
    if ( beam_fixed_fiducial_spectrum )
    {
      for (i = 1; i <= n_beam_modes; ++i )
      {
        for (l = lmin; l <= lmax; ++l )
        {
          mode(l,i) = beam_mode(l,i)*fiducial_cltt(l);
        }
      }
    }
    else
    {
      for (i = 1; i <= n_beam_modes; ++i )
      {
        for (l = lmin; l <= lmax; ++l )
        {
          mode(l,i) = beam_mode(l,i)*cltt(l);
        }
      }
    }
  }

  //! ptsrc mode
  if ( beam_include_ptsrc_mode )
  {
    //print *, 'including beam mode', ptsrc_mode_index
    for ( int l = lmin; l<=lmax; ++l )
      mode( l, ptsrc_mode_index ) = ptsrc_err*ptsrc_mode( l );
    //print *, 'ptsrc_mode(l=1000) = ', mode(1000,ptsrc_mode_index)
  }

  F_mode = 0.0;
  if ( beam_diagonal_sigma )
  {
    for (i = 1; i <= nmodes; ++i )
    {
      for (l = lmin; l <= lmax; ++l )
      {
        F_mode(l,i) = fisher(l,l)*mode(l,i);
      }
    }
  }
  else
  {
    for (i = 1; i <= nmodes; ++i )
    {
      for (l1 = lmin; l1 <= lmax; ++l1 )
      {
        for (l2 = lmin; l2 <= lmax; ++l2 )
        {
          //do l2 = l1-50,l1+50
          if ( l2 < lmin  ||  l2 > lmax )
            continue;
          if ( l2 < l1 )
          {
            F_mode(l1,i) = F_mode(l1,i) + fisher(l2,l1)*mode(l2,i);
          }
          else
          {
            F_mode(l1,i) = F_mode(l1,i) + fisher(l1,l2)*mode(l2,i);
          }
        }
      }
    }
  }

  a = 0.0;
  b = 0.0;
  for (i = 1; i <= nmodes; ++i )
  {
    for (l = lmin; l <= lmax; ++l )
    {
      a(i) = a(i) + (cltt_dat(l)-cltt(l))*F_mode(l,i);
    }
    b(i,i) = b(i,i) + 1.0;
    for (j = i; j <= nmodes; ++j )
    {
      for (l = lmin; l <= lmax; ++l )
      {
        b(i,j) = b(i,j) + mode(l,i)*F_mode(l,j);
      }
      if ( i != j )
        b(j,i) = b(i,j);
    }
  }

  //	print *, 'nmodes = ', nmodes
  //	do i = 1,nmodes
  //		print '("a(",I2,") = ",E)', i, a(i)
  //	end do

  gsl_matrix* gslb = gsl_matrix_alloc( nmodes, nmodes );
  for ( i = 0; i < nmodes; ++i )
    for ( j = 0; j < nmodes; ++j )
      gsl_matrix_set( gslb, i, j, b( i+1, j+1 ) );

  stat = gsl_linalg_cholesky_decomp( gslb );
  if ( stat != 0 )
  {
    gsl_matrix_free( gslb );
    wmap_likelihood_error( "beam/ptsrc: bad dpotrf", stat );
    //		print *, 'bad dpotrf'
    return 0.0;
  }

  c = a;
  gsl_vector* gslc = gsl_vector_alloc( nmodes );
  for ( i =0; i < nmodes; ++i )
    gsl_vector_set( gslc,  i, c( i+1 ) );
  stat = gsl_linalg_cholesky_svx( gslb, gslc );
  if ( stat != 0 )
  {
    gsl_vector_free( gslc );
    gsl_matrix_free( gslb );
    wmap_likelihood_error( "beam/ptsrc: bad dpotrs", stat );
    //         print *, 'bad dpotrs'
    return 0.0;
  }
  dgauss = 0.0;
  for (i = 1; i <= nmodes; ++i )
  {
    dgauss = dgauss + a(i) * gsl_vector_get( gslc, i-1 );
  }

  if ( beam_gaussian_likelihood )
  {
    dlndet = 1.0;
    for (i = 1; i <= nmodes; ++i )
    {
      double bii = gsl_matrix_get( gslb, i-1, i-1 );
      dlndet *= bii * bii;
    }
    dlndet = log(dlndet);

    //print *, 'beam chisq, lndet = ', dgauss, dlndet
    compute_tt_beam_and_ptsrc_chisq_value = -dgauss + dlndet;
  }
  else
  {
    a2 = 0.0;
    for (i = 1; i <= nmodes; ++i )
    {
      for (l = lmin; l <= lmax; ++l )
      {
        a2(i) = a2(i) + (z(l)-zbar(l))*(cltt(l)+neff(l))*F_mode(l,i);
      }
    }
    c = a2;
    for ( i =0; i < nmodes; ++i )
      gsl_vector_set( gslc,  i, c( i+1 ) );
    stat = gsl_linalg_cholesky_svx( gslb, gslc );
    if ( stat != 0 )
    {
      gsl_vector_free( gslc );
      gsl_matrix_free( gslb );
      wmap_likelihood_error( "beam/ptsrc: bad dpotrs", stat );
      //			print *, 'bad dpotrs'
      return 0.0;
    }
    dlnnorm = 0.0;
    for (i = 1; i <= nmodes; ++i )
    {
      dlnnorm += a2( i ) * gsl_vector_get( gslc, i-1 );
    }

    //print *, 'beam chisq, lndet = ', dgauss, dlnnorm, -(dgauss+2d0*dlnnorm)/3d0
    compute_tt_beam_and_ptsrc_chisq_value = -(dgauss + 2.0*dlnnorm)/3.0;
  }
  gsl_vector_free( gslc );
  gsl_matrix_free( gslb );
  return compute_tt_beam_and_ptsrc_chisq_value;
}
} // end namespace wmap_tt_beam_ptsrc_chisq


