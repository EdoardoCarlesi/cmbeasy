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
#include "WMAP_7yr_tt_pixlike.h"

#include "WMAP_7yr_options.h"
#include "WMAP_7yr_util.h"


#include "read_archive_map.h"
#include "read_fits.h"

#include "gsl/gsl_linalg.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

// ===========================================================================
namespace  wmap_tlike
{
//
// This code calculates the likelihood function for TT data at low l 
//
//  1uK per pixel noise at nside=16 is added to regularize the matrix inversion.
//
//  It uses a smoothed ILC map as the temperature map and the Kp2 sky mask.
//
// O. Dore
// D. Spergel, December 2005
//
// << MODIFICATION HISTORY AFTER RELEASE on March 16, 2006 >>
//
// E. Komatsu, June 20, 2006
// -- Use Nside=16.
// -- Use the pre-computed smoothed and degraded ILC and V-band maps. 
//   HEALPIX ROUTINES ARE NO LONGER NEEDED.
// -- Smoothing scale has been increased to 9.1831 degrees.
//   I.e., the input map is band limited at l=40.
// -- White noise (1uK per pixel at nside=16) is added to the degraded maps.
// -- Fore_Marg is set to 2 explicitly.
//
//
// E. Komatsu, December 20, 2009
// -- 7-yr version
//===========================================================================

  int ngood, nsmax, nzzz;
  double_1d wl;
#ifdef FASTERTT
  double_1d  t_map; // ILC
#else
  real_1d  t_map; // ILC
  real_1d f_map; // difference between ILC and V
#endif
  double_2d  C0, C;  // C0 = monopole + dipole marginalized
  double_2d  zzz;
  double_2d  vecs;
  double_1d yyy;
  int_1d pix2x, pix2y; // for pix2vec
  static const int ifore_marg = 2;
  //
  // ifore_marg = 0 (no marginalization)
  // ifore_marg = 1  (assume that the uncertainty is the difference between
  //			V and ILC)
  // ifore_marg = 2  (no prior on amplitude of foreground unc.)
  //
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  namespace
  {
      void pix2vec_nest(unsigned int nside, unsigned int ipix, double_1d& vector, double_2d& vertex);
      void  mk_pix2xy();
      REAL randgauss_boxmuller(int& iseed /*random number state*/);
      REAL ran_mwc(int& iseed /*random number state*/);
  }


  void  setup_for_tt_exact(const int nlhigh, int* tt_ngood)
  {
      pix2x.resize(Range(0,1023));
      pix2y.resize(Range(0,1023));
      pix2x=0; pix2y=0;

    //
    // - READ IN DEGRADED Kp2 MASK (nside=16)
    //
    // - ILC AND V-BAND MAPS SMOOTHED WITH A BEAM, 
    //         exp(-0.5*l*(l+1)*theta_fwhm**2./(8.*alog(2.)))
    //   where
    //         theta_fwhm = 9.1285 degrees for ILC
    //         theta_fwhm = 9.1831 degrees for V-band 
    //                      [ sqrt(9.1285^2+1^2)=9.1831 ]
    //
    // - COMPUTES C0 for fiducial 
    //
    using namespace wmap_util;
    using namespace std;

    int  ires; //= 4           // resolution for low l TT code
    long int np;
    int  ip, jp, ipix, l, ll, jpix;
    string filename[7], ifn, ofn;
    bool  exist;
    int_1d  good;
    real_1d dummy,ilc,vband,mask;
//    REAL, ALLOCATABLE, dimension(:)  vband_clean
//X     REAL_3d  Plin;
    REAL  noise_rms,window_th,lc,theta_fwhm;  //=9.1831
    int  iseed, k, stat;
    double_1d  vec_ip(Range(1,3));
    double x, one_over_l;
    double_1d p(65);

#ifdef TIMING
  wmap_timing_start( "setup_for_tt_exact" );
#endif

	string WMAP_data_dir = WMAP_OPTIONS::self()->WMAP_data_dir();
	ires = WMAP_OPTIONS::self()->lowl_tt_res();
	
        if ( nlhigh != 30 )  {
		cout << "error: set lowl_max=30 when using_gibbs=false" << endl;
		exit(-1);
	}

#ifdef FASTERTT
	if ( ifore_marg != 2 )  {
		cout <<  "*** FASTERTT currently requires ifore_marg=2" << endl ;
		exit(-1);
	}

	ngood = 957;

	C.resize(Range(0,ngood-1),Range(0,ngood-1));
	t_map.resize(Range(0,ngood-1));

	ifn = WMAP_data_dir+"/lowlT/faster/compressed_t_map_f2_7yr.txt";
	ifstream some_lun(ifn.c_str());
	for (ip = 0; ip <= ngood-1; ++ip ) {
		some_lun >> t_map(ip);
	}
        some_lun.close();

	nzzz = (ngood*(ngood+1))/2;
	zzz.resize(Range(1,nzzz),Range(1,nlhigh-1));
	yyy.resize(Range(1,nzzz));

//	do l = 2,nlhigh
//
//		write(ifn,'("/cita/d/raid-sievers2/sievers/WMAP3_cosmomc/faster_TT/wmap_compressed_fproj_",I0.2,".unf")') l
//		open(lun,file=ifn,status='old',action='read',form='unformatted')
//		read(lun) C
//		close(lun)
//
//		k = 1
//		do ip = 0,ngood-1
//		do jp = ip,ngood-1
//			zzz(k,l-1) = C(ip,jp)
//			k = k + 1
//		end do
//		end do
//	end do
//	C = 0d0
//
//	ofn = "compressed_tt_f2.unf"
//	open(lun,file=ofn,action="write",status="unknown",form="unformatted")
//	write(lun) zzz
//	close(lun)

	ifn = WMAP_data_dir+"lowlT/faster/compressed_tt_matrices_f2_7yr.fits";
	Read_FITS_Double_2D(ifn, zzz, stat);
	if (stat  !=  0)  {
		cout <<  "Error " << stat << " while reading " << ifn << endl ;
		exit(-1);
	}

#else

	//filename(0)=trim(WMAP_data_dir)//'wmap_lcdm_pl_model_yr1_v1.dat'
	filename[0]=WMAP_data_dir+"test_cls_v4.dat";

	switch(ires){
// 	case 3:
// 		theta_fwhm=18.3;
// 		filename[1]=WMAP_data_dir+"mask/mask_r3_kp2.fits";
// 		filename[2]=WMAP_data_dir+"maps/low_resolution_map_fs_r3_ilc_smooth_18.3deg.fits"; // ILC map
// 		filename[3]=WMAP_data_dir+"maps/low_resolution_map_fs_r3_vband_smooth_18.3deg.fits"; // raw V-band map for foreground marginalization
// 		filename[4]=WMAP_data_dir+"healpix_data/pixel_window_n0008.txt";
// 		//filename(5)=trim(WMAP_data_dir)//'pl_res3_sp.fits'
// 		//filename(6)=trim(WMAP_data_dir)//'maps/low_resolution_map_fs_r3_vband_clean_smooth_18.3deg.fits' ! template-cleaned V-band map -- not used
//                 break;
	case 4:
		theta_fwhm=9.1831;

		filename[1]=WMAP_data_dir+"lowlT/pixlike/mask_r4_kq_85p0.fits";
		filename[2]=WMAP_data_dir+"lowlT/pixlike/wmap_7yr_r4_ilc_smooth_9.1285deg.fits";
		filename[3]=WMAP_data_dir+"lowlT/pixlike/wmap_7yr_r4_vband_raw_smooth_9.1831deg.fits";

		filename[4]=WMAP_data_dir+"healpix_data/pixel_window_n0016.txt";
		//filename(5)=trim(WMAP_data_dir)//'' ! used only when ires=3
		//filename(6)=trim(WMAP_data_dir)//'maps/wmap_fs_r4_vband_clean_smooth_9.1831deg.fits'! smoothed cleaned V-band -- not used
                break;
	default:
		wmap_likelihood_error( "bad lowl_tt_res", ires );
		exit(-1);
        }
	lc = sqrt(8.*log(2.))/theta_fwhm/3.14159*180.;

    //
    //       READ IN cl's for fiducial model (used for high l)
    //
    double_1d  cl_fiducial;
    cl_fiducial.resize(Range(2,512));
    ifstream lun(filename[0].c_str());
    if ( !lun.good() )  {
       cout << " unable to open " << filename[0] << endl;
       exit(-1);
    }
    double dummy1, dummy2, dummy3;
    for (l = 2; l <= 512; ++l ) {
       lun >> ll >> cl_fiducial(l) >> dummy1 >> dummy2 >> dummy3;
       cl_fiducial(ll) = cl_fiducial(ll)*1.e-6*2.*3.1415927/(ll*(ll+1));
    }
    lun.close();
    //
    // READ IN LOW-RESOLUTION Kp2 MASK
    //
    nsmax = (int)pow(2.,ires);
    nlmax = 4*nsmax;
    npix = 12*nsmax*nsmax;
    dummy.resize(Range(0,npix-1));
    mask.resize(Range(0,npix-1));
    Read_Archive_Map(filename[1],dummy,mask,np,status);
    if (status != 0 || (np != npix))  {
       cout << " error reading mask" << endl;
       exit(-1);
    }
    ngood = (int)mask.sum();
    good.resize(Range(0,ngood-1));
    ip = 0;
    for (ipix = 0; ipix <= 12*nsmax*nsmax-1; ++ipix ) {
       if (mask(ipix) == 1)  {
          good(ip) = ipix;
          ip = ip+1;
       }
    }
    if (ip != ngood)  exit(-1);
    //
    //  READ IN LOW-RESOLUTION ILC AND V-BAND MAPS
    //  - smoothing scale is FWHM = 9.1 degrees
    //
    ilc.resize(Range(0,npix-1)); vband.resize(Range(0,npix-1));
//    ALLOCATE(vband_clean(0:npix-1))
    Read_Archive_Map(filename[2],ilc,dummy,np,status);
    if (status != 0 || (np != npix))  {
       cout << " error reading ILC map" << endl;
       exit(-1);
    }
    Read_Archive_Map(filename[3],vband,dummy,np,status);
    if (status != 0 || (np != npix))  {
       cout << " error reading V-band map" << endl;
       exit(-1);
    }
//    call read_archive_map(filename(6),vband_clean,dummy,np,status)
//    if(status.ne.0.or.np.ne.npix) then
//       write(*,*) ' error reading cleaned V-band map'
//       stop
//    endif
    dummy.free();

    //
    //  CALCULATE COMBINED BEAM
    //
    wl.resize(Range(2,1024));
    ifstream lun2(filename[4].c_str());
    lun2 >> window_th;
    lun2 >> window_th;
    for (l = 2; l <= nlmax; ++l ) {
       lun2 >> window_th;
       wl(l) = window_th*exp(-0.5*l*(l+1)/(lc*lc));
    }
    lun2.close();
    double_1d cl(cl_fiducial(Range(2,nlmax)));
    //FIXME: c++ port: not the fastest way to do this...
    cl *= wl(Range(2,nlmax));
    cl *= wl(Range(2,nlmax));
    cl_fiducial.free();

    C0.resize(Range(0,ngood-1),Range(0,ngood-1));
    C.resize(Range(0,ngood-1),Range(0,ngood-1));
    t_map.resize(Range(0,ngood-1));
    f_map.resize(Range(0,ngood-1));
    vecs.resize(Range(1,3),Range(0,ngood-1));

    C0 = 0.0;

    for (ip = 0; ip <= ngood-1; ++ip ) {
       t_map(ip) = ilc(good(ip));
//!$ t_map(ip) = vband_clean(good(ip)) !!$ don't use cleaned V. use ILC!
       f_map(ip) = vband(good(ip))-ilc(good(ip)); // foreground template = V-band - ilc
      
	//! MRN
        double_2d dummy_vertex;
	pix2vec_nest( nsmax, good(ip), vec_ip, dummy_vertex );
	for (int idx=1; idx <=3; ++idx){
           vecs(idx,ip) = vec_ip(idx);
        }
    }

       //! MRN
       for (ip = 0; ip <= ngood-1; ++ip ) {
          for (jp = ip; jp <= ngood-1; ++jp ) {
              x = 0;
	        for (int idx=1; idx <=3; ++idx)
		    x += vecs(idx,ip)*vecs(idx,jp);
                p(0) = 1.0;
                p(1) = x;
                for (l = 2; l <= nlmax; ++l ) {
			one_over_l = 1.0/l;
                        p(l) = (2.0-one_over_l)*x*p(l-1)-(1.0-one_over_l)*p(l-2);
                }
                for (l = 0; l <= nlmax; ++l ) {
                        p(l) = p(l)*(2.0*l+1.0)/(4.0*3.1415927);
                }
                double_1d sum_tmp = p(Range(nlhigh+1,nlmax));
                sum_tmp *= cl(Range(nlhigh+1,nlmax));
          	C0(ip,jp) = sum_tmp.sum()
			+ cl(2)*(p(0)+p(1));
          }
       }

    good.free(); ilc.free(); vband.free(); cl.free();
    //DEALLOCATE(vband_clean)

    //
    // Add random noise (1uK per pixel) to regularize the covariance matrix
    // The map will be noise-dominated at l>40 (l>20 at res3).
    //
    iseed = -562378; // arbitrary
    noise_rms = 1.e-3; // mK per pixel

    for (ip = 0; ip <= ngood-1; ++ip ) {
       C0(ip,ip) = C0(ip,ip) + noise_rms*noise_rms;
       t_map(ip) = t_map(ip) + noise_rms*randgauss_boxmuller(iseed);
    }
#endif
	if ( tt_ngood != 0 )  {
		*tt_ngood = ngood;
	}
 
#ifdef TIMING
	wmap_timing_end();
#endif

  } //  SETUP_FOR_TT_EXACT


unsigned int tt_pixlike_dof()
{
	return ngood;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  void compute_tt_pixlike(int nlhigh,double_1d& cl_in, double& chisq, double& lndet)
  {

    //
    //  This subroutine returns the likelihood 
    //
    using namespace  wmap_util;
    const double tt_pixlike_lndet_offset = WMAP_OPTIONS::self()->tt_pixlike_lndet_offset();
    const double lowl_tt_res = WMAP_OPTIONS::self()->lowl_tt_res();

    double_1d cl;
    double like;
    REAL  fore_marg,fCinvf;
    double  DDOT;
    int  info,ip,jp, l, k, i;
    double_1d vec_ip(Range(1,3)), p(Range(0,64));
    double  one_over_l, x;
    double  dk, dindex, wt, tot_cl, xxx1, xxx2, max_diff, xxx;
    double  omega_pix;

#ifdef TIMING
    wmap_timing_start( "compute_tt_pixlike" );
#endif

    cl.resize(Range(2,nlhigh));

#ifdef FASTERTT

    cl = cl_in(Range(2,nlhigh));

    //X 	call DGEMV( "N", nzzz, nlhigh-1, 1d0, zzz, 
    //X 		nzzz, cl(2:nlhigh), 1, 0d0, yyy, 1 );
        yyy = 0.;
        //zzz starts from (0,0) as first index, yyy from 1
        //also, we have read in zzz the wrong way around compared to
        //fortran, so we do the indexing and dgemv by hand
        for ( int i = 0; i < nzzz; ++i)
          for ( int j = 0; j < nlhigh-1; ++j)
              //yyy(i+1) += zzz(i,j)*cl(2+j);
              yyy(i+1) += *(zzz.data()+(j*(nzzz)+i))*cl(2+j);

        C = 0.0;
        for (ip = 0; ip <= ngood-1; ++ip ) {
          C(ip,ip) = 1.0;
        }

        k = 1;
        // **this part should not be modified because of the way "yyy" is ordered (EK)**  (also true for c++ port?)
        for (ip = 0; ip <= ngood-1; ++ip ) {
          for (jp = ip; jp <= ngood-1; ++jp ) {
            C(ip,jp) = C(ip,jp) + yyy(k);
            k = k + 1;
          }
        }

#else

    for (l =2; l <= nlhigh; ++l ) {
       cl(l) = cl_in(l)*wl(l)*wl(l)*0.5e-6*(2.*l+1.)/(l*(l+1.));
    }

    // need to fill in only the upper half of C

    C = C0;

       //! MRN
       for (ip = 0; ip <= ngood-1; ++ip ) {
          for (jp = ip; jp <= ngood-1; ++jp ) {
                x = 0;
                for (int idx = 1; idx <= 3; ++idx) {
		  x += vecs(idx,ip)*vecs(idx,jp);
                }
                p(0) = 1.0;
                p(1) = x;
                for (l = 2; l <= nlhigh; ++l ) {
			one_over_l = 1.0/l;
                        p(l) = (2.0-one_over_l)*x*p(l-1)-(1.0-one_over_l)*p(l-2);
                }
                double_1d cl_sum_tmp(cl(Range(2, nlhigh)));
                cl_sum_tmp *= p(Range(2, nlhigh));
          	C(ip,jp) = C(ip,jp) + cl_sum_tmp.sum();
          }
       }

#endif

    cl.free();

#ifdef TIMING
	wmap_timing_checkpoint( "finished C" );
#endif
//X     DPOTRF
//X     CALL DPOTRF("U",ngood,C,ngood,info) ;

    gsl_matrix* gslC = gsl_matrix_alloc( ngood, ngood );
    for (ip=0; ip<ngood; ++ip)
      for (jp=0; jp<ngood; ++jp)
        gsl_matrix_set( gslC, ip, jp, C( ip, jp ) );

    gsl_matrix_transpose( gslC );
    info = gsl_linalg_cholesky_decomp( gslC );

    if (info != 0)  {
       wmap_likelihood_error( "tlike 1st dpotrf failed", info );
       chisq = 0.0;
       lndet = 0.0;
       return;
    }

#ifdef TIMING
	wmap_timing_checkpoint( "finished dpotrf" );
#endif

//    lndet = 0d0
//    omega_pix = 4d0*3.145927d0/dble(12*nsmax**2)
//    do ip=0,ngood-1
//       lndet = lndet + 2.*log(C(ip,ip)*omega_pix)
//    end do

    lndet = 0.0;
    for (ip=0; ip<= ngood-1; ++ip) {
       lndet = lndet + 2.*log(gsl_matrix_get( gslC, ip,ip));
    }

#ifdef TIMING
	wmap_timing_checkpoint( "finished lndet" );
#endif

//X     call DPOTRS("U", ngood, 1, C, ngood, c_inv_t_map, ngood, INFO);
    gsl_vector * c_inv_t_map =  gsl_vector_alloc( ngood );
    for ( int i = 0; i < ngood; ++i )
      gsl_vector_set( c_inv_t_map, i, t_map( i ) );

    info = gsl_linalg_cholesky_svx( gslC, c_inv_t_map );
    if (info != 0)  {
       wmap_likelihood_error( "tlike 2nd dpotrs failed", info );
       chisq = 0.0;
       lndet = 0.0;
       return;
    }

//X     call DPOTRS("U", ngood, 1, C, ngood, c_inv_f_map, ngood, INFO);
#ifndef FASTERTT
    gsl_vector * c_inv_f_map =  gsl_vector_alloc( ngood );
    for ( int i = 0; i < ngood; ++i )
      gsl_vector_set( c_inv_f_map, i, f_map( i ) );

    info = gsl_linalg_cholesky_svx( gslC, c_inv_f_map );
    if (info != 0)  {
       wmap_likelihood_error( "tlike 3rd spotrs failed", info );
       chisq = 0.0;
       lndet = 0.0;
       return;
    }

    fCinvf = 0;
    for ( int i = 0; i < ngood; ++i )
        fCinvf += f_map(i)*gsl_vector_get( c_inv_f_map, i )/ngood;
    //DEALLOCATE(C)
#endif

#ifdef TIMING
    wmap_timing_checkpoint( "finished dpotrs" );
#endif

    chisq = 0;
    for ( int i = 0; i < ngood; ++i )
         chisq += t_map(i)*gsl_vector_get( c_inv_t_map, i );

#ifndef FASTERTT
    if (ifore_marg != 0)  {
        if (ifore_marg == 1)  {
          fore_marg = 1.;
       } else {
          fore_marg = 1.e6;
       }
       double tmpSum = 0;
       for ( int i = 0; i < ngood; ++i )
         tmpSum += t_map(i)*gsl_vector_get( c_inv_f_map, i );
       chisq = chisq - (tmpSum*tmpSum)/ngood
            /( 1./fore_marg + fCinvf );
       lndet = lndet + log(1./fore_marg+fCinvf);
    }
    gsl_vector_free(c_inv_f_map);
#endif
    like = (chisq+lndet)/2.0;
    chisq=chisq/2.;
    lndet=(lndet-tt_pixlike_lndet_offset)/2.;

    gsl_matrix_free(gslC);
    gsl_vector_free(c_inv_t_map);

#ifdef TIMING
    wmap_timing_end();
#endif
} //  COMPUTE_TT_PIXLIKE

namespace
{
//
// The following two subroutines are taken from ran_tools in HEALPix.
//
  //=======================================================================
REAL randgauss_boxmuller(int& iseed /*random number state*/)
{
    static bool empty = true;
    REAL randgauss_boxmuller_retval;
    REAL  fac,rsq,v1,v2;
    static REAL gset; //// test
    
    if (empty  ||  iseed < 0)  { // bug correction, EH, March 13, 2003
       while (true) {
          v1 = 2.*ran_mwc(iseed) - 1.;
          v2 = 2.*ran_mwc(iseed) - 1.;
          rsq = v1*v1+v2*v2;
          if ((rsq<1.)  &&  (rsq>0.)) break;
       }

       fac = sqrt(-2.*log(rsq)/rsq);
       gset = v1*fac;
       randgauss_boxmuller_retval = v2*fac;
       empty = false;
    } else {
       randgauss_boxmuller_retval = gset;
       empty = true ;
    }
   return randgauss_boxmuller_retval;
}

//=======================================================================
REAL ran_mwc(int& iseed /*random number state*/)
{
    REAL  ran_mwc_retval;   //// result
    int  i,iseedl,iseedu,mwc,combined;
    static unsigned int upper; static unsigned int lower; static unsigned int shifter;
    const int mask16=65535; const int mask30=2147483647;
    static REAL small;
    static bool first = true;
 
    if (first || (iseed<=0))  {
       if (iseed==0)  iseed=-1;
       iseed = abs(iseed);
       //small = nearest (1.,-1.)/mask30;
       //small = 4.6566123e-10;
       small = nextafter(1.,-1.)/mask30;

       // Users often enter small seeds - I spread them out using the
       // Marsaglia shifter a few times.
       shifter=iseed;
       for (i=1; i<= 9; ++i) {
          shifter ^= (shifter<<13);
          shifter ^= (shifter>>17);
          shifter ^= (shifter<<5);
       }

       iseedu=(shifter>>16);
       upper=((iseedu+8765)<<16)+iseedu; //This avoids the fixed points.
       iseedl=(shifter&mask16);
       lower=((iseedl+4321)<<16)+iseedl; //This avoids the fixed points.

       first= false;
    }

    do {
       shifter ^= (shifter<<13);
       shifter ^= (shifter>>17);
       shifter ^= (shifter<<5);

       upper=36969*(upper&mask16)+(upper>>16);
       lower=18000*(lower&mask16)+(lower>>16);

       mwc=(upper<<16)+(lower&mask16);

       combined=(mwc&mask30)+(shifter&mask30);

       ran_mwc_retval=small*(combined&mask30);
       if (ran_mwc_retval!=0.)  break;
    }while (true);
    return ran_mwc_retval;
} //  end function ran_mwc;

//
// The following two subroutines are taken from pix_tools in HEALPix.
//
  //=======================================================================
 void  pix2vec_nest(unsigned int nside, unsigned int ipix, double_1d& vector, double_2d& vertex)
  {
    using namespace  wmap_util;

    //=======================================================================
    //     renders vector (x,y,z) coordinates of the nominal pixel center
    //     for the pixel number ipix (NESTED scheme)
    //     given the map resolution parameter nside
    //     also returns the (x,y,z) position of the 4 pixel vertices (=corners)
    //     in the order N,W,S,E
    //=======================================================================
//X     double, INTENT(OUT), dimension(1:)  vector;
//X     double,     INTENT(OUT),dimension(1:,1:), optional  vertex;
    unsigned int oldVectorBaseIndex = vector.firstIndex();
    unsigned int oldVertexFirstRow = vertex.rowsFirstIndex();
    unsigned int oldVertexFirstCol = vertex.colsFirstIndex();
    vector.reindexSelf(1);
    vertex.reindexSelf(1, 1);

    int  npix, npface,
              ipf, ip_low, ip_trunc, ip_med, ip_hi,
              jrt, jr, nr, jpt, jp, kshift, nl4;
    double  z, fn, fact1, fact2, sth, phi;

    int   ix, iy, face_num;
//     common /xy_nest/ ix, iy, face_num ! can be useful to calling routine

    // coordinate of the lowest corner of each face
    int jrll_init[12] = {2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4};  // in unit of nside
    int jpll_init[12] = {1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 }; // in unit of nside/2
    int_1d jrll(Range(1,12));  jrll = jrll_init;
    int_1d jpll(Range(1,12));  jpll = jpll_init;

    double  phi_nv, phi_wv, phi_sv, phi_ev, phi_up, phi_dn;
    double  z_nv, z_sv, sth_nv, sth_sv;
    double  hdelta_phi;
    int  iphi_mod, iphi_rat;
    bool  do_vertex;
    const double  halfpi = 1.570796326794896619231321691639751442099;
    const double  pi = 3.141592653589793238462643383279502884197;
    int  info = -1;
    //-----------------------------------------------------------------------
    if (nside<1  ||  nside>1024) wmap_likelihood_error("nside out of range",info);
    npix = 12 * nside*nside;
    if (ipix <0  ||  ipix>npix-1) wmap_likelihood_error("ipix out of range",info);

    //     initiates the array for the pixel number -> (x,y) mapping
    if (pix2x(1023) <= 0) mk_pix2xy();

    fn = double(nside);
    fact1 = 1.0/(3.0*fn*fn);
    fact2 = 2.0/(3.0*fn);
    nl4   = 4*nside;

    do_vertex =  false;
    if (vertex.size() > 0)  {
       if (vertex.rows() >= 3  &&  vertex.cols() >= 4)  {
          do_vertex =  true;
       } else {
          wmap_likelihood_error(" pix2vec_ring : vertex array has wrong size ",info);
       }
    }

    //     finds the face, and the number in the face
    npface = nside*nside;

    face_num = ipix/npface;  // face number in {0,11}
    ipf = ipix % npface;  // pixel number in the face {0,npface-1}

    //     finds the x,y on the face (starting from the lowest corner)
    //     from the pixel number
    ip_low = ipf % 1024;       // content of the last 10 bits
    ip_trunc =   ipf/1024;        // truncation of the last 10 bits
    ip_med = ip_trunc % 1024;  // content of the next 10 bits
    ip_hi  =     ip_trunc/1024;   // content of the high weight 10 bits

    ix = 1024*pix2x(ip_hi) + 32*pix2x(ip_med) + pix2x(ip_low);
    iy = 1024*pix2y(ip_hi) + 32*pix2y(ip_med) + pix2y(ip_low);

    //     transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy;  // "vertical" in {0,2*(nside-1)}
    jpt = ix - iy;  // "horizontal" in {-nside+1,nside-1}

    //     computes the z coordinate on the sphere
    jr =  jrll(face_num+1)*nside - jrt - 1;   // ring number in {1,4*nside-1}

    nr = nside;                  // equatorial region (the most frequent)
    z  = (2.*nside-(double)jr)*fact2;
    kshift = (jr - nside) % 2;
    if (do_vertex)  {
       z_nv = (2*nside-jr+1)*fact2;
       z_sv = (2*nside-jr-1)*fact2;
       if (jr == nside)  { // northern transition
          z_nv =  1.0 - (nside-1)*(nside-1) * fact1;
       } else if (jr == 3*nside)  {  // southern transition
          z_sv = -1.0 + (nside-1)*(nside-1) * fact1;
       }
    }
    if (jr < nside)  {     // north pole region
       nr = jr;
       z = 1.0 - nr*nr*fact1;
       kshift = 0;
       if (do_vertex)  {
          z_nv = 1.0 - (nr-1)*(nr-1)*fact1;
          z_sv = 1.0 - (nr+1)*(nr+1)*fact1;
       }
    } else if (jr > 3*nside)  { // south pole region
       nr = nl4 - jr;
       z = - 1.0 + nr*nr*fact1;
       kshift = 0;
       if (do_vertex)  {
          z_nv = - 1.0 + (nr+1)*(nr+1)*fact1;
          z_sv = - 1.0 + (nr-1)*(nr-1)*fact1;
       }
    }

    //     computes the phi coordinate on the sphere, in [0,2Pi]
    jp = (jpll(face_num+1)*nr + jpt + 1 + kshift)/2;  // "phi" number in the ring in {1,4*nr}
    if (jp > nl4)  jp = jp - nl4;
    if (jp < 1)    jp = jp + nl4;

    phi = (jp - (kshift+1.)*0.5) * (halfpi / nr);

    sth = sqrt((1.0-z)*(1.0+z));
    vector(1) = sth * cos(phi);
    vector(2) = sth * sin(phi);
    vector(3) = z;

    if (do_vertex)  {
       phi_nv = phi;
       phi_sv = phi;

       phi_up = 0.0;
       iphi_mod = (jp-1) % nr; // in {0,1,... nr-1}
       iphi_rat = (jp-1) / nr;      // in {0,1,2,3}
       if (nr > 1) phi_up = halfpi * (iphi_rat +  iphi_mod   /double(nr-1));
       phi_dn             = halfpi * (iphi_rat + (iphi_mod+1)/double(nr+1));
       if (jr < nside)  {            // North polar cap
          phi_nv = phi_up;
          phi_sv = phi_dn;
       } else if (jr > 3*nside)  {     // South polar cap
          phi_nv = phi_dn;
          phi_sv = phi_up;
       } else if (jr == nside)  {      // North transition
          phi_nv = phi_up;
       } else if (jr == 3*nside)  {    // South transition
          phi_sv = phi_up;
       }

       hdelta_phi = pi / (4.0*nr);

       // west vertex
       phi_wv      = phi - hdelta_phi;
       vertex(1,2) = sth * cos(phi_wv);
       vertex(2,2) = sth * cos(phi_wv);
       vertex(3,2) = z;

       // east vertex
       phi_ev      = phi + hdelta_phi;
       vertex(1,4) = sth * cos(phi_ev);
       vertex(2,4) = sth * sin(phi_ev);
       vertex(3,4) = z;

       // north vertex
       sth_nv = sqrt((1.0-z_nv)*(1.0+z_nv));
       vertex(1,1) = sth_nv * cos(phi_nv);
       vertex(2,1) = sth_nv * sin(phi_nv);
       vertex(3,1) = z_nv;

       // south vertex
       sth_sv = sqrt((1.0-z_sv)*(1.0+z_sv));
       vertex(1,3) = sth_sv * sqrt(phi_sv);
       vertex(2,3) = sth_sv * sqrt(phi_sv);
       vertex(3,3) = z_sv;
    }

    vector.reindexSelf(oldVectorBaseIndex);
    vertex.reindexSelf(oldVertexFirstRow, oldVertexFirstCol);

    return;
  } //  pix2vec_nest

  //=======================================================================
  void  mk_pix2xy()
  {

    //=======================================================================
    //     constructs the array giving x and y in the face from pixel number
    //     for the nested (quad-cube like) ordering of pixels
    //
    //     the bits corresponding to x and y are interleaved in the pixel number
    //     one breaks up the pixel number by even and odd bits
    //=======================================================================
    int   kpix, jpix, IX, IY, IP, ID;

    //cc cf block data      data      pix2x(1023) /0/
    //-----------------------------------------------------------------------
    //      print *, 'initiate pix2xy'
    for (kpix=0 /*pixel number*/; kpix<= 1023; ++kpix) {
       jpix = kpix;
       IX = 0;
       IY = 0;
       IP = 1;               // bit position (in x and y)
//        do while (jpix/=0) ! go through all the bits
       do {
          if (jpix == 0) break; // go through all the bits
          ID = jpix % 2;  // bit value (in kpix), goes in ix
          jpix = jpix/2;
          IX = ID*IP+IX;

          ID = jpix%2;  // bit value (in kpix), goes in iy
          jpix = jpix/2;
          IY = ID*IP+IY;

          IP = 2*IP;         // next bit (in x and y)
       } while (true);
       pix2x(kpix) = IX;     // in 0,31
       pix2y(kpix) = IY;     // in 0,31
    }
  } //  mk_pix2xy
} //anonymous namespace
} // end namespace  wmap_tlike
