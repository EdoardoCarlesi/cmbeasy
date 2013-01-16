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

#include "br_mod_dist.h"

#include "fitsio.h"

#include <iostream>
#include <cmath>

namespace  br_mod_dist {
  namespace {
    int lmin, lmax, numsamples, numchain;
    bool first_eval;
    double offset;
    double_3d sigmas;
  }

  // Initialization routines
  void  initialize_br_mod(const int lmin_in, const double_3d& sigmas_in)
  {
//X     int(i4b),                            intent(in)   lmin_in;

//X     REAL(sp),     dimension(lmin_in:,1:,1:), intent(in)   sigmas_in;

    int i, j, ell;
    lmin       = lmin_in;
    lmax       = sigmas_in.dim1Size() + lmin - 1;
    numchain   = sigmas_in.dim2Size();
    numsamples = sigmas_in.dim3Size();
    first_eval =  true ;

    sigmas = sigmas_in;

    for (ell = lmin; ell <= lmax; ++ell ) {
       for (i = 1; i <= numchain; ++i ) {
          for (j = 1; j <= numsamples; ++j ) {
             if (sigmas(ell, i, j)  <=  0.0)  {
                using namespace std;
                cout <<  "Error: sigma value <= zero." << endl ;
                cout <<  "sigma value is " << sigmas(ell, i, j) << endl ;
                cout <<  "at (ell,chain,sample) = " << ell << " " << i << " " << j << endl ;
                exit(-1);
             }
          }
       }
    }

  } //  initialize_br_mod



  void  clean_up_br_mod()
  {
    sigmas.free();
  } //  clean_up_br_mod




  // Base computation routine
  void  compute_br_estimator(const double_1d& cls, double& lnL)
  {
//X     REAL(dp), dimension(lmin:), intent(in)   cls;
//X     REAL(dp),                   intent(out)  lnL;

    int i, j, l;
    double      subtotal, x;

    if (first_eval)  {
       compute_largest_term(cls);
       first_eval =  false ;
    }

    // Compute the Blackwell-Rao estimator
    lnL = 0.0;
    for (i = 1; i <= numchain; ++i ) {
       for (j = 1; j <= numsamples; ++j ) {

          subtotal = 0.0;
          for (l = lmin; l <= lmax; ++l ) {
             x = sigmas(l,i,j)/cls(l);
             subtotal = subtotal +
                  0.50 * double(2*l+1) * (-x + log(x)) - log(double(sigmas(l,i,j)));
          }

          lnL = lnL + exp(subtotal-offset);

       }
    }

    if (lnL > 1e-20)  {
       lnL = log(lnL);
    } else {
       lnL = log(1e-30);
    }

    // print *, lnL

  } //  compute_br_estimator




  // Routine for reading the Gibbs sigma samples 
  void  read_gibbs_chain(const std::string& filename, long int& lmax, long int& numchains,
                                                      long int& numsamples, FloatToDoubleArray3D& data)
  {
    int          l, status, blocksize, readwrite, i, j, k;
    long int     numspec;
    int          fpixel, group, numargs, anyf;
    float        nullval;
    char*  comment;
    fitsfile*    unit;

    int          naxes[4];
//X     REAL(sp),     pointer, dimension(:,:,:,:)  indata;
    real_3d indata;

    status = 0;
    readwrite = 0;
    nullval = 0.;

    // numargs = 1
    numargs = 0;

    // Open the result file
    fits_open_file(&unit,filename.c_str(),readwrite,&status);

    // Read keywords
    comment = NULL;  //X we don't care about the comment - gr
    fits_read_key_lng(unit,"LMAX",     &lmax,       comment, &status);
    fits_read_key_lng(unit,"NUMSAMP",  &numsamples, comment, &status);
    fits_read_key_lng(unit,"NUMCHAIN", &numchains,  comment, &status);
    fits_read_key_lng(unit,"NUMSPEC",  &numspec,    comment, &status);

    //data.resize(Range(0,lmax), Range(1,numchains), Range(1,numsamples));
    indata.resize(Range(0,lmax), Range(1,numchains), Range(1,numsamples));
    indata = 50.;

//!$    print *, "Allocated arrays"

    // Read the binned power spectrum array
    group  = 1;
    fpixel = 1;
    fits_read_img_flt(unit, group, fpixel, indata.size(), nullval,
                                           indata.data(), &anyf, &status);

//!$    print *, "Read data"

    fits_close_file(unit,&status);

//!$    print *, "Closed file"

    /* no need to port this - gr
    for (i = 0; i <= lmax; ++i ) {
       for (j = 1; j <= numchains; ++j ) {
          for (k = numargs+1; k <= numargs+numsamples; ++k ) {
             //Xdata(i, j, k) = indata(i, 1, j, k);
             // data(:,:,:) = indata(0:lmax,1:1,1:numchains,numargs+1:numargs+numsamples)
          }
       }
    }
    */
    data.fromFloat(indata);

//!$    print *, "Deallocating data"

    indata.free();

//!$    print *, "Leaving subroutine"

  } //  read_gibbs_chain




  // Utility routine for initializing the offset to be subtracted from each term
  // to avoid overflow errors. Only called with the first power spectrum
  void  compute_largest_term(const double_1d& cls)
  {
//X     REAL(dp), dimension(lmin:), intent(in)   cls;

    int  i, j, l;
    double subtotal, x;

    // Compute the Blackwell-Rao estimator
    offset = -1.6375e30;
    for (i = 1; i <= numchain; ++i ) {
       for (j = 1; j <= numsamples; ++j ) {

          subtotal = 0.0;
          for (l = lmin; l <= lmax; ++l ) {
             x = sigmas(l,i,j)/cls(l);
             subtotal = subtotal +
                   0.50 * double(2*l+1) * (-x + log(x)) - log(double(sigmas(l,i,j)));
          }
          offset = std::max(offset,subtotal);
       }
    }

    if (offset < -1.637e30)  {
       using namespace std;
       cout <<  "Error: offset in br_mod_dist not being computed properly" << endl ;
       cout <<  "lmin = " << lmin << endl ;
       cout <<  "lmax = " << lmax << endl ;
       cout <<  "numchain = " << numchain << endl ;
       cout <<  "numsamples = " << numsamples << endl ;
       cout <<  "offset = " << offset << endl ;
       //cout <<  "cls = " << cls(lmin:lmax) << endl ;
       //cout <<  "sigmas(lmin:lmax, 10, 10) = ", sigmas(lmin:lmax, 10, 10) << endl ;
       exit(-1);
    }

  } //  compute_largest_term



} // end namespace  br_mod_dist
