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

// ============================================================================
// "Read_Archive_Map" reads a sky map from a binary table in a FITS file.  The
// map consists of two columns in the table; the first contains temperature
// and the second contains N_Obs.
//
// Arguments:
//   FileName  Character  The FITS filename.
//   Temp      Real*4     The array of temperatures.
//   NObs      Real*4     The array of observation numbers.
//   NPix      Integer    The number of elements in the output arrays.
//   Status    Integer    A status code: 0=success.
//
// Written by Michael R. Greason, SSAI, 07 February 2006.
// ============================================================================

#include "read_archive_map.h"

#include "port_helper.h"

#include "fitsio.h"

#include <iostream>
#include <string>
using namespace std;

void  Read_Archive_Map ( string FileName, real4_1d& Temp, real4_1d& NObs, long int& NPix, int& Status)
{
  int  i, typ;
  int anyf;
  // ----------------------------------------------------------------------------
  //           Open the FITS file.
  //
  NPix   = 0;
  Status = 0;
  i      = 0;
  fitsfile* fptr;
  fits_open_file( &fptr, FileName.c_str(), i, &Status );
  //
  //           Select the HDU.
  //
  i   = 1;
  typ = 0;
  while((typ  !=  2)  &&  (Status == 0))
  {
    i = i + 1;
    fits_movabs_hdu( fptr, i, &typ, &Status);
  }
  if ((typ  !=  2)  ||  (Status  !=  0))
  {
    cout << "Read_Archive_Map: There is no binary FITS table extension in " << FileName << ".* " << endl;
    return;
  }
  //
  //           Extract the number of pixels in the map.
  //
  long int tmpli;
  fits_get_num_rows( fptr, &tmpli, &Status );
  if (Status  !=  0)
    return;
  NPix = tmpli;
  //
  //           Extract the map.
  //
  fits_read_col_flt(fptr, 1, 1, 1, NPix, 0.0E0, Temp.data(), &anyf, &Status);
  fits_read_col_flt(fptr, 2, 1, 1, NPix, 0.0E0, NObs.data(), &anyf, &Status);
  //
  //           Close the FITS file.
  //
  fits_close_file( fptr, &Status);
  //
  return;
  // ----------------------------------------------------------------------------
  //
} //  Read_Archive_Map

