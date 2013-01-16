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
// 'Read_FITS_Real_2D' reads a two-dimensional floating point array from the
// primary header/data unit in a FITS file.  This routine allocates the array
// from the heap if it hasn't already been associated.
//
// If the array has already been created, this routine assumes that it has
// the correct size and dimensions!  If allocated here then it will be indexed
// like a C array:  0:n-1.
//
// Arguments:
//	FileName - The name of the FITS file.
//	Arr      - The array to store the data in.
//	Status   - A status code: 0=success.
//	Dim1     - The first dimension size.  Optional.
//	Dim2     - The second dimension size.  Optional.
//	IndFmt   - A bit-coded optional flag indicating how the array should
//	           be indexed in each axis; 0=0-based, 1=1-based.  defaults
//	           to 0---all 0-based.  This is only relevant if this routine
//	           allocates the array.
//
// Written by Michael R. Greason, SSAI, 07 February 2006.
// ============================================================================

#include "port_helper.h"

#include "fitsio.h"

#include <iostream>
#include <string>

using namespace std;

void  Read_FITS_REAL_2D ( string FileName,  real4_2d& Arr, int& Status, int *Dim1, int *Dim2, int IndFmt)
{
  if ( IndFmt != 0 )  //X c++ port: this is the only case that has been ported
  {
    Status = -1;
    return;
  }

  int I, D1, D2, NEle, Ifmt;
  int AnyF;
  fitsfile *fptr;
  // ----------------------------------------------------------------------------
  Ifmt = 0;
  Ifmt = IndFmt;
  //
  //			Open the FITS file.
  //
  Status = 0;
  I      = 0;

  fits_open_file(  &fptr, FileName.c_str(), I, &Status );
  //
  //			Allocate space for the array in the primary HDU.
  //
  fits_get_img_dim( fptr, &D1, &Status );

  if ((Status  !=  0)  ||  (D1  <  2))
  {
    if (D1  <  2)
      cout << "Read_FITS_REAL_2D:  Status = " << Status << ", Primary HDU has insufficient axes." << endl;
    return;
  }

  float atmp, btmp;
  fits_read_key_flt( fptr, "NAXIS1", &atmp, NULL, &Status );
  fits_read_key_flt( fptr, "NAXIS2", &btmp, NULL, &Status );
  D1 = ( int ) atmp;
  D2 = ( int ) btmp;
  if (Status  !=  0)
    return;
  if (Dim1!=NULL)
    *Dim1 = D1;
  if (Dim2!=NULL)
    *Dim2 = D2;

  //X no need to port... - gr
  // If (.NOT. Associated(Arr)) Then
  //   Select Case (Ifmt)
  //     Case (3)
  //       Allocate (Arr(1:D1,     1:D2    ), Stat=Status)
  //     Case (2)
  //       Allocate (Arr(0:(D1-1), 1:D2    ), Stat=Status)
  //     Case (1)
  //       Allocate (Arr(1:D1,     0:(D2-1)), Stat=Status)
  //     Case Default
  //       Allocate (Arr(0:(D1-1), 0:(D2-1)), Stat=Status)
  //   End Select
  //   If (Status .NE. 0) Return
  // End If
  // if (Status  !=  0)  return;
  Arr.resize( D1, D2 );

  //
  //			Fill the array.
  //
  NEle = D1 * D2;
  fits_read_img_flt( fptr, 0, 1, NEle, 0, Arr.data(), &AnyF, &Status );
  //
  //			Close the FITS file.
  //
  fits_close_file( fptr, &Status);
  //
  return;
  // ----------------------------------------------------------------------------
  //
} //  Read_FITS_REAL_2D

// 'Read_FITS_Real_2D' reads a two-dimensional floating point array from the
// primary header/data unit in a FITS file.  This routine allocates the array
// from the heap if it hasn't already been associated.
//
// If the array has already been created, this routine assumes that it has
// the correct size and dimensions!  If allocated here then it will be indexed
// like a C array:  0:n-1.
//
// Arguments:
//	FileName - The name of the FITS file.
//	Arr      - The array to store the data in.
//	Status   - A status code: 0=success.
//	Dim1     - The first dimension size.  Optional.
//	Dim2     - The second dimension size.  Optional.
//	IndFmt   - A bit-coded optional flag indicating how the array should
//	           be indexed in each axis; 0=0-based, 1=1-based.  Defaults
//	           to 0---all 0-based.  This is only relevant if this routine
//	           allocates the array.
//
// Written by Michael R. Greason, SSAI, 07 February 2006.
// ============================================================================
void  Read_FITS_Double_2D ( string FileName,  double_2d& Arr, int& Status, int *Dim1, int *Dim2, int IndFmt)
{
  if ( IndFmt != 0 )  //X c++ port: this is the only case that has been ported
  {
    Status = -1;
    return;
  }

  int I, D1, D2, NEle, Ifmt;
  int AnyF;
  fitsfile *fptr;
  // ----------------------------------------------------------------------------
  Ifmt = 0;
  Ifmt = IndFmt;
  //
  //			Open the FITS file.
  //
  Status = 0;
  I      = 0;

  fits_open_file(  &fptr, FileName.c_str(), I, &Status );
  //
  //			Allocate space for the array in the primary HDU.
  //
  fits_get_img_dim( fptr, &D1, &Status );

  if ((Status  !=  0)  ||  (D1  <  2))
  {
    if (D1  <  2)
      cout << "Read_FITS_REAL_2D:  Status = " << Status << ", Primary HDU has insufficient axes." << endl;
    return;
  }

  float atmp, btmp;
  fits_read_key_flt( fptr, "NAXIS1", &atmp, NULL, &Status );
  fits_read_key_flt( fptr, "NAXIS2", &btmp, NULL, &Status );
  D1 = ( int ) atmp;
  D2 = ( int ) btmp;
  if (Status  !=  0)
    return;
  if (Dim1!=NULL)
    *Dim1 = D1;
  if (Dim2!=NULL)
    *Dim2 = D2;

  //X no need to port... - gr
  // If (.NOT. Associated(Arr)) Then
  //   Select Case (Ifmt)
  //     Case (3)
  //       Allocate (Arr(1:D1,     1:D2    ), Stat=Status)
  //     Case (2)
  //       Allocate (Arr(0:(D1-1), 1:D2    ), Stat=Status)
  //     Case (1)
  //       Allocate (Arr(1:D1,     0:(D2-1)), Stat=Status)
  //     Case Default
  //       Allocate (Arr(0:(D1-1), 0:(D2-1)), Stat=Status)
  //   End Select
  //   If (Status .NE. 0) Return
  // End If
  // if (Status  !=  0)  return;
  Arr.resize( D1, D2 );

  //
  //			Fill the array.
  //
  NEle = D1 * D2;
  fits_read_img_dbl( fptr, 0, 1, NEle, 0, Arr.data(), &AnyF, &Status );
  //
  //			Close the FITS file.
  //
  fits_close_file( fptr, &Status);
  //
  return;
  // ----------------------------------------------------------------------------
  //
} //  Read_FITS_REAL_2D

// ============================================================================

// ============================================================================
// 'Read_FITS_Real_3D' reads a three-dimensional floating point array from the
// primary header/data unit in a FITS file.  This routine allocates the array
// from the heap if it hasn't already been associated.
//
// If the array has already been created, this routine assumes that it has
// the correct size and dimensions!  If allocated here then it will be indexed
// like a C array:  0:n-1.
//
// Arguments:
//	FileName - The name of the FITS file.
//	Arr      - The array to store the data in.
//	Status   - A status code: 0=success.
//	Dim1     - The first dimension size.  Optional.
//	Dim2     - The second dimension size.  Optional.
//	Dim3     - The third dimension size.  Optional.
//	IndFmt   - A bit-coded optional flag indicating how the array should
//	           be indexed in each axis; 0=0-based, 1=1-based.  defaults
//	           to 0---all 0-based.  This is only relevant if this routine
//	           allocates the array.
//
// Written by Michael R. Greason, SSAI, 07 February 2006.
// ============================================================================
void  Read_FITS_REAL_3D (string FileName, real4_3d& Arr, int& Status, int* Dim1, int* Dim2, int* Dim3, int IndFmt)
{
  //X index fun: we specialize in lndFmt = 0 for the c++ port and adapt the calling routine accordingly    -- gr
  if ( IndFmt != 0 )
  {
    Status = -1;
    return;
  }

  int  I, D1, D2, D3, NEle, Ifmt;
  int  AnyF;
  fitsfile* fptr;
  // ----------------------------------------------------------------------------
  Ifmt = IndFmt;
  //
  //			Open the FITS file.
  //
  Status = 0;
  I      = 0;
  fits_open_file(  &fptr, FileName.c_str(), I, &Status );
  //
  //			Allocate space for the array in the primary HDU.
  //
  fits_get_img_dim( fptr, &D1, &Status );
  if (((Status)  !=  0)  ||  (D1  <  3))
  {
    if (D1  <  3)
      cout << "Read_FITS_REAL_3D:  Status = " << Status << ", Primary HDU has insufficient axes.\n";
    return;
  }

  float atmp;
  fits_read_key_flt( fptr, "NAXIS1", &atmp, NULL, &Status );
  D1 = ( int ) atmp;
  fits_read_key_flt( fptr, "NAXIS2", &atmp, NULL, &Status );
  D2 = ( int ) atmp;
  fits_read_key_flt( fptr, "NAXIS3", &atmp, NULL, &Status );
  D3 = ( int ) atmp;
  if (Status  !=  0)
    return;
  if (Dim1 != NULL)
    *Dim1 = D1;
  if (Dim2 != NULL)
    *Dim2 = D2;
  if (Dim3 != NULL)
    *Dim3 = D3;
  //X not ported, see above - gr
  //X if ( ! Associated(Arr))  {
  //X  Select case (Ifmt);
  //X    case (7);
  //X      Allocate (Arr(1:D1,     1:D2,     1:D3    ), Stat=Status);
  //X    case (6);
  //X      Allocate (Arr(0:(D1-1), 1:D2,     1:D3    ), Stat=Status);
  //X    case (5);
  //X      Allocate (Arr(1:D1,     0:(D2-1), 1:D3    ), Stat=Status);
  //X    case (4);
  //X      Allocate (Arr(0:(D1-1), 0:(D2-1), 1:D3    ), Stat=Status);
  //X     case (3);
  //X       Allocate (Arr(1:D1,     1:D2,     0:(D3-1)), Stat=Status);
  //X     case (2);
  //X       Allocate (Arr(0:(D1-1), 1:D2,     0:(D3-1)), Stat=Status);
  //X     case (1);
  //X       Allocate (Arr(1:D1,     0:(D2-1), 0:(D3-1)), Stat=Status);
  //X     case default;
  //X       Allocate (Arr(0:(D1-1), 0:(D2-1), 0:(D3-1)), Stat=Status);
  //X   End Select;
  //X   if (Status  !=  0)  return
  //X }
  //
  //			Fill the array.
  //
  NEle = D1 * D2 * D3;
  Arr.resize( D1, D2, D3 );
  fits_read_img_flt( fptr, 0, 1, NEle, 0, Arr.data(), &AnyF, &Status );
  //
  //			Close the FITS file.
  //
  fits_close_file( fptr, &Status);
  //
  return;
  // ----------------------------------------------------------------------------
  //
} //  Read_FITS_REAL_3D


//X // ============================================================================
//X // 'Read_FITS_Complex_2D' reads a two-dimensional complex array from the primary
//X // header/data unit in a FITS file.  This routine allocates the array from the
//X // heap if it hasn't already been associated.
//X //
//X // If the array has already been created, this routine assumes that it has
//X // the correct size and dimensions!  If allocated here then it will be indexed
//X // like a C array:  0:n-1.
//X //
//X // FITS does not directly support complex data; it treats each element as a
//X // two-element array.  Therefore, if the complex array has dimensions MxN it is
//X // stored as a 2xMxN array in the FITS file.  This routine will use
//X // 'Read_FITS_Real_3D' to read the file into a 3D array; this array will be
//X // used to fill the output array.
//X //
//X // Arguments:
//X //	FileName - The name of the FITS file.
//X //	Arr      - The array to store the data in.
//X //	Status   - A status code: 0=success.
//X //	Dim1     - The first dimension size.  Optional.
//X //	Dim2     - The second dimension size.  Optional.
//X //	IndFmt   - A bit-coded optional flag indicating how the array should
//X //	           be indexed in each axis; 0=0-based, 1=1-based.  defaults
//X //	           to 0---all 0-based.
//X //
//X // Written by Michael R. Greason, SSAI, 07 February 2006.
//X // ============================================================================
//X void  Read_FITS_Complex_2D (FileName, Arr, Status, Dim1, Dim2, IndFmt)
//X {
//X
//X //
//X string(*),                 Intent(In)   FileName;
//X Complex (Kind=4), Pointer, Dimension(:,:)  Arr;
//X int (Kind=4),             Intent(Out)  Status;
//X int (Kind=4),   Optional, Intent(Out)  Dim1;
//X int (Kind=4),   Optional, Intent(Out)  Dim2;
//X int (Kind=4),   Optional, Intent(In)   IndFmt;
//X //
//X int (Kind=4)  DC, D1, D2, I, I0, I1, J, J0, J1, Ifmt;
//X REAL_4, Pointer, Dimension(:,:,:)  Tmp;
//X // ----------------------------------------------------------------------------
//X Ifmt = 0;
//X if (Present(IndFmt))  Ifmt = IndFmt;
//X //
//X //			Read the file into a temporary array.
//X //
//X I = Ifmt * 2 + 1;
//X Call Read_FITS_REAL_3D (FileName, Tmp, Status, DC, D1, D2, I);
//X if (*Status  !=  0)  return
//X if ((DC  !=  2)  ||  (D1  <=  0)  ||  (D2  <=  0))  {
//X   (*Status) = 1;
//X   return;
//X }
//X if (Present(Dim1))  Dim1 = D1;
//X if (Present(Dim2))  Dim2 = D2;
//X //
//X //			Allocate and fill the output array.
//X //
//X Select case (Ifmt);
//X   case (3);
//X     I0 = 1;
//X     I1 = D1;
//X     J0 = 1;
//X     J1 = D2;
//X   case (2);
//X     I0 = 0;
//X     I1 = D1 - 1;
//X     J0 = 1;
//X     J1 = D2;
//X   case (1);
//X     I0 = 1;
//X     I1 = D1;
//X     J0 = 0;
//X     J1 = D2 -1;
//X   case default;
//X     I0 = 0;
//X     I1 = D1 - 1;
//X     J0 = 0;
//X     J1 = D2 -1;
//X End Select;
//X if ( ! Associated(Arr))  {
//X   Allocate (Arr(I0:I1,J0:J1), Stat=Status);
//X   if ((*Status)  !=  0)  return
//X }
//X //
//X for (J = J0; J <= J1; ++J ) {
//X   for (I = I0; I <= I1; ++I ) {
//X     Arr(I,J) = cmplx(Tmp(1,I,J), Tmp(2,I,J), Kind=4);
//X   }
//X }
//X //
//X Deallocate(Tmp);
//X //
//X return;
//X // ----------------------------------------------------------------------------
//X } //  Read_FITS_Complex_2D

// ============================================================================
// 'Read_FITS_Complex_2D_LM' reads a two-dimensional complex array from the
// primary header/data unit in a FITS file.  This routine reads and fills the
// array one element at a time in an attempt to conserve memory at the cost of
// performance.
//
// If the array has already been created, this routine assumes that it has
// the correct size and dimensions!  If allocated here then it will be indexed
// like a C array:  0:n-1.
//
// FITS does not directly support complex data; it treats each element as a
// two-element array.  Therefore, if the complex array has dimensions MxN it is
// stored as a 2xMxN array in the FITS file.
//
// Arguments:
//	FileName - The name of the FITS file.
//	Arr      - The array to store the data in.
//	Status   - A status code: 0=success.
//	Dim1     - The first dimension size.  Optional.
//	Dim2     - The second dimension size.  Optional.
//	IndFmt   - A bit-coded optional flag indicating how the array should
//	           be indexed in each axis; 0=0-based, 1=1-based.  defaults
//	           to 0---all 0-based.
//
// Written by Michael R. Greason, SSAI, 07 February 2006.
// ============================================================================
void  Read_FITS_Complex_2D_LM ( string FileName, complex_2d& Arr, int& Status, int* Dim1, int* Dim2, int IndFmt)
{
  // index fun: we specialize in lndFmt = 0 and adapt the calling routine accordingly    -- gr
  if ( IndFmt != 0 )
  {
    Status = -1;
    return;
  }

  int  D1, D2, I, I0, I1, J, J0, J1, K;
  int  Ifmt;
  int  AnyF;
  REAL_4  Tmp[2];
  // ----------------------------------------------------------------------------
  Ifmt = 0;
  Ifmt = IndFmt;
  //
  //			Open the FITS file.
  //
  Status = 0;
  I      = 0;
  fitsfile* fptr;
  fits_open_file(  &fptr, FileName.c_str(), I, &Status );
  //
  //			Allocate space for the array in the primary HDU.
  //
  fits_get_img_dim( fptr, &D1, &Status );
  if ((Status  !=  0)  ||  (D1  <  3))
  {
    if (D1  <  3)
      cout << "Read_FITS_Complex_2D_LM:  Status = " << Status << ", Primary HDU has insufficient axes." << endl;
    return;
  }

  float atmp;
  fits_read_key_flt( fptr, "NAXIS2", &atmp, NULL, &Status );
  D1 = (int) atmp;
  fits_read_key_flt( fptr, "NAXIS3", &atmp, NULL, &Status );
  D2 = (int) atmp;
  if (Status  !=  0)
    return;
  if (NULL != Dim1)
    *Dim1 = D1;
  if (NULL != Dim2)
    *Dim2 = D2;

  //X  always default case, see above ... -- gr
  //X  Select case (Ifmt);
  //X    case (3);
  //X      I0 = 1;
  //X      I1 = D1;
  //X      J0 = 1;
  //X      J1 = D2;
  //X    case (2);
  //X      I0 = 0;
  //X      I1 = D1 - 1;
  //X      J0 = 1;
  //X      J1 = D2;
  //X    case (1);
  //X      I0 = 1;
  //X      I1 = D1;
  //X      J0 = 0;
  //X      J1 = D2 - 1;
  //X    case default;
  //X      I0 = 0;
  //X      I1 = D1 - 1;
  //X      J0 = 0;
  //X      J1 = D2 - 1;
  //X  End Select;
  I0 = 0;
  I1 = D1 - 1;
  J0 = 0;
  J1 = D2 - 1;


  if (  Arr.size( ) < ( unsigned int ) D1*D2 )
  {
    cout << "arr.size() < D1*D2 (" << Arr.size() << " < " << (D1*D2) << ") in read_fits.cpp " << endl;
    Status = -1;
    return;
  }


  //
  //			Fill the array.
  //
  K = -1;
  for (J = J0; J <= J1; ++J )
  {
    for (I = I0; I <= I1; ++I )
    {
      K = K + 2;
      fits_read_img_flt( fptr, 0, K, 2, 0, &Tmp[0], &AnyF, &Status );
      if (Status  !=  0)
        goto exitCmplxLoop;
      Arr(I,J) = complex<REAL_4>(Tmp[0], Tmp[1]);
    }
  }
exitCmplxLoop:
  //
  //			Close the FITS file.
  //
  fits_close_file( fptr, &Status);
  //
  return;
  // ----------------------------------------------------------------------------
} //  Read_FITS_Complex_2D_LM

