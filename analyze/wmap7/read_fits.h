#ifndef READ_FITS_H
#define READ_FITS_H
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
// Interface block for the routines in read_fits.f90
// ============================================================================
//

void  Read_FITS_REAL_2D ( std::string FileName, real4_2d& Arr, int& Status, int* Dim1 = NULL, int* Dim2 = NULL, int IndFmt = 0);
void  Read_FITS_Double_2D ( std::string FileName, double_2d& Arr, int& Status, int* Dim1 = NULL, int* Dim2 = NULL, int IndFmt = 0);
void  Read_FITS_REAL_3D ( std::string FileName, real4_3d& Arr, int& Status, int* Dim1 = NULL, int* Dim2 = NULL, int *Dim3 = NULL, int IndFmt = 0);
void  Read_FITS_Complex_2D ( std::string FileName, complex_2d& Arr, int& Status, int* Dim1 = NULL, int* Dim2 = NULL, int IndFmt = 0);
void  Read_FITS_Complex_2D_LM ( std::string FileName, complex_2d& Arr, int& Status, int* Dim1 = NULL, int* Dim2 = NULL, int IndFmt = 0);

#endif // READ_FITS_H
