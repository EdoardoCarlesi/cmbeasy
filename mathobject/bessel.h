#ifndef bESSEL_H
#define bESSEL_H 

#include "global.h"

/*! C++ Wrapper of the Numerical Recipes 
  Bessel routines and handler for bessel functions. */

class Bessel {
 private:
 
  public:  
  Bessel();
  
  
 
  static double  bessj0(double x); //!< J_0(x)
  static double  bessj1(double x); //!< J_1(x)
  static double  bessj(int n, double x);  // J_n(x) with n >= 2
  static double  bessy0(double x); //<! Y_0(x), 

  static double  bessy1(double x); //<! Y_1(x)
  static double  bessy(int n, double x);
  static void bessjy(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp);
  static void sphbes(int n, double x, double *,double*,double*,double*);
  static void beschb(double x, double *gam1, double *gam2, double *gampl,double *gammi);
  static double chebev(double a, double b, double c[], int m, double x);
  static double fractBessj(double nu, double x);
  static double sphJ(int n, double x); //!< convenience wrapper for sphbes, return j_n(x) 
  static float bjl(int l, float x); // CMBFAST' own version of bessj() 
};
#endif

