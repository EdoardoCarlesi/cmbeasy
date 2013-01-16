#ifndef MISCMATH
#define MISCMATH

#include "mathobject.h"

#include <iostream>
#include <map>
#include <list>
#include <cstdlib>

#define MISCMATH_BS_MAXNVAR 1000
#define MISCMATH_BS_KMAXX 8


/*!
  Miscellaneous mathematics class. 

  It mostly consists of C++ adopted Numerical Recipees routines.
  
  If you do not like those, feel encouraged to provide alternative
  implementations or algorithms. I would be happy to include those
  in the code.
 */
class Miscmath : public Mathobject {
  static double splg[20001];        //! splini and splder need it...
  static bool splgValid;
  static bool seedless; //!< random seed
  static double rand_max; // 1.0/RAND_MAX

  
  //static const int burlish;
  //static const int runge;

  static double bs_d[MISCMATH_BS_MAXNVAR+1][MISCMATH_BS_KMAXX+1];

  static double bs_x[MISCMATH_BS_KMAXX+1];
  static double bs_1d[MISCMATH_BS_KMAXX+1];

   double h,a,b;
   double DeltaChi2_Nu;
   double Gamma_z;
 
 public:
   enum keyword {rkutta, bstoer};

  Miscmath();
  static double rombint(moSingle , Mathobject * ,double,double,double); //!< Romberg integration


  // ODE
  //! Ordinary differential equation stepper.
  static double odeint(double*, const int , const double , const double , const double , const  double , const double , moDerivs , const Mathobject&,bool=true,keyword which = rkutta);

  //! Calls odeint() with runge-kutta switched on
  static double rungeInt(double* y, const int n, const double x1, const double x2 , const double eps, const  double h1 , const double hmin, moDerivs d, const Mathobject& mo) { 
    return odeint(y,n,x1,x2,eps,h1,hmin,d,mo);
  }
  //! Calss odeint() with burlish-stoer switched on
  static double burlishInt(double* y, const int n, const double x1, const double x2 , const double eps, const  double h1 , const double hmin, moDerivs d, const Mathobject& mo) {
    return odeint(y,n,x1,x2,eps,h1,hmin,d,mo,true,bstoer); 
  }

double solid_angle(double the1, double the2, double phi1, double phi2);

  // RUNGE KUTTA 
  static void rkck(double *, double *, const int, const double, const double, double*,double*, moDerivs, const Mathobject&,const bool) ; //!< for runge kutta 
  static void rkqs(double*, double*, const int, double*, const double, const double , double*, double *, double *, moDerivs , const Mathobject&,const bool ) ; //!< runge kutta stepper

  // BURLISH STOER

  static void bsstep(double *, double *, const int, double *, const double, const double, double *, double *, double *, moDerivs, const Mathobject  &); //!< burlish-stoer stepper
  static void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep, double yout[], moDerivs , const Mathobject& );
  static void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
  static void rzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);


  // One dimensional version for convolution and
  // integration of f(x) of the odeint and bsstep etc

  static double oneDimOdeint(const double , const double , const double , moSingle , const Mathobject&) ;
  static void oneDimBsstep(double *, const double, double *, const double, const double, const double , double *, double *, moSingle, const Mathobject &);
  static double oneDimMmid(const double , const double dydx, const double xs, const double htot, const int nstep, moSingle , const Mathobject& );
  static void oneDimPzextr(const int , const double, const double, double *, double *);
  

  // FINDING ZEROS OF A FUNCTION

  static double zbrent(moSingle func, const Mathobject &mo,double x1, double x2, double tol=1e-10, double ZeroLevel=0.0); //!< find zero of a function actually, you may shift the zero-level 
  
  
  //  DERIVATIVE OF A FUNCTION
  static double dfridr(moSingle func,Mathobject& mo,double x, double h, double *err);  //!< derive a function



 static void splini(); //! initialize whatever...
 static void splder(double*, double*, const int,bool = true); //! Fits cubic splines to y and interpolates first derivs at grid points dy
 static void splint(double*, double*, const int); //! Spline integration from CMBFAST, please do not confuse with splint() from below

 static void spline(double*, double*, const int, const double, const double, double*); //! Generates interpolating coefficients for subsequent interpolation with splint Nrecipes
 static void spline(double* x, double* y, const int n, double *yp1, double* ypn, double* y2) { spline(x,y,n,*yp1,*ypn,y2); } //!< wrapper
 // void origSpline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
 static double splint(double xa[], double ya[], double y2a[], int n, double x, int *guess = 0) ; //! Spline interpolation Nrecipes + guess for next

 static double rnd(double=1.0); //!< return random number within [-bound, bound]
 static double posRnd(double=1.0); //!< return positive random number within [0, bound]
 static void seed(); //!< seed the random seed
 static double gaussRnd(const double mean = 0.0, const double std =1.0); //!< return gaussian random variable with mean and std standard deviation


 static double sign(double x) { if (x < 0) return -1.0; else return 1.0;} //!< return -1 if x < 0, else return 1.0;

 double fractionCorrespondingToSigma(double sigma); //!< return the fraction under the gaussian distribution function corresponding to so and so many sigma
 double fctsint(const double x) const; //!< helper for fractionCorrespondingToSigma

 static double relError(const double x, const double y);
 static double sigma(const std::map<double,int>&,double s=0.67); //!< for a histogram running from 0... some value, return the value at which the total number exceeds s*m.size();

 static std::pair<double,double> Sigma(const std::map<double,int>&,double s=0.67); //!< return left and right value of sigma percent of the histogram.


 static std::map<double,int> histogram(std::list<double>& l,double a,double b,int k=0);
 static std::map<double,int> histogram(std::list<double>&l, std::list<double>& bins);

 static void printHistogram(const std::map<double, int>&, std::ostream& o=std::cout);

 static double topHatGaussian(double x,double left,double right,double lgauss, double rgauss);

 double DeltaChi2CorrespondingToSigma(double nu, double sigma);
 double DeltaChi2Helper(double chi2);
 double GammaQ(double a, double x);
 double GammaInt(double t);

 static double removePole(double denom, double epsilon); // cut off pole in 1/denom

 //! performs logarithmic addition, i.e. returns log(exp(x)+exp(y))
 static double logAdd(double x, double y);
};



inline double  Miscmath::posRnd(double max) {
  if (seedless)  seed();
  return  rand() * rand_max *max;
}


#endif
