
#ifndef NR_FLAG
#define NR_FLAG

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "CPP.h"
namespace NR {
 
  void myweightedsvdlsq(double **p, double y[],double sig[], int ndata, 
                double a[], int ma,
                double *chisq,double *maxoff, int *nuse, double TOL,
                bool verbose); 
  void mysvdlsq(double **p, double y[], int ndata, double a[], int ma,
                double *chisq,double *maxoff, int *nuse, double TOL,
                bool verbose); 
  double pythag(double a, double b); 
  void svbksb(double **u, double w[], double **v, int m, int n, double b[],
              double x[]);
  void svdcmp(double **a, int m, int n, double w[], double **v,double TOL);
  double chebev(double a, double b, double c[], int m, double x); 
  void linterp(double *xa,double *ya,int n,double x,unsigned long *jlo,
               double *y);
  void hunt(double xx[], unsigned long n, double x, unsigned long *jlo); 

  static double dsqrarg;
  #define DSQR(a) (DEQUAL((dsqrarg=(a)),0.0) ? 0.0 : dsqrarg*dsqrarg)

  #define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
  static double Maxarg1,Maxarg2;
  #define DMAX(a,b) (NR::Maxarg1=(a),NR::Maxarg2=(b),(NR::Maxarg1) > \
                (NR::Maxarg2) ?\
		(NR::Maxarg1) : (NR::Maxarg2))
  static double Minarg1,Minarg2;
  #define DMIN(a,b) (NR::Minarg1=(a),NR::Minarg2=(b),(NR::Minarg1) > \
                (NR::Minarg2) ?\
		(NR::Minarg2) : (NR::Minarg1))
  static int IMinarg1,IMinarg2;
  #define IMIN(a,b) (IMinarg1=(a),IMinarg2=(b),(IMinarg1) > (IMinarg2) ?\
		(IMinarg2) : (IMinarg1))

  void gaussj(double **a, int n, double **b, int m);
  void lubksb(double **a, int n, int *indx, double b[]);
  int *ivector(long nl, long nh);
  void free_ivector(int *v, long nl, long nh);
  int **imatrix(long nrl, long nrh, long ncl, long nch);
  void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
  void ludcmp(double **a, int n, int *indx, double *d);
  double *dvector(long nl, long nh);
  void free_dvector(double *v, long nl, long nh);
  double **dmatrix(long nrl, long nrh, long ncl, long nch);
  void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
  void nrerror(char error_text[]);
  void spline(const vector<double>& x, const vector<double>& y, int n, 
              double yp1, double ypn,vector<double>& y2);
  void splint(const vector<double>& xa, const vector<double>& ya, const vector<double>& y2a, int n, double x, 
              double *y);
  void oldspline(double* x, double* y, int n,double yp1, double ypn,double* y2);
  void oldsplint(double* xa,double* ya,double* y2a, int n, double x, double *y);
  void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a);
  void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n,
        double x1, double x2, double *y);
  void polint(double xa[], double ya[], int n, double x, double *y,
                    double *dy); 
  double ***df3tensor(long nrl, long nrh, long ncl, long nch, long ndl, 
                   long ndh);
  void free_df3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh);
  double plgndr(int l, int m, double x);
  void covsrt(double **covar, int ma, int ia[], int mfit);
}

#endif


