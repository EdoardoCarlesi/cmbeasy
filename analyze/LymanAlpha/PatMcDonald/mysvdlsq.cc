 
#include "myutils.h"
#include "mynrutils.h"

namespace NR {

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#define NRANSI

// y = [p] a + error, except p seems backwards, i.e., really transpose[p]
void myweightedsvdlsq(double **p, double y[],double sig[], int ndata, 
              double a[], int ma, 
              double *chisq, double *maxoff, int *nuse,
              double TOL,bool verbose)
{
	int j,i;
	double wmax,tmp,sum,*b,thresh;
	double **u, **v, *w;
        u=NR::dmatrix(1,ndata,1,ma);
        v=NR::dmatrix(1,ma,1,ma);
        w=NR::dvector(1,ma);
        *maxoff=0.0;
        *nuse=0;
	b=NR::dvector(1,ndata);
	for (i=1;i<=ndata;i++) {
		tmp=1.0/sig[i];
		for (j=1;j<=ma;j++) u[i][j]=p[j][i]*tmp;
		b[i]=y[i]*tmp;
	}
	NR::svdcmp(u,ndata,ma,w,v,DMIN(0.1*TOL,mytinyfloat));
                           //numeric_limits<float>::epsilon()));
	wmax=0.0;
	for (j=1;j<=ma;j++)
		if (w[j] > wmax) wmax=w[j];
	thresh=TOL*wmax;
	for (j=1;j<=ma;j++){
          if(verbose && w[j]>0.0 ) printf("%g %g\n",w[j],thresh);
	  if (w[j] < thresh) { w[j]=0.0; }
          else (*nuse)++;
        }
	NR::svbksb(u,w,v,ndata,ma,b,a);
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*p[j][i];
		*chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
                 if(fabs(y[i]-sum)>fabs(*maxoff)) *maxoff=y[i]-sum; 
                 int ti = (i-1)/125;
                 if(verbose) printf("%d %d %d %g %g %g\n",4*(ti/3)+
                                      int(fmod(ti,3.0))+1,
                                     int(fmod(i-1.0,125.0))/25,
                               int(fmod(i-1.0,25.0)),sum-y[i],y[i],sum);
	}
	NR::free_dvector(b,1,ndata);
        NR::free_dmatrix(u,1,ndata,1,ma);
        NR::free_dmatrix(v,1,ma,1,ma);
        NR::free_dvector(w,1,ma);
}
#undef TOL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software !(^$S_!~. */
void mysvdlsq(double **p, double y[], int ndata, 
              double a[], int ma, 
              double *chisq, double *maxoff, int *nuse,
              double TOL,bool verbose)
{
  double* sig;
  sig=NR::dvector(1,ndata);
  for(int i=1;i<=ndata;i++) sig[i]=1;
  myweightedsvdlsq(p,y,sig,ndata,a,ma,chisq,maxoff,nuse,TOL,verbose);
  NR::free_dvector(sig,1,ndata);
  return;
}

}
