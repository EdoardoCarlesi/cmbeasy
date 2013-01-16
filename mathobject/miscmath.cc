#include "miscmath.h"
#include "spline.h"
#include "global.h"

#include <cstdlib>
#include <cmath>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
//#define PI 3.141592653589
#define PI 4*atan(1.)
 
Miscmath::Miscmath()  {}

/*
 * solid angle returns integral value of dOmega over the four given radial coordinates
 * */

double Miscmath::solid_angle(double the1, double the2, double phi1, double phi2){
the1 *= PI/180; the2 *= PI/180;
phi1 *= PI/180; phi2 *= PI/180;
double x, step; int pts=100;
Spline *sine = new Spline(pts, "sine");
//cout << "pi: " << PI << " t1, t2: " << the1 << " " << the2 << endl;

step = (the2-the1)/(pts-1);
x=the1;

for(int i=0; i<pts; i++){
double y = sin(x);
sine->set(x,y);
//cout << " x: " << x*(180/PI) << " y: " << y << endl;
x += step;
}
sine->arm();
double integral = sine->integrate(the1,the2);
delete sine;
//cout << "Integral solid angle: " << integral*(phi2-phi1) << endl;
return integral*(phi2-phi1);
}


/*!  Rombint returns the integral from a to b of using Romberg integration.
     The method converges provided that f(x) is continuous in (a,b). 
     f must be double precision and must be declared external in the calling
     routine.  tol indicates the desired relative accuracy in the integral. */


double Miscmath::rombint(moSingle f, Mathobject *mo,const double a, const double b, const double tol) {
    // System generated locals 

   double gmax;
   int jmax, nint;
   double g[6], h;
   int i, j, k;
   double fourj, g0, g1, error;

    h = (b - a) * .5;
    gmax = h * ((mo->*f)(a) + (mo->*f)(b));

    g[0] = gmax;
    nint = 1;
    error = 1e20;
    i = 0;
    while(true) {
      ++i;
      if ((i > 40) || (i > 5 && fabs(error) < tol)) break;
      
      //  Calculate next trapezoidal rule approximation to integral. 
      g0 = 0.;
      for (k = 1; k <= nint; ++k)  g0 += (mo->*f)( a + (k + k - 1) * h );
	  

      g0 = g[0] * .5 + h * g0;
      h *= .5;
      nint += nint;
      jmax = min(i,5);
      fourj = 1.;
      
      for (j = 1; j <= jmax; ++j) {
	//  Use Richardson extrapolation. 
	fourj *= 4.;
	g1 = g0 + (g0 - g[j - 1]) / (fourj - 1.);
	g[j - 1] = g0;
	g0 = g1;
      }
      if (fabs(g0) > tol) error = 1. - gmax / g0; else error = gmax;
      gmax = g0;
      g[jmax] = g0;      
    }
   
    if (i > 40 && fabs(error) > tol) {
      cout << "Rombint failed to converge; integral" << g0 << "error = " << error << endl;
      throw Bad_Error("Rombint failed to converge");  
    }
    return g0;
} 



#define NVAR 4000


#define NRANSI
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define TINY 1.0e-30

/*! Modified stepper of NRecipes.
  
  The modification is that only errors that come from variables with
  ysca[] > TINY will contribute. In other words: If derivs is 
  bad behaved and y[] and dxdy[] are both zero on either side 
  of a threshold of x (say as a function of time), then the error 
  would explode by just naively extrapolating the steps. */

void Miscmath::rkqs(double *y, double *dydx, const int n, double *x, const double htry, const double eps,
	double *yscal, double *hdid, double *hnext, moDerivs derivs, const Mathobject& mo, const bool original)  {
	int i;
	double errmax,h,htemp,xnew;
//cout << "Miscmath::rkqs init" << endl;
//	  for(int k=0; k<10; k++) cout << "y[" << k << "]: " << y[k] << endl; 
	double yerr[NVAR], ytemp[NVAR];

//	y[0]=y[5]=y[7]=0;

	h=htry;
	for (;;) {
      rkck(y,dydx,n,*x,h,ytemp,yerr,derivs,mo,original);
      errmax=0.0;
      double maxErrInd = -1;
		if (original) {
		  for (i=1;i<=n;i++) {
            errmax = max(errmax, fabs(yerr[i]/yscal[i]));
          }
		} else {
		  for (i=1;i<=n;i++) {
		    if (  fabs(yerr[i]/yscal[i]) > errmax ) {
		      if (original || (yscal[i] > TINY)) {
                errmax=max(errmax,fabs(yerr[i]/yscal[i])); 
                maxErrInd = i;
              }
		    }
		  }
		}
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? max(htemp,0.1*h) : min(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x) {
          cout << "Miscmath::rkqs() : maxErrInd=" << maxErrInd << endl;
	  for(int k=0; k<10; k++) cout << "y[" << k << "]: " << y[k] << endl; 
		cout << "" << endl;
          throw Bad_Error("stepsize underflow in rkqs; if running crossoverfield quintessence, J may be too low");
        }
	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=1;i<=n;i++) y[i]=ytemp[i];

}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI
/* note #undef's at end of file */
#define NRANSI


void Miscmath::rkck(double *y, double *dydx, const int n, const double x, const double h, double *yout,
	double *yerr, moDerivs derivs,const Mathobject& mo, const bool original)  {
	int i;
	double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
		b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
		dc5 = -277.00/14336.0;
	double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0,dc6=c6-0.25;


	double ak2[NVAR], ak3[NVAR], ak4[NVAR], ak5[NVAR], ak6[NVAR];
	double ytemp[NVAR];
	
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	(mo.*derivs)(x+a2*h,ytemp,ak2);
	// cout << "nummer 10: " << ak2[10] << endl;
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	(mo.*derivs)(x+a3*h,ytemp,ak3);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	(mo.*derivs)(x+a4*h,ytemp,ak4);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	(mo.*derivs)(x+a5*h,ytemp,ak5);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	(mo.*derivs)(x+a6*h,ytemp,ak6);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=1;i<=n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

}

#undef NRANSI

// #include <math.h>
#define NRANSI
// #include "nrutil.h"
//#define MAXSTP 1500000
#define MAXSTP 150000


double Miscmath::odeint(double *ystart, const int nvar, const double x1, const double x2, const double eps, const double h1,  const double hmin, moDerivs derivs, const Mathobject& mo, bool original, keyword which)  {
	 int nstp,i;
	 double x,hdid,h;
	 double hnext = h1;  // to have a proper initialization for bstoer
//cout << "Miscmath::odeint()" << endl;
	 double yscal[NVAR],y[NVAR], dydx[NVAR];
	//cout << "nvar: " << nvar << " x1: " << x1 << " x2: " << x2 << " eps: " << eps << " h1: " << h1 << " hmin: " << hmin << endl;
	if (nvar >= NVAR) throw Bad_Error("odeint()  too many variables, increase NVAR");
	x=x1;
	h=SIGN(h1,x2-x1);
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	for (nstp=1;nstp<=MAXSTP;nstp++) {
	  (mo.*derivs)(x,y,dydx);
//cout << "nvar: " << nvar << " y[" << i << "]: " << y[i] << " x: " << x << endl; 
	  for (i=1;i<=nvar;i++) yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
	 
	  if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
	 
	  switch (which) {
	  case rkutta:
	    rkqs(y,dydx,nvar,&x,h,eps, yscal ,&hdid,&hnext,derivs, mo,original);
	    break;
	  case bstoer:
	    bsstep(y,dydx,nvar,&x,h,eps, yscal ,&hdid,&hnext,derivs, mo);
	    break;
	  default:
	    rkqs(y,dydx,nvar,&x,h,eps, yscal ,&hdid,&hnext,derivs, mo,original);
	  }
	  
	  if ((x-x2)*(x2-x1) >= 0.0) {
	    for (i=1;i<=nvar;i++) ystart[i]=y[i];
	    return hnext;
	  }
	  if (fabs(hnext) <= hmin)  throw Bad_Error("Step size too small in odeint");
	  h=hnext;
	}
	throw Bad_Error("Too many steps in routine odeint");
}
#undef MAXSTP
//#undef TINY
#undef NRANSI

#undef NVAR               //! Undefining NVar




/*! 
  This ist the true (well, a bit modified) NR spline interpolation, sorry for the misnomer 
  of CMBFASTS splint which is an integration routine

  If "x" is out of the range xa[1] xa[n], then the boundary values will be given.
  Meaning:: You will not get any errormessages, if you are out of the boundary,
  just the closest available value 

*/
double Miscmath::splint(double xa[], double ya[], double y2a[], int n, double x, int *guess)  {
     
  int klo,khi,k;
  double h,b,a;
  // cout << "splint" << endl;
  if (x >= xa[n]) { 
    if (guess !=0) *guess = n-1; 
    return xa[n]; 
  }
  if (x <= xa[1]) {
    if (guess !=0)  *guess = 1; 
    return xa[1]; 
  }


  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) throw Bad_Error("Miscmath::splint() Bad xa input to routine splint");
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  if (guess != 0) *guess = klo;
  return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}



/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void Miscmath::splder(double *y, double *dy, const int n, bool adjust) {
    /* System generated locals */
    int  i__1;

    /* Local variables */
     double f[20001];
     int i, n1;

    if (!splgValid) splini();


/*  Splder fits a cubic spline to y and returns the first derivatives at */
/*  the grid points in dy.  Dy is equivalent to a 4th-order Pade */
/*  difference formula for dy/di. */



    /* Parameter adjustments if requested (compatibility mode :-) */
    if (adjust) {
      --dy;
      --y;
    }

    /* Function Body */
    n1 = n - 1;
    if (n1 > 20000) {
      cout << "Spline array overflow!!! n1= " << n1 << "  >20000" <<  endl;
      throw Bad_Error("subroutines::splder() stop");
    }
/*  Quartic fit to dy/di at boundaries, assuming d3y/di3=0. */
    f[0] = (y[1] * -10. + y[2] * 15. - y[3] * 6. + y[4]) / 6.;
    f[n - 1] = (y[n] * 10. - y[n1] * 15. + y[n - 2] * 6. - y[n - 3]) / 6.;
/*  Solve the tridiagonal system */
/*  dy(i-1)+4*dy(i)+dy(i+1)=3*(y(i+1)-y(i-1)), i=2,3,...,n1, */
/*  with dy(1)=f(1), dy(n)=f(n). */
    i__1 = n1;
    for (i = 2; i <= i__1; ++i) {
	f[i - 1] = splg[i - 1] * ((y[i + 1] - y[i - 1]) * 3. - f[i - 2]);
/* L10: */
    }
    dy[n] = f[n - 1];
    for (i = n1; i >= 1; --i) {
	dy[i] = f[i - 1] - splg[i - 1] * dy[i + 1];
/* L20: */
    }
} /* splder_ */


/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc*/
void Miscmath::splini() {
/*  Splini must be called before splder to initialize array g in common. 
*/


    splg[0] = 0.;
    for (int i = 2; i <= 20001; ++i) {
	splg[i - 1] = 1. / (4. - splg[i - 2]);
/* L10: */
    }
    splgValid = true;
} /* splini_ */



/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
void Miscmath::spline(double *x, double *y, const int n, 
	const double yp1, const double ypn, double *y2) {
    /* System generated locals */
    int i__1;

    /* Local variables */
    int i, k;
    double p, u[60000], qn, un, sig;  // 60000 in jlens.cc

    /* Parameter adjustments */
    --y2;
    --y;
    --x;

    /* Function Body */
    if (yp1 > 9.9e29) {
	y2[1] = 0.;
	u[0] = 0.;
    } else {
	y2[1] = -.5;
	u[0] = 3. / (x[2] - x[1]) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1);
    }
    i__1 = n - 1;
    for (i = 2; i <= i__1; ++i) {
	sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
	p = sig * y2[i - 1] + 2.;
	y2[i] = (sig - 1.) / p;
	u[i - 1] = (((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1])
		 / (x[i] - x[i - 1])) * 6. / (x[i + 1] - x[i - 1]) - sig * u[
		i - 2]) / p;
/* L11: */
    }
    if (ypn > 9.9e29) {
	qn = 0.;
	un = 0.;
    } else {
	qn = .5;
	un = 3. / (x[n] - x[n - 1]) * (ypn - (y[n] - y[n - 1]) / (x[n] 
		- x[n - 1]));
    }
    y2[n] = (un - qn * u[n - 2]) / (qn * y2[n - 1] + 1.);
    for (k = n - 1; k >= 1; --k) {
	y2[k] = y2[k] * y2[k + 1] + u[k - 1];
/* L12: */
    }
/*  (C) Copr. 1986-92 Numerical Recipes Software =$j*m,). */
} /* spline_ */


/* note #undef's at end of file */
#define NRANSI
//#include "nrutil.h"





/*!
  Splint integrates a cubic spline, providing the ouput value
 z = integral from 1 to n of s(i)di, where s(i) is the spline fit  to y(i).

  TO ME, IT SEEMS AS IF THIS ROUTINE REQUIRE EQUAL
  SPACING OF x[i], AS x[i] IS NEVER GIVEN, THIS IS 
  WHAT IT HAS TO ASSUME...

  AND IT SAYS SO:  s(i)di , not s(x)dx.... so one has to scale 
  and hope and pray that it is equally spaced. use CAREFULLY

*/ 
 void Miscmath::splint(double *y, double *z, const int n) {
   
    int i, n1;
    double dy1, dyn;

    --y;   // Parameter adjustments 

    n1 = n - 1;
    //  Cubic fit to dy/di at boundaries. */
    /* 	dy1=(-11.0d0*y(1)+18.0d0*y(2)-9.0d0*y(3)+2.0d0*y(4))/6.0d0 */
    dy1 = 0.;
    dyn = (y[n] * 11. - y[n1] * 18. + y[n - 2] * 9. - y[n - 3] * 2.) / 6.;

    *z = (y[1] + y[n]) * .5 + (dy1 - dyn) / 12.;
   
    for (i = 2; i <= n1; ++i) *z += y[i];

} /* splint_ */



double Miscmath::splg[20001];
bool Miscmath::splgValid = false;
bool Miscmath::seedless = true;
double Miscmath::rand_max = 1.0/RAND_MAX;


double Miscmath::rnd(double max) {
  if (seedless)  seed();
  return  2.0*(rand() * rand_max - 0.5)*max;
}

void Miscmath::seed() {
  srand(time(0));
  seedless = false;
}

double Miscmath::gaussRnd(const double mean,const double std) {
  
  double rsq,v1,v2;
  do {
     v1 = rnd(1.0);
     v2 = rnd(1.0);
     rsq = v1*v1 + v2*v2;
  } while (rsq >= 1.0 || rsq == 0.0);

  double fac=sqrt(-2.*log(rsq)/rsq);

  double  gasdev1 = v1*fac;
  //double  gasdev2 = v2*fac;
    
  double gauss1 = gasdev1*std + mean;
  // double gauss2 = gasdev2*std + mean;  idependent could return pair<double,double>

  return gauss1;
}



#include <math.h>
#define NRANSI
//#include "nrutil.h"
#define KMAXX MISCMATH_BS_KMAXX
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
//#define TINY 1.0e-30
#define SCALMX 0.1


void Miscmath::bsstep(double *y, double *dydx, const int nv, double *xx, const double htry, const double eps, double *yscal, double *hdid, double *hnext, moDerivs derivs, const Mathobject& mo) {

  int i,iq,k,kk,km=0;  //  initialize km for the sake of quieting warning. NRecipes do not initialize it
  int first=1,kmax,kopt;
  double epsold = -1.0,xnew;
  double eps1,errmax,fact,h,red=0.0,scale=0.0,work,wrkmin,xest;
  double a[IMAXX+1];
  double alf[KMAXX+1][KMAXX+1];
  int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
  int reduct,exitflag=0;
  
  double err[KMAXX+1];
  
  double *yerr = new double[nv+1];
  double *ysav = new double[nv+1];
  double *yseq = new double[nv+1];
  
  if (nv > MISCMATH_BS_MAXNVAR) throw Bad_Error("Miscmath::bsstep() nvar exceeds MISCMATH_BS_MAXNVAR, please recompile using higher value");
  
  if (eps != epsold)  {
    *hnext = xnew = -1.0e29;
    eps1=SAFE1*eps;
    a[1]=nseq[1]+1;
    for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
    for (iq=2;iq<=KMAXX;iq++) {
      for (k=1;k<iq;k++)
	alf[k][iq]=pow(eps1,(a[k+1]-a[iq+1])/
		       ((a[iq+1]-a[1]+1.0)*(2*k+1)));
    }
    epsold=eps;
    for (kopt=2;kopt<KMAXX;kopt++)
      if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
    kmax=kopt;
  }
  h=htry;
  for (i=1;i<=nv;i++) ysav[i]=y[i];
  if (*xx != xnew || h != (*hnext)) {
    first=1;
    kopt=kmax;
  }
  reduct=0;
  for (;;) {
    for (k=1;k<=kmax;k++) {
      xnew=(*xx)+h;
      if (xnew == (*xx)) throw Bad_Error("step size underflow in bsstep");
      mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq,derivs,mo);
      xest=SQR(h/nseq[k]);
      pzextr(k,xest,yseq,y,yerr,nv);
      if (k != 1) {
	errmax=TINY;
	for (i=1;i<=nv;i++) errmax=max(errmax,fabs(yerr[i]/yscal[i]));
	errmax /= eps;
	km=k-1;
	err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
      }
      if (k != 1 && (k >= kopt-1 || first)) {
	if (errmax < 1.0) {
	  exitflag=1;
	  break;
	}
	if (k == kmax || k == kopt+1) {
	  red=SAFE2/err[km];
	  break;
				}
	else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
	  red=1.0/err[km];
	  break;
	}
	else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
	  red=alf[km][kmax-1]*SAFE2/err[km];
	  break;
					}
	else if (alf[km][kopt] < err[km]) {
	  red=alf[km][kopt-1]/err[km];
	  break;
	}
			}
    }
    if (exitflag) break;
    red=min(red,REDMIN);
    red=max(red,REDMAX);
		h *= red;
		reduct=1;
  }
  *xx=xnew;
  *hdid=h;
  first=0;
  wrkmin=1.0e35;
  for (kk=1;kk<=km;kk++) {
    fact=max(err[kk],SCALMX);
    work=fact*a[kk+1];
    if (work < wrkmin) {
      scale=fact;
      wrkmin=work;
      kopt=kk+1;
    }
  }
  *hnext=h/scale;
  if (kopt >= k && kopt != kmax && !reduct) {
    fact=max(scale/alf[kopt-1][kopt],SCALMX);
    if (a[kopt+1]*fact <= wrkmin) {
      *hnext=h/fact;
      kopt++;
		}
  }
  
  delete[] yseq;
  delete[] ysav;
  delete[] yerr;
}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
//#undef TINY
#undef SCALMX
#undef NRANSI

/* note #undef's at end of file */
#define NRANSI

void Miscmath::mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep,
	double yout[], moDerivs derivs, const Mathobject &mo)
{
	int n,i;
	double x,swap,h2,h;

	double *ym = new double[nvar+1];
	double *yn = new double[nvar+1];

	h=htot/nstep;
	for (i=1;i<=nvar;i++) {
		ym[i]=y[i];
		yn[i]=y[i]+h*dydx[i];
	}
	x=xs+h;
	(mo.*derivs)(x,yn,yout);
	h2=2.0*h;
	for (n=2;n<=nstep;n++) {
		for (i=1;i<=nvar;i++) {
			swap=ym[i]+h2*yout[i];
			ym[i]=yn[i];
			yn[i]=swap;
		}
		x += h;
		(mo.*derivs)(x,yn,yout);
	}
	for (i=1;i<=nvar;i++)
		yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
	delete[] yn;
	delete[] ym;
}

void Miscmath::pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
	int k1,j;
	double q,f2,f1,delta;

	double *c = new double[nv+1];
	bs_x[iest]=xest;
	for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
	if (iest == 1) {
		for (j=1;j<=nv;j++) bs_d[j][1]=yest[j];
	} else {
		for (j=1;j<=nv;j++) c[j]=yest[j];
		for (k1=1;k1<iest;k1++) {
			delta=1.0/(bs_x[iest-k1]-xest);
			f1=xest*delta;
			f2=bs_x[iest-k1]*delta;
			for (j=1;j<=nv;j++) {
				q=bs_d[j][k1];
				bs_d[j][k1]=dy[j];
				delta=c[j]-q;
				dy[j]=f1*delta;
				c[j]=f2*delta;
				yz[j] += dy[j];
			}
		}
		for (j=1;j<=nv;j++) bs_d[j][iest]=dy[j];
	}
	delete[] c;
}
#undef NRANSI



void Miscmath::rzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv) {
	int k,j;
	double yy,v,ddy=0,c,b1,b;

	double *fx = new double[iest+1];
	bs_x[iest]=xest;
	if (iest == 1)
	  for (j=1;j<=nv;j++) {
	    yz[j]=yest[j];
	    bs_d[j][1]=yest[j];
	    dy[j]=yest[j];
	  }
	else {
	  for (k=1;k<iest;k++)
	    fx[k+1]=bs_x[iest-k]/xest;
	  for (j=1;j<=nv;j++) {
	    v=bs_d[j][1];
	    bs_d[j][1]=yy=c=yest[j];
	    for (k=2;k<=iest;k++) {
	      b1=fx[k]*v;
	      b=b1-c;
	      if (b) {
		b=(c-v)/b;
		ddy=c*b;
		c=b1*b;
	      } else
		ddy=v;
	      if (k != iest) v=bs_d[j][k];
	      bs_d[j][k]=ddy;
	      yy += ddy;
	    }
	    dy[j]=ddy;
	    yz[j]=yy;
	  }
	}
	delete[] fx;
}



double Miscmath::bs_d[MISCMATH_BS_MAXNVAR+1][MISCMATH_BS_KMAXX+1];
double Miscmath::bs_x[MISCMATH_BS_KMAXX+1];
double Miscmath::bs_1d[MISCMATH_BS_KMAXX+1];



//#define TINY 1.0e-30
#define MAXSTP 10000
double  Miscmath::oneDimOdeint(const double x1, const double x2, const double eps, moSingle derivs, const Mathobject& mo)  {
	int nstp;
	double x,hnext,hdid,h;

	double ystart = 0;

	double yscal,y, dydx;
	double h1 = fabs(1e-2*(x2 - x1));

	x=x1;
	h=SIGN(h1,x2-x1);
	y = ystart;
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		dydx = (mo.*derivs)(x);
		yscal = fabs(y)+fabs(dydx*h)+TINY;	
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;

		
		oneDimBsstep(&y,dydx,&x,h,eps, yscal ,&hdid,&hnext,derivs, mo);
		
		if ((x-x2)*(x2-x1) >= 0.0)   return y;
		
		if (fabs(hnext) <= 0.0)  throw Bad_Error("Step size too small in odeint");
		h=hnext;
	}
	throw Bad_Error("Too many steps in routine odeint");
}
#undef MAXSTP
//#undef TINY
#undef NRANSI

#undef NVAR               //! Undefining NVar



#define KMAXX MISCMATH_BS_KMAXX
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
//#define TINY 1.0e-30
#define SCALMX 0.1

void Miscmath::oneDimBsstep(double *y, const double dydx, double *xx, const double htry, const double eps, const double yscal, double *hdid, double *hnext, moSingle derivs, const Mathobject& mo) {

	int iq,k,kk,km=0;
	int first=1,kmax,kopt;
	double epsold = -1.0,xnew;
	double eps1,errmax,fact,h,red=0,scale=0,work,wrkmin,xest;
	double a[IMAXX+1];
	double alf[KMAXX+1][KMAXX+1];
	int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
	int reduct,exitflag=0;

	double err[KMAXX+1];

	double yerr,ysav,yseq;


	if (eps != epsold)  {
		*hnext = xnew = -1.0e29;
		eps1=SAFE1*eps;
		a[1]=nseq[1]+1;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=2;iq<=KMAXX;iq++) {
			for (k=1;k<iq;k++)
				alf[k][iq]=pow(eps1,(a[k+1]-a[iq+1])/
					((a[iq+1]-a[1]+1.0)*(2*k+1)));
		}
		epsold=eps;
		for (kopt=2;kopt<KMAXX;kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
	}
	h=htry;
	ysav = *y;

	if (*xx != xnew || h != (*hnext)) {
		first=1;
		kopt=kmax;
	}
	reduct=0;
	for (;;) {
		for (k=1;k<=kmax;k++) {
			xnew=(*xx)+h;
			if (xnew == (*xx)) throw Bad_Error("step size underflow in one dim  bsstep");
			yseq = oneDimMmid(ysav,dydx,*xx,h,nseq[k],derivs,mo);
			xest=SQR(h/nseq[k]);
			oneDimPzextr(k,xest,yseq,y,&yerr);
			if (k != 1) {
				errmax=TINY;
				errmax=max(errmax,fabs(yerr/yscal));
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
			}
			if (k != 1 && (k >= kopt-1 || first)) {
				if (errmax < 1.0) {
					exitflag=1;
					break;
				}
				if (k == kmax || k == kopt+1) {
					red=SAFE2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
						red=1.0/err[km];
						break;
					}
				else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
						red=alf[km][kmax-1]*SAFE2/err[km];
						break;
					}
				else if (alf[km][kopt] < err[km]) {
					red=alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
		if (exitflag) break;
		red=min(red,REDMIN);
		red=max(red,REDMAX);
		h *= red;
		reduct=1;
	}
	*xx=xnew;
	*hdid=h;
	first=0;
	wrkmin=1.0e35;
	for (kk=1;kk<=km;kk++) {
		fact=max(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin) {
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
		}
	}
	*hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct) {
		fact=max(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			*hnext=h/fact;
			kopt++;
		}
	}

}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
//#undef TINY
#undef SCALMX
#undef NRANSI
double Miscmath::oneDimMmid(const double y, const double dydx, const double xs, const double htot, const int nstep, moSingle derivs, const Mathobject &mo) {
	int n;
	double x,swap,h2,h;
	double ym,yn;


	h=htot/nstep;
	ym = y;
	yn =y +h*dydx;

	x=xs+h;
	double yout = (mo.*derivs)(x);
	h2=2.0*h;
	for (n=2;n<=nstep;n++) {
	  swap = ym+h2*yout;
	  ym=yn;
	  yn=swap;
	  x += h;
	  yout = (mo.*derivs)(x);
	}
	return 0.5*(ym+yn+h*yout);
}

void Miscmath::oneDimPzextr(const int iest,const double xest, const double yest, double* yz, double* dy)
{
	int k1;
	double q,f2,f1,delta,c;

	bs_x[iest]=xest;
	*dy=*yz=yest;
	if (iest == 1) bs_1d[1]=yest;
	else {
	  c =yest;
	  for (k1=1;k1<iest;k1++) {
	    delta=1.0/(bs_x[iest-k1]-xest);
	    f1=xest*delta;
	    f2=bs_x[iest-k1]*delta;

	    q=bs_1d[k1];
	    bs_1d[k1]= *dy;
	    delta=c -q;
	    *dy =f1*delta;
	    c =f2*delta;
	    *yz  += *dy;
	    
	  }
	  bs_1d[iest] = *dy;
	}
}
#undef NRANSI


//
//
// Finding zeros of a function
//
//


#include <math.h>
#define NRANSI
//#include "nrutil.h"
//#define ITMAX 100
#define ITMAX 100
#define EPS 3.0e-14

double Miscmath::zbrent(moSingle func, const Mathobject &mo,double x1, double x2, double tol, double ZeroLevel) {
	int iter;
	double a=x1,b=x2,c=x2,d=0,e=0,min1,min2;
	double fa= (mo.*func)(a) - ZeroLevel,fb=(mo.*func)(b) - ZeroLevel,fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
	  cout  << "zbrent(zerolevel= " << ZeroLevel << "): f(a=" << a << ") = "
                << fa << ", and f(b=" << b << ") = "  << fb << endl;
	  throw Bad_Error("Root must be bracketed in zbrent");
	}
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(mo.*func)(b) - ZeroLevel;
	}
	throw Bad_Error("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}
#undef ITMAX
#undef EPS
#undef NRANSI

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0
double Miscmath::dfridr(moSingle func,Mathobject& mo,double x, double h, double *err) {
  int i,j;
  double errt,fac,hh,ans=0;
  bool virgin = true;
  
  //cout << " in" << endl;

  if (h == 0.0) throw Bad_Error("h must be nonzero in dfridr.");
  double a[NTAB+1][NTAB+1];  
  hh=h;
  a[1][1]=( (mo.*func)(x+hh)-(mo.*func)(x-hh))/(2.0*hh);
  //cout << x << " h: " << hh << "  f(x+h)-f(x-h): " << (mo.*func)(x+hh)-(mo.*func)(x-hh) <<"   a[1][1]: " << a[1][1] << endl;
  *err=BIG;
  for (i=2;i<=NTAB;i++) {
    hh /= CON;
    a[1][i]=((mo.*func)(x+hh)-(mo.*func)(x-hh))/(2.0*hh);
    //cout << x << " h: " << hh << "  f(x+h)-f(x-h): " << (mo.*func)(x+hh)-(mo.*func)(x-hh) << "   a[1][i]: " << a[1][i] << endl;
    fac=CON2;
    for (j=2;j<=i;j++) {
      a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
      //cout << "extrapolation: a["<<j<<"]["<<i<<"]: " <<  a[j][i] << endl;
      fac=CON2*fac;
      errt=max(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
      //cout << "errt resulting: "<< errt << endl;
      if (errt <= *err || virgin) {
	*err=errt;
	ans=a[j][i];
	virgin = false;
	    }
    }
    if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
  }
  //if (virgin) throw Bad_Error("Miscmath::dfridr() derivation impossible");
  //cout << " out" << endl;
  //throw Bad_Error("stop");
  return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE



double Miscmath::relError(const double x, const double y) {
  double denom = min(fabs(x) , fabs(y));

  if (denom == 0.0) return 1.0e100;

  return fabs(x-y)/denom;
}

double Miscmath::sigma(const map<double,int>& m,double s) {

  int total=0,sum=0;
  for (map<double,int>::const_reverse_iterator i = m.rbegin(); i != m.rend(); i++) total+= i->second;
  double threshold = (1-s)*total;
  for (map<double,int>::const_reverse_iterator i = m.rbegin(); i != m.rend(); i++) {
    sum += i->second;
    if (sum > threshold) return i->first;
  }
  return 0;
  throw Bad_Error("Miscmath::sigma() reached end of map without crossing threshold.");
}


pair<double,double> Miscmath::Sigma(const map<double,int>& m,double s) {
  
  int sum = 0,total=0;
  for (map<double,int>::const_iterator i = m.begin() ; i != m.end(); i ++) total += i->second;

  double pre = 0.0,left = 0.0,right = 0.0;
  bool left_virgin=true;
  cout << "total : " << total << endl;
  for (map<double,int>::const_iterator i = m.begin(); i != m.end();i++) {
    sum += i->second;
    if (sum > (1-s)*total && left_virgin) { 
      if (i == m.begin()) left = i->first; else left = (pre + i->first)*0.5; 
      left_virgin = false;
    }
    if (sum > s*total) {
      if (i == m.begin()) right= i->first; else right = (pre + i->first)*0.5;
      break;
    }
    pre = i->first;
  }
  return pair<double,double>(left,right);

}


/*!
  Built a histogram.

  Returns a map<double,int> histogram where histogram[x] contains the 
  number of times a value in the list<double> l has been in some bin around x.

  The histogram runs from a to b and has k bins. 

  If you ommit k, then bins will be used such that in average, 200 values are in 
  the bins. If this is impossible, 10 bins will be used.

*/
  
map<double,int> Miscmath::histogram(list<double>& l,double a,double b,int k) {
  if (k == 0) k = l.size()/200;
  if (k == 0) k = 10;
 
  map<double,int> m;
  double delta = (b-a)/k;


  for (int q=0;q<k;q++) m[a + (q+0.5) *delta] = 0;
  for (list<double>::iterator i = l.begin(); i != l.end(); i++) {
    double y = *i;
    if (y >= a && y <= b) {
      double pos = a + (ceil((y-a)/delta) - 0.5)*delta;
      m[pos] += 1;
    } 
  }
  return m;
}


map<double,int> Miscmath::histogram(list<double>& l, list<double>& bins) {
  double MaxRange = -1e100;
  map<double,int> m;
  map<double,double> range;
  //  map<double,double> ToCheck;

  for (list<double>::iterator i = bins.begin(); i != bins.end(); i++) {
    list<double>::iterator j = i;
    if (++j != bins.end()) { 
      double rng = 0.5*(*j - *i);
      m[*i] = 0;  // initialize      
      //      ToCheck[*i] = *i + rng; // mapping of lower bound -> lower bound
      range[*i] = rng;
      if (rng > MaxRange) MaxRange = rng;
      cout << "Central point at: "<< *i + 0.5*(*j - *i) << endl;
    }
  }


  for (list<double>::iterator i = l.begin(); i != l.end(); i++) {
    double x = *i;
    map<double,int>::iterator best = m.lower_bound(x);
    
    if (best == m.end()) {
      cout << "OH WAS SOLL'n DAS: x = " <<x  << endl;
      throw Bad_Error("halt");
    } else {
      if (best != m.begin()) best--;
      double low = best->first;
      // cout << "I bin x: " << x <<  " to : [ " << best->first << " : " << best->second << " ] " << endl;
      if (x < best->first) throw Bad_Error("no bracket low");
      if (x >(++best)->first) throw Bad_Error("no bracket up");
      m[low]++;
    }
  }

  // normalize to MaxRange
  for (map<double,int>::iterator i = m.begin(); i != m.end(); i++) {
    i->second = (int) rint( i->second * MaxRange/range[i->first]);
  }
  return m;
}


void Miscmath::printHistogram(const map<double,int>&m,ostream& o) {

  for (map<double,int>::const_iterator i=m.begin();i != m.end();i++) {
    o << i->first << "\t" << i->second << endl;
  }
    
}


/*!
  Returns number between 0 and 1 following the rule:
  between left  and right 
  it returns 1. From left on leftwards it returns a gaussian
  distribution value continously with standarddevfiation lgauss. The
  same applies to rgauss
*/
double Miscmath::topHatGaussian(double x,double left,double right,double lgauss, double rgauss) {
  if (x >= left && x <=  right) return 1.0;

  if (x < left) return exp(-(x-left)*(x-left) / (2*lgauss*lgauss));
  if (x > right) return exp(-(x-right)*(x-right) / (2*rgauss*rgauss));
  throw Bad_Error("Miscmath::topHatGaussian() can't happen");
}



double Miscmath::fractionCorrespondingToSigma(double sigma) {
  double regio = rombint((moSingle)&Miscmath::fctsint, this, -sigma, sigma, 1e-3);
  double all =  rombint((moSingle)&Miscmath::fctsint, this, -20, 20, 1e-3);
  return regio / all;
}

double Miscmath::fctsint(const double x) const {   return exp(-x*x/2.0); }

double Miscmath::DeltaChi2CorrespondingToSigma(double nu, double sigma) {
  double p = fractionCorrespondingToSigma(sigma);
  DeltaChi2_Nu =  nu;
  return zbrent((moSingle)&Miscmath::DeltaChi2Helper, *this,0,100,1e-5,1.0-p);
}

double Miscmath::DeltaChi2Helper(double chi2) {
  return GammaQ(DeltaChi2_Nu*0.5, chi2*0.5);
}



/*! 
  Quick and dirty incomplete Gamma function
*/
double Miscmath::GammaQ(double a, double x) {
  Gamma_z = a;
  double GammaOfA = Miscmath::rombint((moSingle)&Miscmath::GammaInt,this,0,300,1e-5);
  double Second =  Miscmath::rombint((moSingle)&Miscmath::GammaInt,this,0,x,1e-5);
  return 1.0 - Second/GammaOfA;
}

double Miscmath::GammaInt(double t) {
  if (t > 1.0) return exp( (Gamma_z - 1.0)*log(t) - t);
  return pow(t,Gamma_z-1.0)*exp(-t);
}

/*!
  return min(denom,-epsilon) if denom < 0
  and  max(denom,epsilon) if denom >= 0
  
  So if a quantity diverges as 1 / denom then this
  will cut off at 1/epsilon
*/

double Miscmath::removePole(double denom, double epsilon) {
  if (denom < 0) return min(denom,-epsilon);
  else return max(denom,epsilon);
}


double Miscmath::logAdd(double x, double y)
{
  if (x<y) {
    return y+log(1.+exp(x-y));
  }
  return x+log(1.+exp(y-x));
}

