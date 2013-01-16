#include <math.h>
#include "bessel.h" 
#include <string>
#include <fstream>
#include <iostream>

Bessel::Bessel() {}



/*
double Bessel::spherical(int l,double x) {
  double sj,sy,sjp,syp;
  sphbes(l,x,&sj,&sy,&sjp,&syp); 
  return sj;
}
*/


double Bessel::bessj1(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
			+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
			+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}


 
 
double Bessel::bessj0(double x) 
{
	double ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			-y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}



 

#define ACC 40.0       // 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

double Bessel::bessj(int n, double x)
{
 
	int j,jsum,m;
	double ax,bj,bjm,bjp,sum,tox,ans;
	

	if (n < 2) throw Bad_Error("Index n less than 2 in bessj");
	ax=fabs(x);
	if (ax == 0.0)
		return 0.0;
	else if (ax > (double) n) {
		tox=2.0/ax;
		bjm=bessj0(ax);
		bj=bessj1(ax);
		for (j=1;j<n;j++) {
			bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
		ans=bj;
	} else {
		tox=2.0/ax;
		m=2*((n+(int) sqrt(ACC*n))/2);
		jsum=0;
		bjp=ans=sum=0.0;
		bj=1.0;
		for (j=m;j>0;j--) {
			bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			if (fabs(bj) > BIGNO) {
				bj *= BIGNI;
				bjp *= BIGNI;
				ans *= BIGNI;
				sum *= BIGNI;
			}
			if (jsum) sum += bj;
			jsum=!jsum;
			if (j == n) ans=bjp;
		}
		sum=2.0*sum-bj;
		ans /= sum;
	}
	return x < 0.0 && (n & 1) ? -ans : ans;
}
#undef ACC
#undef BIGNO
#undef BIGNI



double Bessel::bessy0(double x)
{
	double z;
	double xx,y,ans,ans1,ans2;

	if (x < 8.0) {
		y=x*x;
		ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
			+y*(10879881.29+y*(-86327.92757+y*228.4622733))));
		ans2=40076544269.0+y*(745249964.8+y*(7189466.438
			+y*(47447.26470+y*(226.1030244+y*1.0))));
		ans=(ans1/ans2)+0.636619772*bessj0(x)*log(x);
	} else {
		z=8.0/x;
		y=z*z;
		xx=x-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			+y*(-0.934945152e-7))));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}



double Bessel::bessy1(double x)
{
	double z;
	double xx,y,ans,ans1,ans2;

	if (x < 8.0) {
		y=x*x;
		ans1=x*(-0.4900604943e13+y*(0.1275274390e13
			+y*(-0.5153438139e11+y*(0.7349264551e9
			+y*(-0.4237922726e7+y*0.8511937935e4)))));
		ans2=0.2499580570e14+y*(0.4244419664e12
			+y*(0.3733650367e10+y*(0.2245904002e8
			+y*(0.1020426050e6+y*(0.3549632885e3+y)))));
		ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
	} else {
		z=8.0/x;
		y=z*z;
		xx=x-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}



double Bessel::bessy(int n, double x)
{
	int j;
	double by,bym,byp,tox;

	if (n < 2) throw Bad_Error("Index n less than 2 in bessy");
	tox=2.0/x;
	by=bessy1(x);
	bym=bessy0(x);
	for (j=1;j<n;j++) {
		byp=j*tox*by-bym;
		bym=by;
		by=byp;
	}
	return by;
}

double Bessel::chebev(double a, double b, double c[], int m, double x)
{
	double d=0.0,dd=0.0,sv,y,y2;
	int j;

	if ((x-a)*(x-b) > 0.0) throw Bad_Error("x not in range in routine chebev");
	y2=2.0*(y=(2.0*x-a-b)/(b-a));
	for (j=m-1;j>=1;j--) {
		sv=d;
		d=y2*d-dd+c[j];
		dd=sv;
	}
	return y*d-dd+0.5*c[0];
}


/* note #undef's at end of file */
#define NUSE1 5
#define NUSE2 5

void Bessel::beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
	double xx;
	static double c1[] = {
		-1.142022680371168e0,6.5165112670737e-3,
		3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
		3.67795e-11,-1.356e-13};
	static double c2[] = {
		1.843740587300905e0,-7.68528408447867e-2,
		1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
		2.423096e-10,-1.702e-13,-1.49e-15};

	xx=8.0*x*x-1.0;
	*gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
	*gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
	*gampl= *gam2-x*(*gam1);
	*gammi= *gam2+x*(*gam1);
}
#undef NUSE1
#undef NUSE2


#define EPS 1.0e-16
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
#define PI 3.141592653589793


#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

void Bessel::bessjy(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp)
{

	int i,isign,l,nl;
	double a,b,br,bi,c,cr,ci,d,del,del1,den,di,dlr,dli,dr,e,f,fact,fact2,
		fact3,ff,gam,gam1,gam2,gammi,gampl,h,p,pimu,pimu2,q,r,rjl,
		rjl1,rjmu,rjp1,rjpl,rjtemp,ry1,rymu,rymup,rytemp,sum,sum1,
		temp,w,x2,xi,xi2,xmu,xmu2;

	if (x <= 0.0 || xnu < 0.0) throw Bad_Error("bad arguments in bessjy");
	nl=(x < XMIN ? (int)(xnu+0.5) : IMAX(0,(int)(xnu-x+1.5)));
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	w=xi2/PI;
	isign=1;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=1;i<=MAXIT;i++) {
		b += xi2;
		d=b-d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b-1.0/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=c*d;
		h=del*h;
		if (d < 0.0) isign = -isign;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > MAXIT) throw Bad_Error("x too large in bessjy; try asymptotic expansion");
	rjl=isign*FPMIN;
	rjpl=h*rjl;
	rjl1=rjl;
	rjp1=rjpl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--) {
		rjtemp=fact*rjl+rjpl;
		fact -= xi;
		rjpl=fact*rjtemp-rjl;
		rjl=rjtemp;
	}
	if (rjl == 0.0) rjl=EPS;
	f=rjpl/rjl;
	if (x < XMIN) {
		x2=0.5*x;
		pimu=PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,&gam1,&gam2,&gampl,&gammi);
		ff=2.0/PI*fact*(gam1*cosh(e)+gam2*fact2*d);
		e=exp(e);
		p=e/(gampl*PI);
		q=1.0/(e*PI*gammi);
		pimu2=0.5*pimu;
		fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2)/pimu2);
		r=PI*pimu2*fact3*fact3;
		c=1.0;
		d = -x2*x2;
		sum=ff+r*q;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*(ff+r*q);
			sum += del;
			del1=c*p-i*del;
			sum1 += del1;
			if (fabs(del) < (1.0+fabs(sum))*EPS) break;
		}
		if (i > MAXIT) throw Bad_Error("bessy series failed to converge");
		rymu = -sum;
		ry1 = -sum1*xi2;
		rymup=xmu*xi*rymu-ry1;
		rjmu=w/(rymup-f*rymu);
	} else {
		a=0.25-xmu2;
		p = -0.5*xi;
		q=1.0;
		br=2.0*x;
		bi=2.0;
		fact=a*xi/(p*p+q*q);
		cr=br+q*fact;
		ci=bi+p*fact;
		den=br*br+bi*bi;
		dr=br/den;
		di = -bi/den;
		dlr=cr*dr-ci*di;
		dli=cr*di+ci*dr;
		temp=p*dlr-q*dli;
		q=p*dli+q*dlr;
		p=temp;
		for (i=2;i<=MAXIT;i++) {
			a += 2*(i-1);
			bi += 2.0;
			dr=a*dr+br;
			di=a*di+bi;
			if (fabs(dr)+fabs(di) < FPMIN) dr=FPMIN;
			fact=a/(cr*cr+ci*ci);
			cr=br+cr*fact;
			ci=bi-ci*fact;
			if (fabs(cr)+fabs(ci) < FPMIN) cr=FPMIN;
			den=dr*dr+di*di;
			dr /= den;
			di /= -den;
			dlr=cr*dr-ci*di;
			dli=cr*di+ci*dr;
			temp=p*dlr-q*dli;
			q=p*dli+q*dlr;
			p=temp;
			if (fabs(dlr-1.0)+fabs(dli) < EPS) break;
		}
		if (i > MAXIT) throw Bad_Error("cf2 failed in bessjy");
		gam=(p-f)/q;
		rjmu=sqrt(w/((p-f)*gam+q));
		rjmu=SIGN(rjmu,rjl);
		rymu=rjmu*gam;
		rymup=rymu*(p+q/gam);
		ry1=xmu*xi*rymu-rymup;
	}
	fact=rjmu/rjl;
	*rj=rjl1*fact;
	*rjp=rjp1*fact;
	for (i=1;i<=nl;i++) {
		rytemp=(xmu+i)*xi2*ry1-rymu;
		rymu=ry1;
		ry1=rytemp;
	}
	*ry=rymu;
	*ryp=xnu*xi*rymu-ry1;
}
#undef EPS
#undef FPMIN
#undef MAXIT
#undef XMIN
#undef PI

//#include <math.h>
#define RTPIO2 1.2533141

void Bessel::sphbes(int n, double x, double *sj, double *sy, double *sjp, double *syp)
{
 
	double factor,order,rj,rjp,ry,ryp;

	if (n < 0 || x <= 0.0) throw Bad_Error("bad arguments in sphbes");
	order=n+0.5;
	bessjy(x,order,&rj,&ry,&rjp,&ryp);
	factor=RTPIO2/sqrt(x);
	*sj=factor*rj;
	*sy=factor*ry;
	*sjp=factor*rjp-(*sj)/(2.0*x);
	*syp=factor*ryp-(*sy)/(2.0*x);
}
#undef RTPIO2

/*!
  Convenience wrapper for nrecipes bessjy(). Returns the fractional order
  Bessel Function J_nu(x)
*/
double Bessel::fractBessj(double nu, double x) {
  double rj, ry, rjp , ryp;
  bessjy(x, nu, &rj, &ry,&rjp, &ryp);
  return ry;
}


double Bessel::sphJ(int n, double x) {
  double sj,sy,sjp,syp;
  sphbes(n,x,&sj,&sy,&sjp,&syp);
  return sj;
}
  


/*!
  CMBFAST's version for calculating the spherical
  bessel function. Is nicer behaved than the one from
  Numerical recipes
*/

float Bessel::bjl(int l, float x) {
    /* System generated locals */
    double d__1;

    /* Local variables */
    static double beta, secb, cotb, beta2, sec2b, beta4, sec4b, beta6, 
	    sec6b, cot3b, cot6b, gamma1, prefactor, gamma2, deriv1, fl, ax, 
	    pi, cx, nu, sx, llimit, ulimit, ax2, rootpi, nu2, sx2, trigarg, 
	    trigcos, expterm, sum1, sum2, sum3, sum4, sum5;

    
    double jl=0; // result

/*  Calculates the spherical bessel function j_l(x) */
/*  and optionally its derivative for float x and int l>=0. */
/*  Asymptotic approximations 8.11.5, 8.12.5, and 8.42.7 from */
/*  G.N.Watson, A Treatise on the Theory of Bessel Functions, */
/*  2nd Edition (Cambridge University Press, 1944). */
/*  Higher terms in expansion for x near l given by */
/*  Airey in Phil. Mag. 31, 520 (1916). */
/*  This approximation is accurate to near 0.1% at the boundaries */
/*  between the asymptotic regions; well away from the boundaries */
/*  the accuracy is better than 10^{-5}. The derivative accuracy */
/*  is somewhat worse than the function accuracy but still better */
/*  than 1%. */
/*  Point *jlp initially to a negative value to forego calculating */
/*  the derivative; point it to a positive value to do the derivative */
/*  also (Note: give it a definite value before the calculation */
/*  so it's not pointing at junk.) The derivative calculation requires */
/*  only arithmetic operations, plus evaluation of one sin() for the */
/*  x>>l region. */
/*  Original code by Arthur Kosowsky   akosowsky@cfa.harvard.edu */
/*  This fortran version only computes j_l(x) */
    pi = (float)3.1415926536;
    rootpi = (float)1.772453851;
    gamma1 = (float)2.6789385347;
/* /+ Gamma function of 1/3 +/ */
    gamma2 = (float)1.3541179394;
/* /+ Gamma function of 2/3 +/ */
    ax = fabs(x);
    fl = (l);
    beta = pow(fl, .325);
    llimit = beta * 1.31;
/* /+ limits of asymptotic regions; fitted +/ */
    ulimit = beta * (float)1.48;
    nu = fl + (float).5;
    nu2 = nu * nu;
    if (l < 0) {
      cout << "Bessel function index < 0" << endl;
      throw Bad_Error("bjl() :: stop");
    }
    //         /+************ Use closed form for l<6 *********+/ 
    if (l < 6) {
	sx = sin(ax);
	cx = cos(ax);
	ax2 = ax * ax;
	if (l == 0) {
	    if (ax > (float).001) {
		jl = (float) (sx / ax);
	    } else {
		jl = (float) (1. - ax2 / 6.);
	    }
/*   /+ small x trap +/ */
	}
	if (l == 1) {
	    if (ax > (float).001) {
		jl = (float) ((sx / ax - cx) / ax);
	    } else {
		jl = (float) (ax / 3.);
	    }
	}
	if (l == 2) {
	    if (ax > (float).001) {
		jl = (float) ((cx * -3. / ax - sx * (1. - 3. / ax2)) / ax);
	    } else {
		jl = (float) (ax2 / 15.);
	    }
	}
	if (l == 3) {
	    if (ax > (float).001) {
		jl = (float) ((cx * (1. - 15. / ax2) - sx * (6. - 15. / ax2) /
			 ax) / ax);
	    } else {
		jl = (float) (ax * ax2 / 105.);
	    }
	}
	if (l == 4) {
	    if (ax > (float).001) {
		jl = (float) ((sx * (1. - 45. / (ax * ax) + 105. / (ax * ax * 
			ax * ax)) + cx * (10. - 105. / (ax * ax)) / ax) / ax);
	    } else {
		jl = (float) (ax2 * ax2 / 945.);
	    }
	}
	if (l == 5) {
	    if (ax > (float).001) {
		jl = (float) ((sx * (15. - 420. / (ax * ax) + 945. / (ax * ax 
			* ax * ax)) / ax - cx * ((float)1. - 105. / (ax * ax) 
			+ 945. / (ax * ax * ax * ax))) / ax);
	    } else {
		jl = (float) (ax2 * ax2 * ax / 10395.);
	    }
	}
/*          /+********************* x=0 *********************+/ */
    } else if (ax < 1e-30) {
	jl = (float)0.;
/*          /+************** Region 1: x << l ***************+/ */
    } else if (ax <= fl + (float).5 - llimit) {
      /*       beta=acosh(nu/ax) */
      if (nu / ax < 1.) cout << "trouble with acosh" << endl;
/* Computing 2nd power */
	d__1 = nu / ax;
	beta = log(nu / ax + sqrt(d__1 * d__1 - 1.));
/* (4.6.21) */
	cotb = nu / sqrt(nu * nu - ax * ax);
/* /+ cotb=coth(beta) +/ */
	cot3b = cotb * cotb * cotb;
	cot6b = cot3b * cot3b;
	secb = ax / nu;
	sec2b = secb * secb;
	sec4b = sec2b * sec2b;
	sec6b = sec4b * sec2b;
	sum1 = sec2b * (float)3. + (float)2.;
	expterm = sum1 * cot3b / (nu * (float)24.);
	sum2 = sec2b + (float)4.;
	expterm -= sum2 * sec2b * cot6b / (nu2 * (float)16.);
	sum3 = (float)16. - sec2b * (float)1512. - sec4b * (float)3654. - 
		sec6b * (float)375.;
	expterm -= sum3 * cot3b * cot6b / (nu * (float)5760. * nu2);
	sum4 = sec2b * (float)288. + (float)32. + sec4b * (float)232. + sec6b 
		* (float)13.;
	expterm -= sum4 * sec2b * cot6b * cot6b / (nu2 * (float)128. * nu2);
	expterm = exp(-nu * beta + nu / cotb - expterm);
	prefactor = sqrt(cotb / secb) / (nu * (float)2.);
	jl = (float) (prefactor * expterm);
/*          /+*************** Region 2: x >> l ***************+/ */
    } else if (ax >= fl + (float).5 + ulimit) {
	beta = acos(nu / ax);
	cotb = nu / sqrt(ax * ax - nu * nu);
/* /+ cotb=cot(beta) +/ */
	cot3b = cotb * cotb * cotb;
	cot6b = cot3b * cot3b;
	secb = ax / nu;
	sec2b = secb * secb;
	sec4b = sec2b * sec2b;
	sec6b = sec4b * sec2b;
	trigarg = nu / cotb - nu * beta - pi / (float)4.;
	sum1 = sec2b * (float)3. + (float)2.;
	trigarg -= sum1 * cot3b / (nu * (float)24.);
	sum3 = (float)16. - sec2b * (float)1512. - sec4b * (float)3654. - 
		sec6b * (float)375.;
	trigarg -= sum3 * cot3b * cot6b / (nu * (float)5760. * nu2);
	trigcos = cos(trigarg);
	sum2 = sec2b + (float)4.;
	expterm = sum2 * sec2b * cot6b / (nu2 * (float)16.);
	sum4 = sec2b * (float)288. + (float)32. + sec4b * (float)232. + sec6b 
		* (float)13.;
	expterm -= sum4 * sec2b * cot6b * cot6b / (nu2 * (float)128. * nu2);
	expterm = exp(-expterm);
	prefactor = sqrt(cotb / secb) / nu;
	jl = (float) (prefactor * expterm * trigcos);
/*          /+**************** Region 3: x near l ***************+/ */
    } else {
	beta = ax - nu;
	beta2 = beta * beta;
	beta4 = beta2 * beta2;
	beta6 = beta2 * beta4;
	sx = (float)6. / ax;
	sx2 = sx * sx;
	cx = sqrt(sx);
	secb = pow(sx, 0.333333333 );
	sec2b = secb * secb;
	deriv1 = gamma1 * secb;
	deriv1 += beta * gamma2 * sec2b;
	sum1 = (beta2 / (float)6. - (float).066666666666666666) * beta;
	deriv1 -= sum1 * sx * secb * gamma1 / (float)3.;
	sum2 = beta4 / (float)24. - beta2 / (float)24. + (float)
		.0035714285714285713;
	deriv1 -= sum2 * (float)2. * sx * sec2b * gamma2 / (float)3.;
	sum3 = beta6 / (float)720. - beta4 * (float)7. / (float)1440. + beta2 
		/ (float)288. - (float)2.7777777777777778e-4;
	deriv1 += sum3 * (float)4. * sx2 * secb * gamma1 / (float)9.;
	sum4 = (beta6 / (float)5040. - beta4 / (float)900. + beta2 * (float)
		19. / (float)12600. - (float)4.1269841269841269e-4) * beta;
	deriv1 += sum4 * (float)10. * sx2 * sec2b * gamma2 / (float)9.;
	sum5 = (beta4 * beta4 / (float)362880. - beta6 / (float)30240. + 
		beta4 * (float)71. / (float)604800. - beta2 * (float)121. / (
		float)907200. + (float)3.4095203738060878e-5) * beta;
	deriv1 -= sum5 * (float)28. * sx2 * sx * secb * gamma1 / (float)27.;
	jl = (float) (deriv1 * cx / (rootpi * (float)12.));
    }
    return jl;
} /* bjl_ */
