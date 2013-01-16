   
#include "mynrutils.h"
#include "myutils.h"
namespace NR {

#include <math.h>
#define NRANSI

double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
	else return (DEQUAL(absb,0.0) ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software !(^$S_!~. */

}
