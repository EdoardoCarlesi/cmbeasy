#include "mynrutils.h"
#include "myutils.h"
namespace NR {

#define NRANSI

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
	int jj,j,i;
	double s,*tmp;

	tmp=dvector(1,n);
	for (j=1;j<=n;j++) {
		s=0.0;
		if (!DEQUAL(w[j],0.0)) {
			for (i=1;i<=m;i++) s += u[i][j]*b[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=1;j<=n;j++) {
		s=0.0;
		for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
		x[j]=s;
	}
	free_dvector(tmp,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software !(^$S_!~. */

}
