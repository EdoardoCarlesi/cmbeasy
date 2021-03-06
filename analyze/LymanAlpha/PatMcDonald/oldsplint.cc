
#include "mynrutils.h"

namespace NR {

void oldsplint(double* xa, double* ya, double* y2a, int n, double x, double *y)
{
        int klo,khi,k;
        double h,b,a;
                                                                                
        klo=1;
        khi=n;
        while (khi-klo > 1) {
                k=(khi+klo) >> 1;
                if (xa[k] > x) khi=k;
                else klo=k;
        }
        h=xa[khi]-xa[klo];
        if (h <= 0.0) {
          printf("Bad xa input to routine splint\n");
          exit(1);
        }
        a=(xa[khi]-x)/h;
        b=(x-xa[klo])/h;
        *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

}
