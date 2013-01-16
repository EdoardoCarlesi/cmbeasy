#include "miscmath.h"
#include "bessel.h"
#include "spline.h"
#include "controlpanel.h"
Miscmath miscmath;

int main()  {
/*     Calculates tables for computing Bessel functions in Lesing code. */  
  double lmax = 3500;
  lmax += 300;
  double xmax = M_PI*lmax/10;
  double dx =0.05;
  Spline b0(20000, "bessel0");
  Spline b2(20000, "bessel2");
  Spline b4(20000, "bessel4");
  Spline b6(20000, "bessel6");
  for (double x = 0; x <= xmax; x += dx) {
    b0.set(x, Bessel::bessj0(x));
    b2.set(x, Bessel::bessj(2,x));
    b4.set(x, Bessel::bessj(4,x));
    b6.set(x, Bessel::bessj(6,x));      
  }
  ofstream file(ControlPanel::cmbeasyDir("/resources/jlens.dat").c_str());  
  b0.save(file);
  b2.save(file);
  b4.save(file);
  b6.save(file);
}
