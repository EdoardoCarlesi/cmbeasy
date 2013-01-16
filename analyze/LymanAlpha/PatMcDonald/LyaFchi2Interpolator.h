#ifndef LYAFCHI2INTERPOLATOR_FLAG
#define LYAFCHI2INTERPOLATOR_FLAG

#include "CPP.h"
#include "GlassInterpolator.h"

class LyaFchi2Interpolator {
 
  private:

    GlassInterpolator* gi;
    static const double lkcent; //log(0.009);
    static const double lkstep; //log(1.1);
    string tablefile;

  public:
   
    LyaFchi2Interpolator(string tablefile="default",double scale0=(0.43625-0.418642)*100,
                          double scale1=(2.38504-2.36142)*100,double scale2=0.03333333*100);
    double chi2_old(int N, double *k, double* D2);
    double chi2(double D2, double neff, double dndlnk); //!< chi2 of LYA data 
    ~LyaFchi2Interpolator();
    
};

#endif

