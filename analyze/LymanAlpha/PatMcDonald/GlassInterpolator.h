#ifndef GLASSINTERPOLATOR_FLAG
#define GLASSINTERPOLATOR_FLAG

#include "CPP.h"
#include "myutils.h"
#include "mynrutils.h"
#include "ObjectInSpace.h"

struct GlassInterpolatorOutOfBound {};

class GlassInterpolator {
 
  private:

  vector<ObjectInSpace> voP;
  vector<double> scales;
  double SVDtol;
  int NDim;
  int Nd;
  
  public:
   
    GlassInterpolator(vector<ObjectInSpace> voP,const vector<double> scales,                      double SVDtol);
    void Interpolate(ObjectInSpace* P);
    ~GlassInterpolator();
    
};

#endif

