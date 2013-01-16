
#ifndef MYUTILS_FLAG
#define MYUTILS_FLAG

#include <stdio.h>
#include <stdlib.h>
//#include <limits>
#include "mynrutils.h"
#include "global.h"

using namespace std;

namespace MU {
  #define mytinyfloat 1.0e-7
  #define mytinydouble 1.0e-16
  #define DEQUAL(a,b) (fabs((a)-(b))<mytinydouble)
  #define FEQUAL(a,b) (fabs((a)-(b))<mytinyfloat)
//  #define DEQUAL(a,b) (fabs((a)-(b))<numeric_limits<double>::epsilon())
//  #define FEQUAL(a,b) (fabs((a)-(b))<numeric_limits<float>::epsilon())
  #define EQUAL(a,b,c) (fabs((a)-(b))<(c))

  #define PI (4.0*atan(1.0))
  #define cmpermeter (100.0)
  #define kmpermeter (0.001)
  #define gramperkg (1000.0)
  #define Angstrompermeter (1.0e10)

  #define radperdegree (2.0*PI/360.0)
  #define degreeperrad (1.0/radperdegree)
  #define radperarcmin (radperdegree/60.0)
  #define radperarcsec (radperarcmin/60.0)

  #define hbar_mks (1.05457266e-34) 
  #define G_mks (6.67259e-11)
  #define c_mks (2.99792458e8) 
  #define c_kms (c_mks*kmpermeter) 
  #define GeVperJoule (6.24150636e9)
  #define Mpcpermeter (3.24075574e-23)
  #define lywav_Ang (1215.67)
  #define lywav_mks (lywav_Ang/Angstrompermeter)
  #define kB_mks (1.380658e-23) //Joule/Kelvin 
  #define blackbodyconstant_mks (PI*PI*pow(kB_mks,4)/15/pow(hbar_mks*c_mks,3))
  #define Riemannzeta (1.202)

  #define H0overh (100.0/kmpermeter*Mpcpermeter) // inverse seconds
  #define LPlanck_mks (sqrt(hbar_mks*G_mks/pow(c_mks,3)))
  #define TPlanck_mks (LPlanck_mks/c_mks)
  #define MPlanck_mks (sqrt(hbar_mks*c_mks/G_mks))
  #define EPlanck_mks (MPlanck_mks*c_mks*c_mks)
  #define EDPlanck_mks (EPlanck_mks/pow(LPlanck_mks,3))

  #define inverseGeVpersecond (1.0/(hbar_mks*GeVperJoule))
  #define inverseGeVpermeter (inverseGeVpersecond/c_mks)
  #define EPlanck_GeV (EPlanck_mks*GeVperJoule)
  #define yearspersecond (1.0/(365.25*24*60*60))
  #define kgperproton (1.6726231e-27)
  #define kgperelectron (9.1093897e-31) 
  void reversevectororder(double *v,int N);
  void myexit(char error_text[]);

}

#endif

