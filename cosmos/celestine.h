#ifndef CELESTINE_H
#define CELESTINE_H

#include "anchor.h"
#include "arbitrary.h"
#include "mathobject.h"
#include <cmath>
#include "spline.h"

/*!
  Arbitrary quintessence with the following functional form:

  w(a) = w_0 + w_1 (1 - a)
  

*/
  


class Celestine : public Arbitrary {

 public:
  
  Celestine(QuintCosmos &c) : Arbitrary(c) {
    param.resize(2);
    param[0] = -0.95;
    param[1] = 0;
  }
  

  virtual double w(double a) const  { return param[0] + param[1]*(1-a);}  //!< equation of state as function of scale factor
  virtual double dwda(double a) const { return -param[1];} //!< first  derivative wrt a:  d w(a)/ da
  virtual double d2wda2(double a) const { return 0;}  //!< second derivative wrt a:  d^2 w(a)/ da^2
  
  virtual Type type() { return celestine;}

};


#endif 
