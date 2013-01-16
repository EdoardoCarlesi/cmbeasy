#ifndef CROSSOVER_H
#define CROSSOVER_H

#include "anchor.h"
#include "arbitrary.h"
#include "mathobject.h"
#include <cmath>
#include "spline.h"

/*!
  Obsolete Crossover-Quintessence class.
  Use CrossoverField instead.
*/
class Crossover : public Arbitrary {

 private:
  Anchor CrossAnchor;
  Spline *CW, *CW1, *CW2;
  double a_catch;

 public:
  
  Crossover(QuintCosmos &c);
 
  double w0() const { return param[0];}
  double w_ls() const { return param[1];}
  double A() const { return param[2];}
  double x_ls() const { return log(1.0/1100.0);}


  Type type() { return crossover;}
  void prepare();
  
  virtual double w(double a) const  {
    return CW->fastY(a);
  } 

 
  double CrossW(double a, bool catchit) {
    if (a <= a_catch && catchit) a = a_catch;
    if (a > 1.0 && catchit) a = 1.0;
    double x = log(a);
    return w0() + 2*x/x_ls()*(w_ls() - w0()) + A()*x*(3*x - 2*x_ls());
  }  
  
  virtual double dwda(double a) const {
    return CW1->fastY(a);
  }
  virtual double d2wda2(double a) const {
    return CW2->fastY(a);
  }
};


#endif 
