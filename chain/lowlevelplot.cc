#include "lowlevelplot.h"
#include <math.h>
#include "global.h"
#include <iostream>
#include <sstream>
#include <string>
#include <stdlib.h>

#ifdef MACOSX_PANTHER
extern "C" int isnan(double); 
extern "C" int isinf(double); 
#endif 


LowLevelPlot::LowLevelPlot(float w,float h, float lft, float btm, float rght, float top) : Width(w), Height(h), LeftMargin(lft), BottomMargin(btm), RightMargin(rght), TopMargin(top)  {
  StandardLineWidth = 2e-3;
  MajorTickLength = 1e-2;

  FrameLineWidth = 100;
  CurveLineWidth = 100;
  TicksLineWidth = 100;

  XLabelOffset = 0.06;
  YLabelOffset = 0.11;

  TickLabelSize = 100;
  AxisLabelSize = 100;
  XTickLabelStyle = decimal;
  YTickLabelStyle = decimal;

  Width = 80;
  Height = 100;

}

void LowLevelPlot::prepare() {
  
  InvWidth = 1.0/Width;
  InvHeight = 1.0/Height;

  InvSkala = 0.5*(InvWidth + InvHeight);

  //  ScaleX2View = 1.0 - XShiftLow - XShiftUp;
  //  ScaleY2View = 1.0 - YShiftLow - YShiftUp;

  ScaleX2View = 1.0 - LeftMargin - RightMargin;
  ScaleY2View = 1.0 - BottomMargin - TopMargin;
  
}

string LowLevelPlot::int2string(int x) {
  std::ostringstream o;
  if (o << x) return o.str();
  return "nan";
}


string LowLevelPlot::float2string(float x, int significant) {  return number2string<float>(x,significant); }
string LowLevelPlot::double2string(double x, int significant) {  return number2string<double>(x,significant); }


template<class T> string LowLevelPlot::number2string(T x,int significant) {
  if (rint(x) == x) return int2string((int)rint(x));
  std::ostringstream o;
  if (fabs(x) >= 0.99999*pow(10.0,-significant))  o.setf(ios::fixed); else o.setf(ios::scientific);
  o.precision(significant);
  if (x == 0) return "0";
  if (o << x) return o.str();
  return "nan";
} 

double LowLevelPlot::splitExponent(double x) {
  double e = log10(fabs(x));
  return floor(e);
}

double LowLevelPlot::splitPrefactor(double x) {
  double e = splitExponent(x);
  double c = pow(10.0,e);
  return x / c;
}


TickLabel LowLevelPlot::createLabel(double x, double axis,LabelStyle style,int significant) {
  string d,e;
  switch (style) {
  case exponent: 
    d = double2string(splitPrefactor(x),significant);
    //cout << " Split prefactor   " << double2string(splitPrefactor(x)) << "  SplitExponent: " <<  double2string(splitExponent(x)) << endl;
    e = double2string(splitExponent(x));
    return TickLabel(axis,d,e);   
  default: 
    d = double2string(x,significant);
    return TickLabel(axis,d); 
  } 
}

string LowLevelPlot::createGraceLabel(double x,LabelStyle style,int significant) {
  string d,e="";
  switch (style) {
  case exponent: 
    if (x == 0) {
      d = "0";
    } else {
      d = double2string(splitPrefactor(x),significant);
      e = double2string(splitExponent(x));
      if (d == "1") d = ""; else d+= "\\S.\\N";
    }
    break;
  default: 
    d = double2string(x,significant);
  } 
  if (e != "") {
    d += "10\\S " + e;
  }
  return d;
}
