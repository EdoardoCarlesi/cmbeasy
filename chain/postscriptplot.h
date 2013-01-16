#ifndef POSTSCRIPTPLOT_H
#define  POSTSCRIPTPLOT_H

#include "lowlevelplot.h"
#include <string>
#include <fstream>

using namespace std;

class PostscriptPlot : public LowLevelPlot {
  
 public:
  
  string Resource;
  //PostscriptPlot(string resourcefile);

  PostscriptPlot(const LowLevelPlot&);

  void setResourceFile(string s) { Resource = s;}
  void plot(ofstream &out);
  string parseGraceString(string s,bool show,float LabelSize, float XStauch,bool Landscape=true);
  string parseTexLabel(string);
  int mac(float,float=10000);
  static void InterpolateColor(Color &c,Color &up, Color& down, float z, float Min,float Max);  //!< color c such that given value z within range Min,Max it is up for 1 and down for 0
};

#endif
