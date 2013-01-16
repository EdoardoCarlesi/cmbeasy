#ifndef LOWLEVELPLOT_H
#define  LOWLEVELPLOT_H

#include <vector>
#include <list> 
#include <string>
#include <map>
#include <cmath>
//#include <iostream>

using namespace std;
/*!
  Small helper for keeping track
  of the coordinate boundaries
*/
struct CoordRect {
  double x,X,y,Y; //!< min (lower case)  and max (upper case) coordinates
  CoordRect(double minX,double maxX, double minY, double maxY) : x(minX), X(maxX), y(minY), Y(maxY) {}
};

struct CoordPoint {
  float  x,y;
  CoordPoint() : x(0), y(0) {}
  CoordPoint(float X,float Y) : x(X), y(Y) {}
  CoordPoint(pair<float,float> p) : x(p.first), y(p.second) {}
  CoordPoint(const CoordPoint& c) : x(c.x), y(c.y) {} 
};

inline CoordPoint operator+(const CoordPoint& p1, const CoordPoint& p2) {
  return CoordPoint(p1.x+p2.x, p1.y+p2.y);
}

inline CoordPoint operator-(const CoordPoint& p1, const CoordPoint& p2) {
  return CoordPoint(p1.x-p2.x, p1.y-p2.y);
}

inline bool operator==(const CoordPoint& p1, const CoordPoint& p2) {
  return (p1.x==p2.x && p1.y==p2.y);
}


struct Color {
  float r,g,b;
  Color(float R=0, float G=0, float B=0) : r(R), g(G), b(B) {}
  Color(vector<float> c) : r(c[0]), g(c[1]), b(c[2]) {}
  void multiply(float f) { r *= f; g*=f ; b*=f;}
};

struct Rectangle {
  CoordRect C;
  // vector<float> FillColor;  // r.g.b.
  Color FillColor;
  Rectangle(double x,double X,double y, double Y) : C(x,X,y,Y)  {};
};

struct Curve {
  list<CoordPoint> Points;
  Color color; 
  int dash;
  /*! Empty Curve constructor */
  Curve() : dash(0) {}
  /*! Convenience constructor for one-straight-line-only curves. Takes a two pairs 
    of numbers (x1,y1) (x2,y2) and fills the Points list with these two points */
  Curve(pair<double,double> p1, pair<double,double> p2) : dash(0) { Points.push_back(CoordPoint(p1.first,p1.second)); Points.push_back(CoordPoint(p2.first,p2.second)); } 
  void thinOut(float ruler);
};

/*!
  Thins out a curve, i.e. drops points that are not needed for a smooth appearance.
  Parameter ruler sets the length scale on which it should be smooth. It drops points,
  if the remaining points suffice to approximate the y-position of each intermediate (dropped)
  point to better than "ruler". 
*/
inline void Curve::thinOut(float ruler) {
  float ruler2 = ruler*ruler;
  if (Points.size() < 5) return; // sensless for small curves, caution: algorithm crashes for < 4 
  list<CoordPoint> cpy;
  list<CoordPoint>::iterator i=Points.begin(); 
  cpy.push_back(*i);  // first point always necessary
  list<CoordPoint>::iterator j,nxt,sneak,right;
  nxt = i; nxt++;
  right = nxt; right++; // initially, i,nxt,right are next to each other
  sneak=right; sneak++; // sneak is one additional to the right
  for (;;) {
    bool drop=true; // can we drop the point ? 
    float delta_x =  (sneak->x- i->x);
    float delta_y = (sneak->y - i->y);
    for (j = nxt; j != sneak; j++) { // for all points between j+1 and right ...
      if (delta_x != 0) {
	float m = delta_y/delta_x;	
	float ymid = i->y + m*(j->x - i->x);	
	//cout << "m: " << m << "  ymid:  " << ymid << "  really it is at: " << j->y << endl;
	if ( (ymid - j->y)*(ymid-j->y) > ruler2 ) { // midpoint is necessary
	  //cout << "( " << j->x << " , " << j->y << " ) differs by " << fabs(ymid - j->y) << "  in between: " <<  "( " << i->x << " , " << i->y << " )  and ( " << sneak->x << " , " << sneak->y << " ) ";
	  drop = false; 
	  break;
	}
      } else drop = false;
    }
    // so we have checked all intermediate points
    if (drop) {
      right++; sneak++; // increase sneak and right
      if (sneak == Points.end()) { // if sneak is out of range, we're done
	cpy.push_back(*right);
	break;
      }
    } else {  // don't drop!
      cpy.push_back(*right);
      i=right;
      nxt=sneak;
      right=nxt;
      right++;
      if (right == Points.end()) break;
      sneak = right;
      sneak++; 
      if (sneak == Points.end()) {
	cpy.push_back(*right);
	break;	
      }
    }
  }
  Points = cpy;
}

/*!
  Points that can have a z-information for
  color encoding. z runs from 0..1 is relative to Min and Max, i.e.
  
  value = z*(Max-Min) + Min
   
*/
struct ColorPoints {
  //! a vector of n points, each having x,y,z information 
  vector< vector<float> > v;
  float Min,Max; //!< the minimum and Maximum value for which z gives the relative value
  float Radius; //!< Radius in units 0..1 of each point
  CoordPoint LeftLower; //!< LeftLower Corner in 0..1 units for the color scale
  float Width,Height; //!< Width and Height in 0..1 units for the color scale
  Color MinColor, MiddleColor, MaxColor; //!< Colors the Min and Max and 0.5 value
};

/*!
  Elementary confidence region.
 */
struct ConfidenceRegion {
  Color FillColor;
  Color PenColor;
  bool Fill;
  list<CoordPoint>  region;
  ConfidenceRegion() : Fill(true) {};	       
};


struct Block {
  float Width;
  Color FillColor;
  Color PenColor;
  bool UsePen;
  map<float, list<CoordPoint> > lines;
  Block() : UsePen(false) {};
  //  Block() : lines(0) {};
  //~Block() { if (lines) delete lines; }
};

struct TickLabel {
  float axis; //   position on axis (either x or y)
  float size; // (0..1)
  float offset; // the offset perpendicular to the axis
  string decimal;
  string exponent;
  string grace; 
  TickLabel(float a,string d="", string e="", string g="",float s=0.05) : axis(a), size(s), offset(0.01), decimal(d), exponent(e), grace(g) {}
};

struct TickMark {
  float axis;  // position on axis (either x or y)
  // float other;
  TickMark(float a) : axis(a) {}
};

struct AxisLabel { 
  float axis;
  float size;
  float offset;
  string text;
  AxisLabel(string t="", float a=0.5, float o=0.08, float s=0.08) : axis(a), size(s), offset(o), text(t) {}
} ;

struct SpecialPoint {
  enum  Type {Star, Circle} type;
  CoordPoint position;
  float size;
  SpecialPoint( Type t, CoordPoint p, float s ): type( t ), position( p ), size( s ) {}
};

struct ThingsToPlot {
  list<Rectangle> Rectangles;
  list<Curve> Curves;
  list<ConfidenceRegion> Regions;
  list<Block> Blocks;
  list<ColorPoints> CPs;
  list<SpecialPoint> SpecialPoints;
  list<TickMark> XTickMarks;
  list<TickMark> YTickMarks;
  list<TickLabel> XTickLabels;
  list<TickLabel> YTickLabels;
  AxisLabel XAxisLabel, YAxisLabel;
};

/*!
  Small struct to store data for axis label parsing
*/
struct GraceLabelHelper {
  double scale, shift;
  string family;
  string text;
  bool stopwagon;
};



class LowLevelPlot {
  
 public:
  enum LabelStyle {decimal, exponent};
  LowLevelPlot(float w=100,float h=70, float xsl=0.15, float ysl=0.15, float xsu=0.05, float ysu=0.05);

  void prepare();
  double splitExponent(double x);
  double splitPrefactor(double x);
  template<class T> string number2string(T x, int significant=4);
  string double2string(double x,int significant=4); 
  string float2string(float x, int significant=4); 

  string int2string(int x);
  TickLabel  createLabel(double x, double axis, LabelStyle style,int significant );
  string createGraceLabel(double x,LabelStyle stylem,int significant);
  ThingsToPlot Things;

  /*  ===========================================
      == These variabels are the ones that should be specified
      == before printing
      =========================================== */
  float Width, Height;
  float LeftMargin,BottomMargin,RightMargin,TopMargin;
  
  float StandardLineWidth; //!< overall linewidth (in a 0..1 range, 1 filling the view)
  int FrameLineWidth, CurveLineWidth, TicksLineWidth; // !< relative width (in %) w.r.t StandardLineWidth
  float MajorTickLength; //! Length of Tick-Mark
  
  float Aspect, InvAspect;
  
  float InvWidth, InvHeight;  // 1/width, 1/height
  float InvSkala; // mean of InvWidth, InvHeight


/*  ===========================================
    == The following variables have "control" character: They
    == are not needed for this class, but we use this class to  
    == store them
    =========================================== */
  
  int TickLabelSize; //!< Relative (in %) to some sensible size 
  LabelStyle XTickLabelStyle, YTickLabelStyle;
  bool AutomaticTick;
  float StartTick_x, StartTick_y;
  float StepTick_x, StepTick_y;  
  int Significant_x, Significant_y;
  map<float,string> XTickMarkList, YTickMarkList;  
  
  int AxisLabelSize; //!< Relative (in %) to some sensible size
  float XLabelOffset, YLabelOffset;
  string XLabelText,YLabelText; 

  string Name; 
  
 protected:
  

  // the sheet has coordinate 0..1
  // but the view (which) is smaller will
  // also have coordinate 0..1
  // the variables below will scale from the world 
  // to the view 

  float ScaleX2View,ScaleY2View;
  //float ScaleXFromView, ScaleYFromView;  // inverse variables
  


 
};


/*!
  Small struct controlling some parameters for postscript output of
  2-dimensional marginalized likelihood plots
*/
struct Printing  {
  float xtick, ytick;
  float StartTick_x, StartTick_y;
  string XLabel, YLabel;
  int labelsize, ticksize;
  float xlabeloffset,ylabeloffset;
  float scaleX,scaleY; // paper scaling to 0..1, 0..1
  
  float TopMargin,BottomMargin,LeftMargin,RightMargin;

  LowLevelPlot::LabelStyle XTickLabelStyle, YTickLabelStyle;

  // cmbeasy gui stuffs
  bool ProtectedSize;
  bool AutomaticTick;
  Printing(int ls=100,int ts=100,float xo=0.08, float yo=0.08, float sx=593.3, float sy=841.9) : labelsize(ls), ticksize(ts), xlabeloffset(xo), ylabeloffset(yo),scaleX(sx), scaleY(sy) , XTickLabelStyle(LowLevelPlot::decimal), YTickLabelStyle(LowLevelPlot::decimal), ProtectedSize(true) , AutomaticTick(true) {
    TopMargin=0.05;
    BottomMargin=0.05;
    LeftMargin = 0.15;
    RightMargin=0.15;
  };
};
 

#endif
