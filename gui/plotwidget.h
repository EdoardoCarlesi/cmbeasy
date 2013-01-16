#ifndef PLOTWIDGET_H
#define PLOTWIDGET_H

#include "cleanvector.h"
#include "analyzethis.h"

#include <QWidget>
#include <QPainter>
#include <QMouseEvent>
#include <QLabel>
#include <QKeyEvent>
#include <QPaintEvent>
#include <QColor>
#include <QPixmap>
#include <QRect>
#include <QPoint>

#include <math.h>
#include <list>
#include <vector>
#include <stack>
#include <map>
#include <string>
#include <fstream>


class QPaintEvent;
class QPainter;
class Spline;
class Data;
class DataEntry;
class QMouseEvent;
class QPoint;
class QLabel;
class SplineWeb;
class QKeyEvent;
class QTime;
class QTimer;
class Printing;

#include "postscriptplot.h"


/*!
  A group to hold several confidence regions.
  Say the regions for 0.1 .... 2 sigma for drawing
  the confidence limits on a given model
*/
struct ConfidenceRegionGroup {
  ConfidenceRegionGroup(  list<ConfidenceRegion> &r) : PinPos(0,0) ,Clip(0,0,1,1), Highlighted(false) , Sticky(false) ,Regions(r)   {}
  ConfidenceRegionGroup() : PinPos(0,0), Clip(0,0,1,1), Highlighted(false), Sticky(false) {} 
  CoordPoint PinPos;
  CoordRect Clip;
  bool Highlighted;
  bool Sticky;
  list<ConfidenceRegion>  Regions;
  list<QRegion> QRegions; //!< QRegions that belong to the Regions that's for cursor over region detection
  ColorPoints cp; //!< optional color points for z-information
};

class PlotWidget : public QWidget {
  Q_OBJECT

  QPainter *paint;

  bool LogX,LogY;
  QString PrinterName;

  bool AntiAlias;
  bool BadDraw;
  bool MidMouseDown, MidMouseMove,LeftMouseMove;
  int MidMouseX, MidMouseY;
  int LeftMouseX,LeftMouseY, NewLeftMouseX, NewLeftMouseY;
  QPoint mouseMovedPos;
  QRect OldMidMouse, NewMidMouse;
  QTime *LastDraw;
  QTimer *EnsureBlockDrawing;
  enum RepaintReasonType {
       LeftButtonPressed, LeftButtonMoved,
       MidButtonMoved, MouseMoved, NoSpecialReason
  };

  CoordPoint LastDetectedGroupCollision;
  QPixmap Pin;
  RepaintReasonType RepaintReason;
  bool DrawBestFitModel;
  bool DrawModelToMark;


 protected:
  virtual void mousePressEvent ( QMouseEvent * e );
  virtual void mouseMoveEvent ( QMouseEvent * e );
  virtual void mouseReleaseEvent( QMouseEvent * e );
  virtual void mouseDoubleClickEvent ( QMouseEvent * e );
  virtual void leaveEvent(QEvent *e);

 public:

  PlotWidget(QWidget * parent=0, const char * name=0, Qt::WFlags f=0 );
  ~PlotWidget() {};


  pair<float,float> BestFitModel;
  pair<float,float> ModelToMark;
  CoordRect Co;
  int xoffset,yoffset; //! holding the number of points that are offsets from the left and bottom border of the widget until the graph - accesible area begins. In other words: this leaves space for axis and labels 


  //QPainter &p() { return *paint; } 

  int SplineId;
  map<int,Spline*> splines;
  SplineWeb* Web; //!< if this plot is for contourplotting, you give a web
  //  RasterizeReturn *Blocks;
  list<Block> *Blocks;
  //  list<ConfidenceRegion>* Regions;
  map<int,ConfidenceRegionGroup> RegionGroups;
  CleanVector <QColor*> colors;
  CleanVector <QColor*> dataColors;
  list<QPoint > DrawPath;
  stack<CoordRect> Zoom;
  list<Rectangle> *Rectangles; 
  list<Data*> *data;  //!< pointer to data list (points to cmbmainwindow.h)

  void setDataList(list<Data*> * d) {data = d;}  //!< set the data pointer
  void paintEvent ( QPaintEvent * ); 
  void paintOnLeftButtonPress();
  void paintOnLeftMouseMove();
  void paintOnMidMouseMove();
  void paintOnMouseMove();

  QColor blend(float);
  void checkPainter();

  void incSplineId() { SplineId++;} //!< add one (i.e. skip one id number);

  bool IsPrinting; 
  int PrintingHeight, PrintingWidth; //!< used by Height() and set by printQt()
  int Height(); //!< usually height(), except for printing
  int Width(); //!< see Height()

  void drawAxis();
  double axisTicks(double,double,list<double>&); //!< fill map with sensible ticks from min to max
  void logAxisTicks(double,double,list<double>&,bool yaxis); //!< fill map with sensible ticks from min to max this is the log version and it will return the exponent with base 10 as a string belonging to the relative 0...1 coordinate which you may use for all functions taking these (like e.g. drawTextBase() etc).
  void autoTick(LowLevelPlot&);
  void semiAutoTick(LowLevelPlot&);

#if 0
  //  QColor ToColor(vector<float> c) { return QColor((int)(255*c[0]), (int)(255*c[1]),(int)(255*c[2])); } //!< convert 3 dim vector c to qcolor. vector is r-g-b and range 0..1
  QColor ToColor(Color& c) { return QColor((int)(255*c.r), (int)(255*c.g),(int)(255*c.b),(int)(255*c.alpha)); } //!< convert Color for printing to Qt Color class
  Color FromColor(QColor& c) { return Color(c.red()/255.0, c.green()/255.0,c.blue()/255.0, c.alpha()/255.);} //! convert Qt Color to printing Color class (from lowlevelplot.h)
#endif

  //  QColor ToColor(vector<float> c) { return QColor((int)(255*c[0]), (int)(255*c[1]),(int)(255*c[2])); } //!< convert 3 dim vector c to qcolor. vector is r-g-b and range 0..1
  QColor ToColor(Color& c) { return QColor((int)(255*c.r), (int)(255*c.g),(int)(255*c.b)); } //!< convert Color for printing to Qt Color class
  Color FromColor(QColor& c) { return Color(c.red()/255.0, c.green()/255.0,c.blue()/255.0);} //! convert Qt Color to printing Color class (from lowlevelplot.h)
  Color highlightColor(Color& c,float factor); //!< highlight color. Not used currently

  void drawRectangles(QRect);
  void drawSplines(QRect); //!< draw all splines of plot within recange
  void drawSplines(); //!< draw all splines of plot
  void drawSpline(int id, const QColor&, Qt::PenStyle, QRect&);
  void drawSpline(int id, const QColor&, Qt::PenStyle);
  void drawData(); //!< draw all data sets 
  void drawDataSet(Data* , int); //!< Draw a whole data set (e.g. boomerang)
  void drawSingleDataEntry(const DataEntry &);  //!< Draw one Entry of a Data-Set
  //  void drawWeb();
  
  void drawBlockArea(map<float, list<CoordPoint> >& lines, QColor& color,int linewidth,bool bound);
  vector< vector<double> >* drawWebLayer(double frac, const QColor&, bool store=false); 

  void drawRegions();

  void autoScale(bool doX = true,bool notify=true);
  void setLog(bool x,bool y) { LogX = x; LogY = y;}//!< set logarithmic scale for x and y axis
  bool logX() { return LogX;}
  bool logY() { return LogY;}

  int tx(double x); //!< transform x- range 0...1 to graph-world 
  int ty(double y); //!< transform y-range 0...1  to graph-world  

  pair<double,double> translate(double x,double y); //!< traxnsform data (x,y) coordinate to relative 0 .. 1 range (which in turn you may scale using tx(), ty() to the widget world
  CoordPoint translateCoordPoint(CoordPoint& p); //!< convenience wrapper for translate using coordpoint


  double reverse_x(int x); //!< transfrom from pixel back to physical x
  double reverse_y(int y); //!< transfrom from pixel back to physical y

  int sx(double x) { return (int)rint((Width()-xoffset)*x);}
  int sy(double y) {return (int)rint((Height()-yoffset)*y);}
  void drawRect(double x,double y,double w,double h) {
    //cout << "drawRect: " << (tx(x)) << "  "  <<  ty(y) << "  " << sx(w) << "  " << sy(h) << endl;
    paint->drawRect(tx(x),ty(y),sx(w),sy(h));
    //cout << "finished " << endl;
  }
  void fillRect(int x, int y, int w, int h, QBrush&); //!< fill with brush but clip at axis !

  void keyPressEvent(QKeyEvent*);

  void drawText(double x,double y,const QString &s) {
    paint->drawText(tx(x),ty(y),s);
  }

  void drawTextBase(double x,double y,const QString &s,int XOffset=0,int YOffset=0) {
    //int X = (int)rint(tx(x)) + XOffset;
    //int Y = (int)rint(ty(y)) + YOffset;
    //cout << "DTB: " << X << "  " << Y << "  :: " << s << endl;
    paint->drawText(tx(x) + XOffset,ty(y) + YOffset,s);
  }

  void drawLine(double x1,double y1,double x2,double y2) {
    paint->drawLine(tx(x1),ty(y1),tx(x2),ty(y2));
  }

  int Rint(double x) { return (int)rint(x);}

  void drawDataPoint(double x,double y, bool store=false);
  void drawDataLine(double x1,double y1,double x2,double y2); //!< Draw a line from Data point (x1,y1) to data point (x2,y2) 
  void drawFixedXDataLine(double x,double y,int,int);
  void drawFixedYDataLine(double x,double y,int,int);
  void altDrawPath();
  void drawPath();
  void startDrawing(QRect);

  QString toStr(double,int post=0);

  int setSpline(const Spline& s); 
  Spline *getSpline(int);
  void ridOf(int id, bool doRepaint=true); //!< kill plot with that id 
  void clearDrawPath() { DrawPath.clear();}
  int detectCollision(const QPoint &pos,int=5);

  int detectRegionGroup(const QPoint &pos);
  void unHighlightRegionGroups();



  //  void setWeb(SplineWeb* w);  //!< make a 2d-Splineweb the player of the game :-). 

  void setBlocks(list<Block> *b); 
  void setRectangles(list<Rectangle> *);
  void addRegionGroup(list<ConfidenceRegion>* r);
  void clearRegionGroups() { RegionGroups.clear(); }
  void deleteRegionGroup(int id);
  void colorRegionGroup(int id);
  void toFront(int id); //!< bring region to front
  void toBack(int id); //!< send region back
  void setRegionGroupSticky(int id); //!< make region sticky
  void attachColorPoints(ColorPoints&); //!< attach ColorPoints to the last RegionGroup


  void printStatus();
  void printGnuplot();  //!< print using gnuplot 
  void saveToFile(ostream&o); //!< save all splines in one file, separated by double blank line
  void saveToFile(int i,ostream&o); //!< save spline with index int to ostream

  bool antiAlias() { return AntiAlias; }


  map<int,bool> *validIds;  //!< keep track of ids that are in the game 
  int idCount(int id); //!< given the  id, get the position if all id's are ordered (they are ordered :-), i.e. lowest id means 0 etc

  void displayCoords(double,double); //!< draw coordinates e.g. [12.34 ,1e-4]


 public slots:
   void setDrawBestFitModel(bool b) { DrawBestFitModel = b; }
   void setDrawModelToMark(bool b) { DrawModelToMark = b; }
   void launchLPR(); //!< do lpr /tmp/o.ps
   void setAntiAlias(bool b) { AntiAlias = b; repaint();}
   void printQt(); //!< print using Qt painter
   void printPostscript(PostscriptPlot&, ThingsToPlot&, const QString fileName = QString());
   void revertZoom(); //!< get stored zoom in
   void drawBlocks();

 signals:
  void pressedRMB(int, const QPoint&); //!< sends id (or 0) for RMB pressed event
  void leftButtonDoubleClicked(const QPointF&); //!< sends _physical_ coordinates on left mouse button press
  void highlighted(int);
  void scaleChanged(); //!< whenever autoscale is called 
};


#endif
