#include "postscriptplot.h"
#include "controlpanel.h"

#include "data.h"
#include "plotwidget.h"
#include "spline.h"

#include <QPaintEvent>
#include <QPainterPath>
#include <QKeyEvent>
#include <QPolygon>
#include <QMessageBox>
#include <QMouseEvent>
#include <QPrinter>
#include <QPrintDialog>
#include <QProcess>
#include <QTime>
#include <QTimer>
#include <QSettings>
#include <QColorDialog>
#include <QDebug>


#include <fstream>
#include <cmath>



PlotWidget::PlotWidget(QWidget * parent, const char * name, Qt::WFlags flags) : QWidget(parent,flags) , paint(0), LogX(false), LogY(false), AntiAlias(true), BadDraw(false), MidMouseDown(false), Pin(QString::fromStdString(ControlPanel::cmbeasyDir("/pix/pin4c.png"))), RepaintReason(NoSpecialReason), DrawBestFitModel(false), Co(2,2000,1e-5,1e-4), xoffset(40), yoffset(20), SplineId(0) , Web(0), Blocks(0), Rectangles(0), data(0) , IsPrinting(false) {
  setObjectName( name );
  colors.resize(100);

  colors[0] = new QColor(0,0,0);
  colors[1] = new QColor(180,0,0);
  colors[2] = new QColor(0,0,180);
  colors[3] = new QColor(0,180,0);
  
  for (unsigned int i=4;i < colors.size(); i++) colors[i]=new QColor(0,0,0);
  
  dataColors.resize(4);
  dataColors[0] = new QColor(100,0,200);
  dataColors[1] = new QColor(20,180,0);
  dataColors[2] = new QColor(180,20,0);
  QFont f( "Helvetica", 10);
  setFont( f );
  setMouseTracking(true);

  QFont g( "Helvetica", 12);

  QPalette palette;
  palette.setColor( backgroundRole(), Qt::white );
  setPalette( palette );

  setFocusPolicy(Qt::ClickFocus);

  LastDraw = new QTime();
  LastDraw->start();
  EnsureBlockDrawing = new QTimer();
  connect(EnsureBlockDrawing,SIGNAL(timeout()),this,SLOT(repaint()));

}

void PlotWidget::paintEvent(QPaintEvent *e) {
  // cout << "Paintevent for : " << this << "RepaintReason: " << RepaintReason << endl;

    checkPainter();
  // paint->fillRect(0,0,width(),Height(),QBrush(Qt::white));
  paint->setPen( Qt::black);
  //cout << "PaintEvent: " << e->rect().left() << " ::  " << e->rect().right() << endl;
  if (BadDraw) {
    paint->setPen( Qt::red);
    paint->drawText(width()/2-100, Height()/2,"Bad Data. Repaint impossible");
  } else {
    try {
      startDrawing(e->rect());
    } catch  (Bad_Error x) {
      QMessageBox::critical(this,"Bad Error", x.s.c_str(),  QMessageBox::Ok,QMessageBox::NoButton);
      cerr << "BAD: ERROR" << endl;
      cerr << x.s << endl;
      BadDraw = true;
    }
  }

  switch (RepaintReason) {
    case LeftButtonPressed:
          RepaintReason = NoSpecialReason;
          paintOnLeftButtonPress();
          break;
    case LeftButtonMoved:
          RepaintReason = NoSpecialReason;
          paintOnLeftMouseMove();
          break;
    case MidButtonMoved:
          RepaintReason = NoSpecialReason;
          paintOnMidMouseMove();
          break;
    case MouseMoved:
          RepaintReason = NoSpecialReason;
          paintOnMouseMove();
          break;
    case NoSpecialReason:
    default:
         break;
  }

  delete paint; paint=0;
}

void PlotWidget::printQt() {
  //printPostscript();
  return;

  QPrinter *prt = new QPrinter();
  QSettings s;
  //prt->setOutputToFile(s.value("/cmbeasy/plot/printqt/tofile",false).toBool());
  prt->setOutputFileName(s.value("/cmbeasy/plot/printqt/name","print.ps").toString());
  prt->setOrientation((QPrinter::Orientation)s.value("/cmbeasy/plot/printqt/orientation",QPrinter::Landscape).toInt());

  QPrintDialog dialog(prt);
  if (dialog.exec()) {
    s.setValue("/cmbeasy/plot/printqt/tofile",prt->outputFileName().isEmpty());
    s.setValue("/cmbeasy/plot/printqt/name",prt->outputFileName());
    s.setValue("/cmbeasy/plot/printqt/orientation",prt->orientation());
    bool aa = AntiAlias;
    AntiAlias = false;
    PrintingHeight = prt->height();
    PrintingWidth = prt->width();

    if (paint) { delete paint; paint=0;}
    paint = new QPainter(prt);
    if (paint == 0) cout << ":::::::::::::::::::::::::::: PAINT is 0 " << endl;
  
    //    cout << "PRINTING HEIGHT " << PrintingHeight <<  "  :: " << PrintingWidth << endl;
    paint->setPen( Qt::black);
    IsPrinting = true;
    startDrawing(QRect(0,0,Width(),Height())); 
    IsPrinting = false;
    delete paint; paint=0;
    AntiAlias = aa;
  } 
  delete prt;
}

void PlotWidget::printPostscript(PostscriptPlot& post, ThingsToPlot& things, const QString fileName)
{
  bool fileNameGiven = !fileName.isEmpty();
  QString outputFile;
  if (fileNameGiven)
   outputFile = fileName;
  else
   outputFile = "/tmp/o.ps";

  ofstream out;

  static QPrinter *prt = new QPrinter();
  QSettings s;
  prt->setOutputFileName(s.value("/cmbeasy/plot/printqt/name",QString()).toString());
  prt->setOrientation((QPrinter::Orientation)s.value("/cmbeasy/plot/printqt/orientation",QPrinter::Landscape).toInt());

  QPrintDialog dialog(prt);
  if (fileNameGiven || dialog.exec()) {
    s.setValue("/cmbeasy/plot/printqt/tofile",prt->outputFileName().isEmpty()); // shouldn't be needed anymore
    s.setValue("/cmbeasy/plot/printqt/name",prt->outputFileName());
    s.setValue("/cmbeasy/plot/printqt/orientation",prt->orientation());

    if (!fileNameGiven && !prt->outputFileName().isEmpty()) out.open(prt->outputFileName().toLatin1().data());
    else out.open(outputFile.toAscii());

    // a4 paper:   x = 595.3    y = 841.9 
    double scaleX=595.3,scaleY=841.9;
    double t;
    switch (prt->pageSize()) {
    case QPrinter::A4: break;
    case QPrinter::A5: t=scaleX; scaleX=scaleY /2; scaleY=t; break;
    case QPrinter::A3: t=scaleX; scaleX=scaleY; scaleY=2*t; break;
    case QPrinter::Letter: scaleX*=216.0/210.0; scaleY *= 279.0/297.0; break;
      //case QPrinter::Custom: scaleX *= post.Width/210; scaleY *= post.Height/297.0;break;
    default: 
      scaleX *= post.Width/210; scaleY *= post.Height/297.0;
    }
    
    post.Width = scaleX;
    post.Height = scaleY;


    // ==============================================
    // == Initialize Confidence Regions
    // ==============================================
  
    if (! RegionGroups.empty()) {
      for (map<int,ConfidenceRegionGroup>::iterator g = RegionGroups.begin(); g != RegionGroups.end(); g++) {
	// first the z-dots, if any
	if (! g->second.cp.v.empty() ) { 
	  ColorPoints cp(g->second.cp);
	  int size = cp.v.size();
	  for (int i = 0; i < size; i++) {
	    pair<double,double> tr = translate(cp.v[i][0],cp.v[i][1]);  // die Punkte nach 0..1 uebersetzen
	    cp.v[i][0] = tr.first;  
	    cp.v[i][1] = tr.second;
	  }
	  things.CPs.push_back(cp);
	  // Also, we need to push back labels for Min and Max values for the color scale
	  TickLabel MinLabel(cp.LeftLower.y);
	  MinLabel.grace =  post.createGraceLabel(cp.Min,LowLevelPlot::decimal,2);
	  MinLabel.size =  post.TickLabelSize*0.01*0.05;
	  MinLabel.offset = -cp.LeftLower.x + 0.01;
	  things.YTickLabels.push_back(MinLabel);

	  TickLabel MaxLabel(cp.LeftLower.y + cp.Height);
	  MaxLabel.grace =  post.createGraceLabel(cp.Max,LowLevelPlot::decimal,2);
	  MaxLabel.size =  post.TickLabelSize*0.01*0.05;
	  MaxLabel.offset = -cp.LeftLower.x + 0.01;
	  things.YTickLabels.push_back(MaxLabel);

	}
	// now the regions
	list<ConfidenceRegion> &Regions = g->second.Regions;	
	for (list<ConfidenceRegion>::iterator i = Regions.begin(); i != Regions.end(); i++) {
	  ConfidenceRegion r =*i; // make a copy
	  for (list<CoordPoint>::iterator k = r.region.begin(); k != r.region.end(); k++) {
	    *k = translateCoordPoint(*k);
	  }
	  things.Regions.push_back(r); 
	}
      }
    }

    // ==============================================
    // == Initialize Special Points (like marker for best fit model)
    // ==============================================
    if ( !things.SpecialPoints.empty() ) {
      list<SpecialPoint>::iterator it = things.SpecialPoints.begin();
      list<SpecialPoint>::iterator end = things.SpecialPoints.end();
      for ( ; it != end; ++it ) {
        it->position = translateCoordPoint( it->position);
      }
    }



    // ==============================================
    // == Initialize Blocks
    // ==============================================
  
    if (Blocks) {
      for (list<Block>::reverse_iterator i = Blocks->rbegin(); i != Blocks->rend(); i++) {
	Block b = *i; // make a copy
	b.lines.clear();
	b.Width /= (Co.X - Co.x);
	for (map<float, list<CoordPoint> >::iterator j = i->lines.begin(); j != i->lines.end(); j++) {
	  float x = translate(j->first,j->first).first;
	  for (list<CoordPoint>::iterator k = j->second.begin(); k != j->second.end(); k++) {
	    float y1 = translate(k->x,k->x).second;
	    float y2 = translate(k->y,k->y).second;
	    b.lines[x].push_back(CoordPoint(y1,y2));
	  }
	}
	things.Blocks.push_back(b); 
      }
    }
    
   
    // ==============================================
    // == Initialize Rectangles
    // ==============================================
  
    if (Rectangles) {
      for (list<Rectangle>::iterator i = Rectangles->begin(); i != Rectangles->end(); i++) {
	  pair<double,double> left_up  = translate(i->C.x,i->C.Y);
	  pair<double,double> left_down  = translate(i->C.x,i->C.y);
	  pair<double,double> right_down  = translate(i->C.X,i->C.y);
	  
	  double x = left_up.first;
	  double Y = left_up.second;
	  double X = right_down.first;
	  double y =  right_down.second ;
	  
	  Rectangle r(x,X,y,Y);	
	  r.FillColor = i->FillColor;
	  things.Rectangles.push_back(r);
	  
      }
    }

    // ============================================
    // ==== Make Line segments out of  the Splines
    // ============================================
    int count = 0;
    for (map<int,Spline*>::iterator i=splines.begin(); i != splines.end();i++,count++) {
      Curve curve;
      double x = Co.x;
      pair<double,double> rel,lastrel;
      x = max(x, i->second->start());
      if (i->second->isWithinBounds(x)) {
	rel = translate(x,i->second->fastY(x));
	lastrel = rel;
	curve.Points.push_back(CoordPoint(rel.first,rel.second));
      } 
      curve.dash = count;
      while (x <= Co.X) {
	if (i->second->isWithinBounds(x)) {
	  rel = translate(x,i->second->fastY(x));
	  float distance = pow(fabs(rel.first - lastrel.first),2) + pow(fabs(rel.second - lastrel.second),2);	  
	  distance = sqrt(distance);
	  if (distance > 0.003 || x > Co.X - (Co.X - Co.x)*3e-3 ) {
	    lastrel = rel; 
	    curve.Points.push_back(CoordPoint(rel.first,rel.second));
	  }	  
	  //	  curve.Points.push_back(CoordPoint(rel.first,rel.second));  // we sample densly, below we thin this out
	}
	x += (Co.X - Co.x)*2e-5;
      }
      /*cout << "BEFORE THIN: " << curve.Points.size() << endl;	cout << "LAST POINT: " << curve.Points.rbegin()->x << " : " <<  curve.Points.rbegin()->y << endl; */


      // curve.thinOut(5e-4);
     

 /*      cout << "AFTER THIN: " << curve.Points.size() << endl;
	      cout << "LAST POINT: " << curve.Points.rbegin()->x << " : " <<  curve.Points.rbegin()->y << endl; */

      curve.color = FromColor(*colors[idCount(i->first)]);  // get the color of this curve and convert from qt-color to my color
      things.Curves.push_back(curve);
    }

    // ============================================
    // === Construct data points
    // ============================================

    if (data) {
      count =0;
      for (list<Data*>::iterator i = data->begin(); i != data->end(); i++,count++) { 
	Data *d = *i;
	for (list<DataEntry>::iterator j = d->points.begin(); j != d->points.end(); j++) { 
	  DataEntry e = *j;
	  // calculate center of data point in 0..1 coordinates 
	  // and get length of error - bar lines in 0..1 coordinates
	  //double x1 = e.x, y1 = e.y;
	  double x2 = e.xRight(), y2 = e.yMiddle();
	  double x3 = e.xLeft(), y3 = e.yMiddle();
	  double x4 = e.xMiddle(), y4 = e.yUpper();
	  double x5 = e.xMiddle(), y5 = e.yLower();
	  //pair<double,double> center = translate(x1,y1);
	  pair<double,double> right_horizontal = translate(x2,y2);
	  pair<double,double> left_horizontal = translate(x3,y3);
	  pair<double,double> up_vertical = translate(x4,y4);
	  pair<double,double> down_vertical  = translate(x5,y5);
	  
	  Curve vertical(up_vertical,down_vertical);
	  Curve horizontal(left_horizontal,right_horizontal);

	  vertical.color = FromColor(*dataColors[count]);
	  horizontal.color = FromColor(*dataColors[count]);
	  things.Curves.push_back(vertical);
	  things.Curves.push_back(horizontal);

	  double TipWidth = 4; // Width of the little bar at the tip of the errorbars
	  pair<double,double> uptip1(up_vertical);
	  pair<double,double> uptip2(up_vertical);
	  uptip1.first -= TipWidth/scaleY;
	  uptip2.first += TipWidth/scaleY;
	  Curve uptip(uptip1,uptip2);
	  uptip.color = horizontal.color;
	  things.Curves.push_back(uptip);
	  
	  pair<double,double> downtip1(down_vertical);
	  pair<double,double> downtip2(down_vertical);
	  downtip1.first -= TipWidth/scaleY;
	  downtip2.first += TipWidth/scaleY;
	  Curve downtip(downtip1,downtip2);
	  downtip.color = horizontal.color;
	  things.Curves.push_back(downtip);

	  
	  pair<double,double> lefttip1(left_horizontal);
	  pair<double,double> lefttip2(left_horizontal);
	  lefttip1.second -= TipWidth/scaleX;
	  lefttip2.second += TipWidth/scaleX;
	  Curve lefttip(lefttip1,lefttip2);
	  lefttip.color = horizontal.color;
	  things.Curves.push_back(lefttip);
	  
	  pair<double,double> righttip1(right_horizontal);
	  pair<double,double> righttip2(right_horizontal);
	  righttip1.second -= TipWidth/scaleX;
	  righttip2.second += TipWidth/scaleX;
	  Curve righttip(righttip1,righttip2);
	  righttip.color = horizontal.color;
	  things.Curves.push_back(righttip);
		  
	  /*
	  cout << left_horizontal.first << " " << down_vertical.second << endl;  // lower left
	  cout << left_horizontal.first << " " <<left_horizontal.second << endl; // left mid
	  cout << right_horizontal.first << " " << right_horizontal.second << endl; // right mid
	  cout << up_vertical.first << " " << up_vertical.second << endl;  // upper mid
	  cout << down_vertical.first << " " << down_vertical.second << endl;  // lower mid 
	  */
	  
	  
	  

	}
      }
    }


    // =============================================
    // the tickmarks  --- first x-Axis
    // =============================================
    /*
    list<double> tick;

    if (logX()) logAxisTicks(Co.x,Co.X,tick,false);   
    else   axisTicks(Co.x,Co.X,tick);

    for (list<double>::iterator i = tick.begin(); i != tick.end(); i++) {
      double x = *i;
      double pos = translate(x,Co.y).first;
      things.XTickMarks.push_back(TickMark(pos));
      //TickLabel label = post.createLabel(x,pos,post.XTickLabelStyle);
      things.XTickLabels.push_back(label);
    }
    */
    for (map<float,string>::iterator i = post.XTickMarkList.begin(); i != post.XTickMarkList.end(); i++) {
      things.XTickMarks.push_back(TickMark(i->first));
      TickLabel label(i->first);
      label.grace = i->second;
      label.size = post.TickLabelSize*0.01*0.05;
      things.XTickLabels.push_back(label);
    }

 

    // =============================================
    // the tickmarks  --- now y-Axis
    // =============================================
    
    /*
    tick.clear();

    if (logY()) logAxisTicks(Co.y,Co.Y,tick,false);   
    else   axisTicks(Co.y,Co.Y,tick);

    for (list<double>::iterator i = tick.begin(); i != tick.end(); i++) {
      double y = *i;
      double pos = translate(Co.x,y).second;
      things.YTickMarks.push_back(TickMark(pos));
      TickLabel label = post.createLabel(y,pos,post.YTickLabelStyle);
      things.YTickLabels.push_back(label);
    } */ 

    for (map<float,string>::iterator i = post.YTickMarkList.begin(); i != post.YTickMarkList.end(); i++) {
      things.YTickMarks.push_back(TickMark(i->first));
      TickLabel label(i->first);
      label.grace = i->second;
      label.size = post.TickLabelSize*0.01*0.05;
      things.YTickLabels.push_back(label);
    }
   
    post.Things = things;
    post.plot(out);  // Go!
  }

  // =================== DONE CREATING THE PS-FILE =================

  if (fileNameGiven || !prt->outputFileName().isEmpty()) return;

  PrinterName = prt->printerName();
  launchLPR();
}

void PlotWidget::startDrawing(QRect rect) {

  // cout << "Stratdrawing in "  << this << endl;

  bool stop=false;
  
  if (Co.x >= Co.X) stop = true;
  if (logX() && ( Co.x <= 0 || Co.X <=0)) stop = true;
  if (logY() && ( Co.y <= 0 || Co.Y <=0)) stop = true; 
  if (Web) {
    if (Co.x < Web->x(0)) stop = true;
    if (Co.y < Web->y(0)) stop = true;
    if (Co.X > Web->x(Web->igrid)) stop = true;
    if (Co.Y > Web->y(Web->jgrid)) stop = true;
  }  
  
  if (stop) { 
    drawAxis();
    paint->setPen( Qt::red);
    paint->drawText((Width()*2)/5,(Height()*2)/5,"x-range invalid (you may consider autoscale)");  
  } else {
    
    bool goon = true;
    if (Blocks) { drawBlocks(); goon=false;}
    if (!RegionGroups.empty()) { drawRegions(); goon=false;}
    
    if (goon) {
      paint->setClipRect(xoffset,0,Width()-xoffset,Height()-yoffset);
      drawData();
      if (Rectangles) drawRectangles(rect);
      // cout << "SPLINES" << endl;
      drawSplines(rect);
      // cout << "done splines" << endl;
    }
    paint->setClipRect(0,0,Width(),Height());
    drawAxis();
  }
  // cout << " finished drawing in "  << this << endl;
}


void PlotWidget::drawAxis() {
  QColor grau(200,200,200);
  paint->setPen(grau);
  drawRect(0.0,1.0,1.0,1.0);
  
  list<double> tick;
  
  QFont f(font());
  QFont big( "Helvetica", 12);
  QFont small("Helvetica",8);

  LowLevelPlot help;

  // cout << "Drawaxis" << endl;
  double ticklength = 2.0/Height();
  double labeloffset = 2.0/Height();

  if (logX()) {
    logAxisTicks(Co.x,Co.X,tick,false);   
    for (list<double>::iterator i = tick.begin(); i != tick.end(); i++) {
      double x =*i; 
      double pos = translate(x,Co.y).first;
      TickLabel label = help.createLabel(x,pos,LowLevelPlot::exponent,1);

      drawLine(pos,-ticklength,pos,ticklength); 
      paint->setPen(Qt::black);
      int expwidth = fontMetrics().width(label.exponent.c_str());
      int basewidth =  fontMetrics().width("10");
      int middle = (expwidth + basewidth)/2;
      drawTextBase(pos,0,label.exponent.c_str(), -middle+basewidth+2,14);
      paint->setFont(big); 
      drawTextBase(pos,0,"10", -middle,19); 
      paint->setFont(f); 
      paint->setPen(grau);
    }
  } else {
    axisTicks(Co.x,Co.X,tick);
    for (list<double>::iterator i = tick.begin(); i != tick.end(); i++) {
      double x =*i; 
      double pos = translate(x,Co.y).first;
      TickLabel label = help.createLabel(x,pos,LowLevelPlot::decimal,3);
      //      cout << "label: " << x << "  " << label.decimal << endl;
      drawLine(pos,0.98,pos,1.02);
      paint->setPen(Qt::black);
      drawTextBase(pos,0,label.decimal.c_str(), -paint->fontMetrics().width(label.decimal.c_str())/2,14);
      paint->setPen(grau);
    }
  }


  tick.clear();
  ticklength = 2.0/Width();
  labeloffset = 2.0/Width();
  if (logY()) {
    logAxisTicks(Co.y,Co.Y,tick,true);
    for (list<double>::iterator i = tick.begin(); i != tick.end(); i++) {
      double y = *i; 
      double pos = translate(Co.x,y).second;
      TickLabel label = help.createLabel(y,pos,LowLevelPlot::exponent,1);
      drawLine(-ticklength,pos,ticklength,pos);
      paint->setPen(Qt::black);
      drawTextBase(-labeloffset,pos,label.exponent.c_str(), -fontMetrics().width(label.exponent.c_str()),14);
      paint->setFont(big); 
      drawTextBase(-labeloffset,pos,"10", -paint->fontMetrics().width(label.exponent.c_str())-  paint->fontMetrics().width("10"),19); 
      paint->setFont(f); 
      paint->setPen(grau);
    }
  } else {
    axisTicks(Co.y,Co.Y,tick);
    for (list<double>::iterator i = tick.begin(); i != tick.end(); i++) {
      double y = *i; 
      double pos = translate(Co.x,y).second;
      TickLabel label = help.createLabel(y,pos,LowLevelPlot::decimal,2);
      drawLine(-ticklength,pos,ticklength,pos); 
      paint->setPen(Qt::black);
      drawTextBase(-labeloffset,pos,label.decimal.c_str(), -paint->fontMetrics().width(label.decimal.c_str()),14);
      paint->setPen(grau);
    }
  }
  // cout << "~drawaxis" << endl;
}


void PlotWidget::autoTick(LowLevelPlot& plot) {
  // =============================================
  // the tickmarks  --- first x-Axis
  // =============================================
  list<double> tick;
  qDebug() << "autoticking" << endl;
  if (logX()) logAxisTicks(Co.x,Co.X,tick,false);   
  else   axisTicks(Co.x,Co.X,tick);
  
  for (list<double>::iterator i = tick.begin(); i != tick.end(); i++) {
    double x = *i;
    double pos = translate(x,Co.y).first;
    LowLevelPlot::LabelStyle style = plot.XTickLabelStyle;
    if (logX()) style = LowLevelPlot::exponent;
    plot.XTickMarkList[pos] = plot.createGraceLabel(x,style,plot.Significant_x);
    cout << "x label for: " << x << "  " <<   plot.XTickMarkList[pos] << endl;
  }


  // =============================================
  // the tickmarks  --- now y-Axis
  // =============================================
    
  tick.clear();
    
  if (logY()) logAxisTicks(Co.y,Co.Y,tick,false);   
  else   axisTicks(Co.y,Co.Y,tick);
    
  for (list<double>::iterator i = tick.begin(); i != tick.end(); i++) {
    double y = *i;
    double pos = translate(Co.x,y).second;
    LowLevelPlot::LabelStyle style = plot.YTickLabelStyle;
    if (logY()) style = LowLevelPlot::exponent;
    plot.YTickMarkList[pos] = plot.createGraceLabel(y, style,plot.Significant_y); 
  }
}

void PlotWidget::semiAutoTick(LowLevelPlot& plot) {
  // =============================================
  // the tickmarks  --- first x-Axis
  // =============================================

 
  int cnt =0;
  for (float x = plot.StartTick_x; x <= Co.X*1.00001; x +=  plot.StepTick_x) {
    double pos = translate(x,Co.y).first;
    if (cnt++ > 20) throw Bad_Error("SemiAuto:: Too many x-axis ticks");
    plot.XTickMarkList[pos] = plot.createGraceLabel(x,plot.XTickLabelStyle,plot.Significant_x);
  }


  // =============================================
  // the tickmarks  --- now y-Axis
  // =============================================
   
  cnt=0;
  for (float y =  plot.StartTick_y; y <= Co.Y*1.00001; y +=  plot.StepTick_y) {
    if (cnt++ > 20) throw Bad_Error("SemiAuto:: Too many y-axis ticks");
    double pos = translate(Co.x,y).second;
    plot.YTickMarkList[pos] = plot.createGraceLabel(y, plot.YTickLabelStyle,plot.Significant_y); 
  }
}

void PlotWidget::logAxisTicks(double min,double max, list<double> &tick,bool yaxis) {

  if (min && max) {
    double lgMin = ceil(log(min)/log(10.0));
    double lgMax = floor(log(max)/log(10.0));
    
    if (lgMax - lgMin > 20) {
//X       cout << "logAxisTicks:  lgMin: " << lgMin << "  lgMax: " << lgMax << endl;
      throw Bad_Error("PlotWidget::logAxisTicks(): Too many ticks");
    }

    for (double lg = lgMin; lg <= lgMax; lg++) {
      tick.push_back(pow(10,lg));
    }
  }
}

double  PlotWidget::axisTicks(double min, double max, list<double> &tick) {
  LowLevelPlot low;
  double range = max - min;
  //qDebug() << "axisticks: min: " << min <<  "  " << max << endl;
  if (range == 0) return 1;

  double sgn=1.0;
  if (min < 0) sgn = -1.0;
  int order = (int)floor(log10(range))-1;
  //  if (order >= 0) order = floor(order); else order = ceil(order);
  
  double start=0, startexp=0,startfactor=0;
  if (min != 0) {
    start = fabs(min);
    startexp = low.splitExponent(start);
    startfactor = low.splitPrefactor(start);
  }

  //cout << "start: " << start << " startexp: " << startexp << "  startfact: " << startfactor << endl;
  //  cout << "order: " << order << endl;
  
  double pre[] = {1,2,2.5,5};
  
  for (int od = order; od <= order+2; od ++) {
    for (int i = 0; i < 4; i++) {
      double skip = pre[i];
      double unit = skip*pow(10.0,od);
      //      cout << "skip: " << skip << " unit: "<<unit << endl;
      double x = (rint(min / unit)-1)*unit;  // start one unit before the start
      tick.clear();
      do {
	if (fabs(x/unit) < 1e-5) x = 0;
	//	cout << "strt: " << x << endl;
	if (x >= min && x <= max) tick.push_back(x);
	x += unit;
      } while (x <= max && tick.size() < 12);
      if (tick.size() < 12) break;
    }
    if (tick.size() < 12) break;
  }
  
  for (list<double>::iterator i = tick.begin(); i != tick.end(); i++) {
    //    cout << "final: " << *i << endl;
  }

  return 1.0;



  




  double factor=0;
  if (range > 0) {
    double order = (int)floor(log(range)/log(10.0));
    if (fabs(order) > 3) {factor = order; order =0;}
    
    double unit = pow(10,factor+order);
    double pack = range/unit;

    //cout << "at first: unit " << unit << endl;
    //cout << "at first: pack " << pack << endl;
    while (pack < 5 && pack*2 <12) { unit*=0.5; pack*=2;}

    int k = (int) rint(min / unit);
  
    
    // cout << "range: " << range << endl;
//     cout << "order:" << order << endl;
//     cout << "faktor: " << factor << endl;
//     cout << "pack: " << pack << endl;
//     cout << "unit: " << unit << endl;
//     cout << "k: " << k << endl;
    int cnt=0;
    while (k * unit < max) {
      if (k*unit > min) {
	if (cnt++ > 20) throw Bad_Error("PlotWidget::axisTicks(): Too many ticks");
	tick.push_back(k*unit * pow(10,-factor));
      }
      k++;
    }
  }  
  return factor;
}

void PlotWidget::drawData() {
  if (data) {
    // cout << "PAINTING DATA" << endl;
    //if (paint ==0) paint = new QPainter(this);
    int count =0;
    for (list<Data*>::iterator i= data->begin(); i != data->end(); i++,count++) drawDataSet(*i,count);
  
    //if (paint) delete paint; paint=0;
  }
}

void PlotWidget::drawDataSet(Data* d, int count ) {
  //  cout << "drawDataSet" << endl;
  paint->setPen(*dataColors[count]);
  //  cout << "PAINTING set " << endl;
  for (list<DataEntry>::iterator i = d->points.begin(); i != d->points.end(); i++) drawSingleDataEntry(*i);
  //  cout << "done" << endl;
}

void PlotWidget::drawSingleDataEntry(const DataEntry &point) {
  //cout << "try to paint point: " << point.x << "  " << point.y   << endl;

  if (point.x()  < Co.x && point.x() > Co.X)  return;

  drawDataLine(point.xMiddle(),point.yLower(), point.xMiddle(), point.yUpper());
  drawDataLine(point.xRight(),point.yMiddle(),point.xLeft(),point.yMiddle());
  drawFixedYDataLine(point.xMiddle(), point.yUpper(), -3, +3);
  drawFixedYDataLine(point.xMiddle(), point.yLower(), -3, +3);
  drawFixedXDataLine(point.xRight(), point.yMiddle(),-3,+3);
  drawFixedXDataLine(point.xLeft(), point.yMiddle(),-3,+3);
}


QString PlotWidget::toStr(double x,int post) {
  QString tmp;
  tmp.setNum(x);
  int k = tmp.indexOf('.');
  int e = tmp.indexOf('e');

  if (k == -1 && e == -1) return tmp;   // not . not e
  if (e == -1) return tmp.left(k+1+post);  //  .  but no e
  
  // both
  
  return tmp.left(k+1+post) + tmp.right(tmp.length()-e);
}




int  PlotWidget::setSpline(const Spline &s) {
  Spline *neu = new Spline(s.size(), "plotwidget"); 
  neu->getData(s);
  neu->arm();
  splines[++SplineId] = neu;
  (*validIds)[SplineId] = true;
  return SplineId;
}

Spline* PlotWidget::getSpline(int idx) {
  Spline *s = 0;
  if (splines.find(idx) != splines.end()) s = splines[idx];
  return s;
}

void PlotWidget::ridOf(int id, bool doRepaint) {
  if (splines.find(id) != splines.end()) splines.erase(splines.find(id));  //erase the spline
  if (validIds->find(id) != validIds->end())  validIds->erase(validIds->find(id));  // erase the id 
  if (doRepaint) repaint();
}

void PlotWidget::checkPainter() {
  if (paint) return;
  paint = new QPainter(this);
}


void PlotWidget::drawSpline(int id, const QColor &c,  Qt::PenStyle  style) {
  QRect rect(0,0,Width(),Height());
  drawSpline(id,c,style,rect);
}

void PlotWidget::drawSpline(int id, const QColor &c,  Qt::PenStyle  style, QRect &rect) {
  //cout << "drawSpline: " << id << endl;
  //if (paint ==0) paint = new QPainter(this);
  
  //cout  << "Plotwidget: " << this << "   drawSpline: " << splines[id] << endl;

  checkPainter();
  paint->setPen(c);
  paint->setClipRect(xoffset,0,Width(),Height()-yoffset);
  //p().setPen(style);
  clearDrawPath();
 
  Spline &s = *splines[id];
  
  int left = max(xoffset, rect.left()-3);
  int right=min(Width(),rect.right()+3);
  int top = max(0,rect.top()-3);
  int bottom = min(Height(),rect.bottom()+3);
  double ymax = reverse_y(top);
  double ymin = reverse_y(bottom);
  double store_x,store_y;
  bool stored=false; // if there is a point that could not be drawn (i.e. y out of bound)
  bool drawn=false; // if at least one point is in the path, true
  //cout << "draw spline in " << left << "  " << right << endl;
  double y;
  for (int pixel=left; pixel <= right; pixel++) { 
    double x = reverse_x(pixel);
    //    cout << " reverse x of " << pixel << "  =  " << x << endl;
    if ( x <=  s.stop()  && x >= s.start())  { 
      y = s(x);
      //      cout << "x: " << x << " y: "  << y << endl;
      if (y <= ymax && y>=ymin) {
	if (stored) {  // one point left to the boundary
	  drawDataPoint(store_x,store_y,true); 
	  stored=false;
	}
	drawDataPoint( x,  y, true );
	drawn=true;
      } else {
	store_x = x; store_y=s(x); stored=true;
	if (drawn) { // so this is one point right to the boundary
	  drawn = false;
	  drawDataPoint(x,s(x),true); 
	}
      }
    } 
  }
  drawPath();
  //if (paint) delete paint; paint=0;
} 


void PlotWidget::drawRectangles(QRect rect) {
  checkPainter();
  paint->setClipRect(xoffset,0,Width(),Height()-yoffset);
  for (list<Rectangle>::iterator i = Rectangles->begin(); i != Rectangles->end(); i++) {
    
    pair<double,double> left_up  = translate(i->C.x,i->C.Y);
    pair<double,double> left_down  = translate(i->C.x,i->C.y);
    pair<double,double> right_down  = translate(i->C.X,i->C.y);
    
    double x = left_up.first;
    double Y = left_up.second;
    double w = right_down.first - x;
    double h =  Y - right_down.second ;

    int ix = tx(x);
    int iY = ty(Y);
    int iw = sx(w);
    int ih = sy(h);
    
    QBrush black(QColor(0,0,0));
    QBrush brush(ToColor(i->FillColor));

    fillRect(ix-1,iY-1,iw+1,ih+1,black);
    fillRect(ix,iY,iw,ih,brush);

  }  
}

void PlotWidget::drawSplines() {
  drawSplines(QRect(0,0,Width(),Height()));
}

void PlotWidget::drawSplines(QRect rect) {
  int count = 0;
  for (map<int,Spline*>::iterator i=splines.begin(); i != splines.end();i++,count++) 
    drawSpline(i->first, *colors[idCount(i->first)] , Qt::SolidLine,rect);
}  

/*!
  Draw Data point, meaning that all data is  translated 
  with respect to minimal and maximal x and y values.
  If store is true then no point is drawn, instead it is added to 
  the DrawPath list 
*/
void PlotWidget::drawDataPoint(double x,double y, bool store) {
  double X,Y;
  if (logX()) X = (log(x)-log(Co.x))/(log(Co.X) - log(Co.x)); else X = (x-Co.x)/(Co.X - Co.x);
  if (logY()) Y = (log(y)-log(Co.y))/(log(Co.Y) - log(Co.y)); else Y = (y-Co.y)/(Co.Y-Co.y);
 
  if (store) DrawPath.push_back(QPoint(tx(X),ty(Y))); else 
    paint->drawPoint(tx(X),ty(Y));
}


/*! translate data coordinates into relative 0 ... 1 coordinates  
  actually, if it exceeds the coordinate boundaries, also values > 1 and < 0 are possible
 */
pair<double,double> PlotWidget::translate(double x,double y) {
  double X,Y;
  //cout << "Y, Co.y : " << Y << "  " << Co.y << "  max " << max(Y,Co.y) << "  ";
  //y = max(y,Co.y);

  if (logX()) X = (log(x)-log(Co.x))/(log(Co.X) - log(Co.x)); else X = (x-Co.x)/(Co.X - Co.x);
  if (logY()) Y = (log(y)-log(Co.y))/(log(Co.Y) - log(Co.y)); else Y = (y-Co.y)/(Co.Y-Co.y);
  return pair<double,double>(X,Y);
}

CoordPoint PlotWidget::translateCoordPoint(CoordPoint& p) {
  pair<double,double> paar = translate(p.x,p.y);
  return CoordPoint(paar.first,paar.second);
}




void PlotWidget::drawDataLine(double x1,double y1, double x2,double y2) {
  // cout << y1 << "  " << Co.y << "  ";
  pair<double,double> a = translate(x1,y1);
  pair<double,double> b = translate(x2,y2);

  int X1= Rint(tx(a.first)), Y1 = Rint(ty(a.second));
  int X2= Rint(tx(b.first)), Y2 = Rint(ty(b.second));
   
  if (max(X2,X1) <= xoffset) return;
  if (min(Y1,Y2) >= Height() - yoffset) return; 

  X1 = max(X1,xoffset);
  X2 = max(X2,xoffset);
  if (min(Y1,Y2) >= Height() - yoffset) return; 
  paint->drawLine(X1,Y1,X2,Y2);
}

void PlotWidget::drawFixedXDataLine(double x,double y,int top,int bottom) { 
  pair<double,double> a = translate(x,y);
  int X1 = Rint(tx(a.first)), Y1 = Rint(ty(a.second)) + top;
  int Y2 = Y1 - top + bottom;
  
  if (X1 <= xoffset) return;
  if (min(Y1,Y2) >= Height() - yoffset) return; 
  Y2 = min(Y2,Height()-yoffset);
  Y1 = min(Y1,Height()-yoffset);
  
  paint->drawLine(X1,Y1, X1,Y2);
}

void PlotWidget::drawFixedYDataLine(double x,double y,int left,int right) { 
  pair<double,double> a = translate(x,y);
  int X1= Rint(tx(a.first)) + left, Y1 = Rint(ty(a.second));
  int X2= X1 - left + right;
  if (X1 <= xoffset) return;
  if (Y1 >= Height() - yoffset) return; 
  paint->drawLine(X1,Y1,X2,Y1);
}

int PlotWidget::tx(double x) { return (int)rint(xoffset +  (Width()-xoffset)*x); }
int PlotWidget::ty(double y) { return (int)rint( (Height()-yoffset)*(1.0-y) ); }



void PlotWidget::altDrawPath() {
  if (DrawPath.size() < 2 ) return;
  
  QPolygon a(100000);
  QPen save = paint->pen();
  int step = 1,point = -1;
  double distance;
  QPoint last;
  for (list<QPoint >::iterator j=DrawPath.begin(); j != DrawPath.end();j++,step++ ) {
    if (step == 1) distance =1e5;    
    if (step != 1) distance = pow((double)last.x()- j->x(),2) + pow((double)last.y() - j->y(),2);
    if (distance > 9) {
      a.setPoint(++point,*j);
      last  = *j;
    }
  }
  // cout << "*********************************" << endl;
  //cout << step << "  and point  " << point << endl;
  a.resize(point);
  paint->setPen(save);
  if (IsPrinting) {
    paint->drawPolyline(a);
  } else paint->drawPoints(a);
}

/*!
  just plots one continous line through all the points in DrawPath 
*/

void PlotWidget::drawPath() {
  if (IsPrinting) {
    altDrawPath(); return;
  }
  if (DrawPath.size() < 2 ) return;

 
  QPolygon neu(DrawPath.size());
  unsigned int k =0;
  for  ( list<QPoint >::iterator i = DrawPath.begin();  i != DrawPath.end(); i++,k++) {
    neu.setPoint(k,*i);
  }
  paint->drawPolyline(neu);
  return;
  
  QPolygon a(100000);
  QPolygon alias(100000), dark(100000);
  
  list<QPoint >::iterator i = DrawPath.begin(); 

  QPen save = paint->pen();
  //int step = -1;
  int point = -1;
  int ap = -1,dc =-1;

  a.setPoint(++point,*i);
  
  list<QPoint >::iterator j;
  for(;;) {
    j = i;
    if (++j == DrawPath.end()) break;
    int x1 = i->x(); 
    int y1= i->y();

    int x2 = j->x();
    int y2 = j->y();

   
    int dx = x2 - x1, dy = y2-y1;

    int total =0;
    if (abs(dx) > abs(dy))  {  // rather horizontal
      int count_y = 0;
      //cout << "horizontal " << endl;
      if (dy) count_y = (int)rint(0.5 * (float)abs(dx)/(float)abs(dy));
      while (total <= abs(dx)) {
	if (a.point(point).x() != x1 ||  a.point(point).y() != y1) {
	  a.setPoint(++point,x1,y1);
	  
	  alias.setPoint(++ap, x1,y1+1);
	  alias.setPoint(++ap, x1,y1-1);
	}
	int yold = y1;
	count_y += abs(dy);
	if (count_y > abs(dx)) {
	  count_y -= abs(dx); 
	  if (dy > 0) y1++; else y1--;
	}
	dark.setPoint(++dc, x1,y1);
	if (dx > 0) x1++; else x1--;
	dark.setPoint(++dc, x1,yold);
	total++;
      }
    } else {
      int count_x = 0;
      if (dx) count_x = (int)rint(0.5 * (float)abs(dy) / (float)abs(dx));
      while (total <= abs(dy)) {
	if (a.point(point).x() != x1 ||  a.point(point).y() != y1) {
	  a.setPoint(++point,x1,y1);
	 
	  alias.setPoint(++ap, x1+1,y1);
	  alias.setPoint(++ap, x1-1,y1);
	}

	count_x += abs(dx);
	int xold = x1;
	if (count_x > abs(dy)) {
	  count_x -= abs(dy); 
	  if (dx > 0) x1++; else  x1--;	    
	}
	dark.setPoint(++dc, x1,y1);
	if (dy > 0) y1++; else y1--;
	dark.setPoint(++dc, xold,y1);
	total++;
      }
    }
    i++;
  }
  //cout << "insgesamt: " << point << endl;
  a.resize(point);

  if (antiAlias()) {

    alias.resize(ap);
    dark.resize(dc);
   
    paint->setPen(blend(0.4));
    paint->drawPoints(alias);
    
    paint->setPen(blend(1.4));
    paint->drawPoints(dark);
  }

  paint->setPen(save);
  if (IsPrinting)  {
    uint j=0;
    QPolygon b(a.size());
      //      cout << a.point(i).x() << "  " << a.point(i).y();

   for ( int i = 0; i < a.size(); i++ ) {
      //      cout << a.point(i).x() << "  " << a.point(i).y();
      if (  ( i&1 ) == 1 ) {
        //cout << "   in b ";
        b.setPoint( j++,a.point( i ).x(),a.point( i ).y() );
      }
      //cout << endl;
    } 
    b.resize(j);    
    paint->drawPolyline(b);
  } else paint->drawPoints(a);
  //  for( list<QPoint >::iterator i = DrawPath.begin();  i != DrawPath.end(); i++) {
  //p().setPen(QColor(255,255,255));
  //p().drawPoint(*i);
  //}
}

QColor PlotWidget::blend(float ratio) {
  const QColor &b = palette().background().color();
  const QColor &c = paint->pen().color();
  float norm = 1.0 + ratio;
  float red = b.red() + ratio* c.red();
  float green = b.green() + ratio * c.green();
  float blue = b.blue() + ratio * c.blue();
  return QColor((int)rint(red/norm), (int)rint(green/norm), (int)rint(blue/norm));
}

/*!
  autoscales the plot such that data and splines nicely fit.
  If notify is true (default) the SIGNAL scaleChanged() is 
  emitted. Set notify to false, if you do not want this 
*/

void PlotWidget::autoScale(bool doX, bool notify)
{
  if (Web) {
    if (doX) {
      Co.x = Web->x(0);
      Co.X = Web->x(Web->igrid);
    }
    Co.y = Web->y(0);
    Co.Y = Web->y(Web->jgrid);
  }


  if (! Blocks && !Web)  {
    Co.Y = -1e100;
    Co.y=1e100;
    if (doX) {
      Co.X = -1e100;
      Co.x=1e100;
    }

    for (map<int,Spline*>::iterator i=splines.begin(); i != splines.end();i++) {
      Spline &s = * i->second;
      //s->printStatus();
      for (int j = 0; j < s.size(); j++) {
        if (doX) {
          if (s.x(j) > Co.X) Co.X = s.x(j);
          if (s.x(j) < Co.x) Co.x = s.x(j);      
        }
        if (s.x(j) >= Co.x && s.x(j) <= Co.X) {  // for honoring x range in doX false case
          if (s.y(j) > Co.Y) Co.Y = s.y(j);
          if (s.y(j) < Co.y) Co.y = s.y(j); 
        }
      }


      for (double x = s.start(); x <= s.stop(); x += s.range(10000)) {
        if (x >= Co.x && x  <= Co.X) {  // for honoring x range in doX false case
          double w = s(x);
          if (w > Co.Y) Co.Y = w;
          if (w < Co.y) Co.y = w;
        }
      }	
    }
//     cout << "AUTOSCALE : Co.Y : " << Co.Y << endl;
//     cout << "AUTOSCALE: Co.y : " << Co.y << endl;
    if (data) {
      for (list<Data*>::iterator i= data->begin(); i != data->end(); i++) {
        for (list<DataEntry>::iterator j =(*i)->points.begin(); j != (*i)->points.end(); j++) {
          if (doX) {
            if (j->xRight() > Co.X) Co.X = j->xRight();
            if (j->xLeft()  < Co.x) Co.x = j->xLeft();
          }
          if (j->x()  >= Co.x && j->x() <= Co.X) {  // for honoring x range in doX false case
            if (j->yUpper() > Co.Y) Co.Y = j->yUpper();
            if (j->yLower() < Co.y && j->yLower() > 0) Co.y = j->yLower();
          }
        }
      }
    }

    Co.Y += (Co.Y  - Co.y)*0.05;
  }
  //  cout << "PLOTWIDGET: " << this << "  has min_y: " << Co.y << "  Max_y: " << Co.Y << endl;
  if (notify) emit(scaleChanged());
}


void PlotWidget::keyPressEvent(QKeyEvent *e) {
  double shift_x = 0.1*(Co.X - Co.x);
  double shift_y = 0.1*(Co.Y - Co.y);
  //cout << "key press even " << endl;
  if (e->key() == Qt::Key_Left) {Co.x -= shift_x; Co.X -=shift_x;}
  if (e->key() == Qt::Key_Right) {Co.x += shift_x; Co.X += shift_x;}
  if (e->key() == Qt::Key_Down) {Co.y -= shift_y; Co.Y -=shift_y;}
  if (e->key() == Qt::Key_Up) {Co.y += shift_y; Co.Y += shift_y;}
  delete paint; paint=0;
  repaint();
  emit (scaleChanged());
}

void PlotWidget::mousePressEvent(QMouseEvent *e) {
  if (e->button() == Qt::MidButton) {
    //cout << "MP: " << e->x() << endl;
    MidMouseDown = true;
    MidMouseMove = false;
    MidMouseX = e->x();
    MidMouseY = e->y();
  }
  if (e->button() == Qt::LeftButton) {
    LeftMouseX = e->x();
    LeftMouseY=e->y();
    RepaintReason = LeftButtonPressed;
    update();
  }

  if (e->button() == Qt::RightButton) {
    if (RegionGroups.empty())
      emit pressedRMB(detectCollision(e->pos()),e->pos());
    else
      emit pressedRMB(detectRegionGroup(e->pos()),e->pos());
  }
}

void PlotWidget::mouseDoubleClickEvent(QMouseEvent *e) {
 if (e->button() == Qt::LeftButton) {
    LeftMouseX = e->x();
    LeftMouseY=e->y();
    emit leftButtonDoubleClicked(QPointF(reverse_x(LeftMouseX), reverse_y(LeftMouseY)));
  }
}

void PlotWidget::colorRegionGroup(int id) {
  QColor A= QColorDialog::getColor();
  if ( !A.isValid() ) return; // dialog was cancelled
  QColor a= QColorDialog::getColor();
  if ( !a.isValid() ) return; // dialog was cancelled
#if 0
  bool ok;
  QColor A= QColorDialog::getRgba(Qt::white, &ok);
  if ( !ok ) return;
  QColor a= QColorDialog::getRgba(Qt::white, &ok);
  if ( !ok ) return;
#endif

  ConfidenceRegionGroup &group  = RegionGroups[id];
  double x = group.Regions.size()-1; // how many need to be shaded ?
  x = max(x,1.0);

  int r = a.red(); int g = a.green(); int b = a.blue(); int t = a.alpha();
  int R = A.red(); int G = A.green(); int B = A.blue(); int T = A.alpha();

  //cout << "r: " << r << "  g: " << g << " b: " << b << endl;
  //cout << "R: " << R << "  G: " << G << " B: " << B << endl;


  double count = 0;
  for (list<ConfidenceRegion>::iterator i = group.Regions.begin(); i != group.Regions.end(); i++) {
    float v_r =  (r + count*(R-r)/x) / 255.;
    float v_g =  (g + count*(G-g)/x) / 255.;
    float v_b =  (b + count*(B-b)/x) /255.;
    float v_t =  (t + count*(T-t)/x) /255.;
    cout << "count: " << count << " v_r: " << v_r   << " v_g: " << v_g  << " v_b: " << v_b << " v_t: "<< v_t << endl;
    count++;
//X     i->FillColor = Color(v_r,v_g,v_b,v_t);
    i->FillColor = Color(v_r,v_g,v_b);
    // check if the pencolor is black or white indicating 1,2,3 sigma boundaries
    float cr = i->PenColor.r;
    float cg =  i->PenColor.g;
    float cb =  i->PenColor.b;
    bool isspecial = false;
    if (cr == 1 && cg == 1 && cb == 1) isspecial = true;
    if (cr == 0 && cg == 0 && cb == 0) isspecial = true;
    if (! isspecial) i->PenColor = Color(v_r,v_g,v_b);
  }
}

void PlotWidget::deleteRegionGroup(int id) {
  map<int,ConfidenceRegionGroup> replace;	       
  int k =0;
  for (unsigned int i =0; i < RegionGroups.size(); i++) {
    if ((int) i != id) replace[k++] = RegionGroups[i];
  }
  RegionGroups = replace;
  drawRegions();
}

void PlotWidget::toFront(int id) {
  cout << "to front for : " << id << endl;
  ConfidenceRegionGroup tmp = RegionGroups[id];
  for ( unsigned int j = (unsigned int)id + 1; j < RegionGroups.size(); j++) RegionGroups[j-1] = RegionGroups[j];
  RegionGroups[RegionGroups.size()-1] = tmp;
  cout << "done to front" << endl;
}
void PlotWidget::toBack(int id) {
  cout << "to back for : " << id << endl;
  ConfidenceRegionGroup tmp = RegionGroups[id];
  for (unsigned int j = RegionGroups.size()-1; j > 0 ; j--) RegionGroups[j] = RegionGroups[j-1];
  RegionGroups[0] = tmp;
}



void PlotWidget::paintOnLeftButtonPress()
{
    checkPainter();
    paint->setPen(Qt::DashLine);
    paint->setPen(QColor(150,150,150));
    paint->drawLine(LeftMouseX,0,LeftMouseX,Height()-yoffset);
    paint->drawLine(xoffset,LeftMouseY,Width(),LeftMouseY);
    paint->setPen(Qt::black);
    displayCoords(reverse_x(LeftMouseX),reverse_y(LeftMouseY));
    paint->setPen(Qt::SolidLine);
}

void PlotWidget::mouseReleaseEvent(QMouseEvent *e) {
  if (e->button() == Qt::MidButton) {
    MidMouseDown = false; 
    if (MidMouseMove) {
      // repaint(OldMidMouse);
      // we got to store this, otherwise reverse_x gets corrupted
      double d2  = pow((double)MidMouseX-e->x(),2) + pow((double)MidMouseY-e->y(),2); 
      // if box is too small, maybe this is only an accident then do
      // nothing
      //cout << "D2: " << d2 << endl;
      if (d2 > 200) {  
	Zoom.push(Co);  // store old
	double mx = reverse_x(OldMidMouse.x()); 
	Co.X = reverse_x(OldMidMouse.x()+OldMidMouse.width());
	Co.x = mx;
	mx = reverse_y(OldMidMouse.y());
	Co.y = reverse_y(OldMidMouse.y()+OldMidMouse.height());
	Co.Y = mx;
      }
    } else  revertZoom();
    delete paint; paint=0;
    repaint();
    emit(scaleChanged());
  }
  if (e->button() == Qt::LeftButton) {
    LeftMouseMove=false;
    delete paint; paint=0;
    repaint();
  }
}

void PlotWidget::leaveEvent(QEvent *e) {
  if ( RegionGroups.empty() ) return;

  unHighlightRegionGroups();
  RepaintReason = MouseMoved;
  repaint();
}

/*!
  Whenever the mouse moves, this handler will be called.
  Depending on wether the left button, the mid button or no 
  button at all is pressed, it will draw a cross hair, an rectangular
  region or try to detect wether the mouse is close to a curve
*/
void PlotWidget::mouseMoveEvent(QMouseEvent *e) {
  //if (e->button() == QMouseEvent::LeftButton) {
  //checkPainter(); 
  //cout << "MM" << e->x() << endl;
  if (MidMouseDown) {
    if (!MidMouseMove) {
      delete paint;paint=0;
      repaint(); // in this case, get rid of the coordinate crosshair
    }
    //cout << "MM and LEFTMOUSE"<<endl;
    MidMouseMove = true;
    int x = MidMouseX, y = MidMouseY;
    int w = x - e->x();
    int h = y - e->y();
    if (w > 0) x = e->x();
    if (h > 0) y = e->y();
    w = abs(w); h = abs(h);
    NewMidMouse = QRect( x, y, w, h );
    RepaintReason = MidButtonMoved;
    update();
  } else if (e->buttons() == Qt::LeftButton) {
    RepaintReason = LeftButtonMoved;
    NewLeftMouseX = e->x();
    NewLeftMouseY = e->y();
    repaint();
  } else {
    mouseMovedPos = e->pos();
    RepaintReason = MouseMoved;
    update();
 }
}

void PlotWidget::paintOnMidMouseMove()
{
  if (OldMidMouse.isValid()) {
 /* Georg: it seems we don't need this?
    QRect &o = OldMidMouse;
    checkPainter(); 
    paint->setPen(Qt::white);
    paint->drawRect(OldMidMouse);
    delete paint; paint=0;
    //cout << "old mid mouse repaint" << endl;
    update(o.x()+1,o.y()+1,2,o.height());
    update(o.right(),o.y(),2,o.height());
    update(o.x(),o.y(),o.width(),2);
    update(o.x(),o.bottom(),o.width(),2);
*/
  }
  checkPainter(); 
  paint->setPen(Qt::black);
  //cout << "drawRect: " << x << "  " << y << "  :: " << abs(w) << "   " << abs(h) << endl;
  OldMidMouse = NewMidMouse;
  paint->drawRect(OldMidMouse);
  //paint->drawRect(x, y, abs(w),abs(h));
}


void PlotWidget::paintOnLeftMouseMove()
{
    // first erase the old stuff :-)
    checkPainter();
/*
    paint->setPen(Qt::white);
    if (LeftMouseX != NewLeftMouseX) paint->drawLine(NewLeftMouseX,0,NewLeftMouseX,Height()-yoffset);
    if (LeftMouseY != NewLeftMouseY) paint->drawLine(xoffset,NewLeftMouseY,Width(),NewLeftMouseY);
*/
    displayCoords(reverse_x(NewLeftMouseX),reverse_y(NewLeftMouseY));
     delete paint;paint=0;
/*
    update(0,NewLeftMouseY,Width(),2);
    update(NewLeftMouseX,0,2,Height());
*/
    // paint->fillRect(0,LeftMouseY,Width(),2,QBrush(QColor(255,0,0)));
    //paint->fillRect(LeftMouseX,0,2,Height(),QBrush(QColor(0,0,255)));
    delete paint;paint=0;
    LeftMouseX = NewLeftMouseX;
    LeftMouseY = NewLeftMouseY;

    checkPainter();
    paint->setPen(Qt::DashLine);
    paint->setPen(QColor(150,150,150));    
    paint->drawLine(LeftMouseX,0,LeftMouseX,Height()-yoffset);
    paint->drawLine(xoffset,LeftMouseY,Width(),LeftMouseY);
    paint->setPen(Qt::black);
    displayCoords(reverse_x(LeftMouseX),reverse_y(LeftMouseY));
    paint->setPen(Qt::SolidLine);
}

void PlotWidget::paintOnMouseMove()
{
  static int clean = 0;
  int id = detectCollision(mouseMovedPos,3);
  if (id) {
    if (id != clean && clean != 0) drawSplines();
    clean = id;
    drawSpline(id, QColor(200,200,0), Qt::SolidLine);
    emit highlighted(id);
  } else {
    if (clean) drawSplines();
    clean = 0;
  }
  // after checking for spline collisions, let's check 
  // for confidenceregions
  id = detectRegionGroup(mouseMovedPos);
  //cout << "id from detectRegion: " << id << endl;
  if (id != -1) {
    if (! RegionGroups[id].Highlighted) { //not yet set for highlighting
      for (map<int,ConfidenceRegionGroup>::iterator g = RegionGroups.begin(); g != RegionGroups.end(); g++) {
        g->second.Highlighted=false;
      }
      RegionGroups[id].Highlighted = true;   
      drawRegions();
    }
  } else {
    // so no region is under cursor, but maybe was ?
    for (map<int,ConfidenceRegionGroup>::iterator g = RegionGroups.begin(); g != RegionGroups.end(); g++) {
      if (g->second.Highlighted) {
        g->second.Highlighted=false;
        drawRegions();
        break;
      }
    }
  } 
}

void PlotWidget::unHighlightRegionGroups()
{
    for (map<int,ConfidenceRegionGroup>::iterator g = RegionGroups.begin(); g != RegionGroups.end(); g++) {
      if (g->second.Highlighted) {
        g->second.Highlighted=false;
        break;
      }
    }
}

int PlotWidget::detectRegionGroup(const QPoint &pos) {
  if ( !visibleRegion().contains( pos ) ) return -1;

  for (map<int,ConfidenceRegionGroup>::reverse_iterator g = RegionGroups.rbegin(); g != RegionGroups.rend(); g++) {
    //    ConfidenceRegionGroup &group = g->second;
    list<QRegion> &QRegions = g->second.QRegions;	
    //    cout << "Size of Regions: " << g->second.Regions.size() << "  QRegions: " << g->second.QRegions.size()  << endl;
    //    cout << endl;
    for (list<QRegion>::iterator i = QRegions.begin(); i != QRegions.end(); i++) {
      if (i->contains(pos)) {
	LastDetectedGroupCollision = CoordPoint(reverse_x(pos.x()),reverse_y(pos.y()));
	return g->first;
      }
    }
  }
  return -1;
}

/*!
  return id of the spline that is closest to the QPoint pos. If there is none, return 0
*/
int  PlotWidget::detectCollision(const QPoint &pos, int accept) {
  // we misuse the drawPath as a list<QPoint> that tell us the position of the curves
  //cout << "detectcollission" << endl;
  double least = 10000;
  int id = 0;
  for (map<int,Spline*>::iterator i=splines.begin(); i != splines.end();i++) {
    clearDrawPath();

    int left = max(xoffset, pos.x() - accept);
    int right = min(Width(), pos.x() + accept);
    left = min(Width(), left);
    right = max(xoffset, right);

    // cout << "left: " << left << "  " << right << endl;
    
    for (int pixel =left; pixel  <=  right; pixel ++) {
      double x = reverse_x(pixel);
      //cout << " reverse x of " << pixel << "  =  " << x << endl;
      if ( x <=  i->second->stop()  && x >= i->second->start() )  drawDataPoint( x,  i->second->fastY(x), true );     
    }
    double min_d = 1000;
    for (list<QPoint>::iterator j = DrawPath.begin(); j != DrawPath.end(); j++) {
      double dx = j->x() - pos.x();
      double dy = j->y() - pos.y();
      double d = sqrt( dx*dx + dy*dy); 
      if (d  < min_d) min_d = d;
    }
    if (min_d < least) { 
      least = min_d; 
      if (least <= accept) id = i->first; 
    }
  }
  //cout << "best match: " << id << "  with " << least << endl;
  return id;
}

double PlotWidget::reverse_x(int x) {
  double frac = ((double)(x - xoffset)) / ((double)(Width() - xoffset)); // frac is between 0 ... 1
  
  if (logX()) return  exp ( log(Co.x) + frac*(log(Co.X) - log(Co.x)) );
  return Co.x + frac*(Co.X - Co.x);
}

double PlotWidget::reverse_y(int y) {
  if (y > Height() - yoffset) return Co.y;
  double frac = ((double)y) / ((double)(Height() - yoffset)); // frac is between 0 ... 1
  frac = 1.0 - frac;

  if (logY()) return  exp ( log(Co.y) + frac*(log(Co.Y) - log(Co.y)) );
  return Co.y + frac*(Co.Y - Co.y);
}

void PlotWidget::printStatus() {
  cout << "PLOTWIDGET:: PRINTSTATUS: " << this << endl;
  cout << "splines.size() " << splines.size() << endl;
  for (map<int,Spline*>::iterator i = splines.begin(); i != splines.end(); i++)
    cout << "content[]: " << i->first << "    points to: "  << i->second << endl;
  cout << "::::: ENDSTATUS" << endl;




}

void PlotWidget::printGnuplot() {
  static QPrinter *prt = new QPrinter();
  cout << "PRINTER AT: " << prt << endl;
  QPrintDialog dlg(prt);
  if (dlg.exec()) {
    bool aa = AntiAlias;
    AntiAlias = false;
    // first, we save all splines in /tmp/cmbeasyprt

    ofstream o("/tmp/cmbeasyprt");
    saveToFile(o);
    o.close();
    
    // next, we generate the gnuplot file

    o.open("/tmp/cmbeasy_gnuplot");
    cout << "plot:print:1" << endl;

    QString output = "/tmp/o.ps";
    QProcess *proc = new QProcess();
    if (!prt->outputFileName().isEmpty()) {
      output = prt->outputFileName();
      // DO SOMETHING ABOUT IT
    }  else {
      PrinterName = prt->printerName();
      connect(proc,SIGNAL(finished(int)),this,SLOT(launchLPR()));
    }
     cout << "plot:print:2" << endl;


     o << "set terminal postscript" << endl;
     o << "set output \"" << output.toStdString() << "\""<<endl;
     if (LogX || LogY) {
       o << "set logscale ";
       if (LogX) o << "x";
       if (LogY) o << "y";
       o << endl;
     }
     o << "plot [" << Co.x << ":" << Co.X << "] ";
     int count =0;
     for (map<int,Spline*>::iterator i=splines.begin(); i != splines.end(); i++) {
       if (count) o << ", ";
       o << "\"/tmp/cmbeasyprt\" ";
       o << "i " << count++;
       o << " w l";
     }
     if (data) {
      for (list<Data*>::iterator i= data->begin(); i != data->end(); i++,count++) {
	if (count) o << ", ";
	o << "\"" << (*i)->FileName << "\" u 1:(" << (*i)->Normal << "*$2):3:(" << (*i)->Normal << "*$4) ";
	o << "title \""<< (*i)->Name << "\"  with xyerrorbars";	
      }
    }

    o << endl;
    cout << "plot:print:3" << endl;

    proc->start("gnuplot", QStringList()<<"/tmp/cmbeasy_gnuplot");
     cout << "plot:print:4" << endl;
     AntiAlias=aa;

  } 
}

void PlotWidget::saveToFile(int i, ostream& o) {
  if (splines.find(i) == splines.end()) return;
  Spline &s = *splines[i];
  for (int k = 0; k < s.size()-1; k++) {
    for (double  X = s.x(k); X <= s.x(k+1); X += 0.05*(s.x(k+1) - s.x(k)) )
      o << X << "\t" << s(X) << endl;
  }
}

void PlotWidget::saveToFile(ostream& o) {
  for (map<int,Spline*>::iterator i=splines.begin(); i != splines.end(); i++) {
    saveToFile(i->first,o);
    o << endl << endl;
  }
}

void PlotWidget::launchLPR() {
  //cout << "LAUNCHLPR: " << PrinterName << endl;
  QProcess *proc = new QProcess();
  QStringList args;
  args << "-P";
  args << PrinterName;
  args << "/tmp/o.ps";
  proc->start("lpr", args);
}

/*
void PlotWidget::drawWeb() { 
  if (!BadDraw) {
    vector< vector<double> >* v;
    v = drawWebLayer(0.95,QColor(200,200,200), true);
    drawWebLayer(0.67,QColor(150,150,150));  
    drawWebLayer(0.1,QColor(100,100,100));  
    delete v; // trust me, this makes sense: the vector created and stored in drawWebLayer is now deleted
    pair<double,double> max = Web->maximum();
    cout << "THE MAXIMUM OF THE WEB IS AT: " << max.first << "  :  " << max.second << endl;
  }
}
*/


void PlotWidget::drawBlocks() {
  //qDebug() << "DRAWBLOCKS" << endl;
  int stepregion = 1; // draw each #stepregion stepregion
  int elapsed = LastDraw->elapsed();      // ask this once
  if (elapsed < 1000) {   // if we redrew recently
    if (! EnsureBlockDrawing->isActive()) { // so if the timer is not already busy
      EnsureBlockDrawing->setSingleShot(true);
      EnsureBlockDrawing->start(1000 - elapsed); // call this
    }
    stepregion = 10; // less regions
  } 
  
  list<Block> &b = *Blocks;
  //  cout << "PlotWidget::drawBLOCKS()"<< endl;
  //  cout << "Coord: " << Co.x << " " << Co.X << "  :: " << Co.y << "  " << Co.Y << endl;
  for (list<Block>:: reverse_iterator i = b.rbegin(); i!=b.rend(); i++) {
    int linewidth = (int)ceil( (Width()-xoffset)*(  i->Width / (Co.X - Co.x) ));
    //for (int i = r.regions-1; i >= 0; i-=stepregion) {
    //    cout << "i: " << i << endl;

    QColor c(ToColor(i->PenColor)); 
    if (i->UsePen) drawBlockArea(i->lines,c,linewidth,true);
    c = ToColor(i->FillColor); 
    drawBlockArea(i->lines,c,linewidth,false);
    
    //    drawBlockArea(i->lines,ToColor(i->FillColor),linewidth,false);
    //    cout << "done draw block area" << endl;
  }
  //  paint->fillRect(20,20,100,100,QBrush(QColor(100,100,100)));
  if (stepregion == 1) LastDraw->start();
}


void PlotWidget::drawBlockArea(map<float, list<CoordPoint> >& lines, QColor& color,int linewidth, bool bound) {
  QBrush brush(color);

  for (map<float, list<CoordPoint > >::iterator i = lines.begin(); i != lines.end(); i++) {
    float x;
    for (list<CoordPoint>::iterator j = i->second.begin(); j != i->second.end(); j++) {
      x = i->first;
      float y1 = j->x;
      float y2 = j->y;

      y1 = (y1 - Co.y)/(Co.Y - Co.y);
      y2 = (y2 - Co.y)/(Co.Y - Co.y);
      x =  (x - Co.x)/(Co.X - Co.x);

      if (bound) {
        fillRect(tx(x)-1,ty(y2)-1,linewidth+2,sy(y2-y1)+2,brush);
      } else {
        fillRect(tx(x),ty(y2),linewidth,sy(y2-y1),brush);
      }
    }
  }
}

void PlotWidget::fillRect(int x,int y,int w,int h, QBrush& brush) {
  
  if (x < xoffset) { x=xoffset; w -= xoffset; }
  if (y < 0) { h += y; y=0; }

  w = min(Width() - x, w);
  h = min(Height() - yoffset - y, h);

  //  h = min(Height()-yoffset - y,h);
  //  if (y+h > Height()-yoffset) h = Height()-y-yoffset;

  if (x > Width() || y > Height()-yoffset) return;
  if (w <0 || h <0) return;

  //  if (x < 0 || y < 0 || w < 0 || h < 0) cout << "AAAH: " << x << "  " << y << "  " << w << "   " << h << endl;


  paint->fillRect(x,y,w,h,brush);
}

void PlotWidget::drawRegions() {
  checkPainter();
  paint->setClipRect(xoffset,0,Width(),Height()-yoffset);
  // For all RegionsGroups
  for (map<int,ConfidenceRegionGroup>::iterator g = RegionGroups.begin(); g != RegionGroups.end(); g++) {
    CoordRect Clip(g->second.Clip);
    CoordPoint ClipLeftUpper(Clip.x,Clip.y);
    CoordPoint ClipRightLower(Clip.X,Clip.Y);

    CoordPoint LeftUpper = translateCoordPoint(ClipLeftUpper);
    CoordPoint RightLower = translateCoordPoint(ClipRightLower);


    int cx = tx(LeftUpper.x);
    int cX = tx(RightLower.x);
    int cy = ty(LeftUpper.y);
    int cY=ty(RightLower.y);

//X     cout << "CLIPRECT: " << xoffset << "  :  " << 0 << "  :  " << Width() << " :   "  << Height()-yoffset << endl;
//X     cout << "CLIPRECT: " << cx << "  :  " << cY << "  :  " << cX << " :   "  << cy << endl;

    cx = max(xoffset,cx);
    cY = max(0,cY);
    cX = min(Width(),cX);
    cy = min(Height()-yoffset,cy);
    cout << "CLIPRECT: " << cx << "  :  " << cY << "  :  " << cX << " :   "  << cy << endl;
    //paint->setClipRect(xoffset,0,Width(),Height()-yoffset);
    paint->setClipRect(cx,cY,cX,cy);  



    list<ConfidenceRegion> & Regions = g->second.Regions;
    g->second.QRegions.clear();
    for (list<ConfidenceRegion>::iterator i = Regions.begin(); i != Regions.end();i++) {
      list<CoordPoint> &current = i->region;
      QPolygon a(current.size());
      int count=0;
      for (list<CoordPoint>::iterator k = current.begin(); k != current.end(); k++,count++) {
	CoordPoint p = translateCoordPoint(*k);
	a[count] = QPoint(tx(p.x),ty(p.y));
	//cout << "( " << k->first << " , " <<  k->second << " )  --> ( " << p.first << ", " << p.second << ")  == ( " << tx(p.first) << ", " << tx(p.second) << ") " << endl;
      }
      Color c = i->FillColor;     
      Color p = i->PenColor;
      //      if (g->second.Highlighted) c=highlightColor(c,1.5);
      //if (g->second.Highlighted) p=highlightColor(p,1.5);
      paint->setPen(ToColor(p));
      QBrush brush(ToColor(c)); 
      if (i->Fill) {
	if (g->second.Highlighted) {
	  paint->setBrush(QColor(255,255,255));
	  paint->drawPolygon(a); 
	  brush.setStyle(Qt::Dense3Pattern);
	  paint->setBrush(brush);
	}
	paint->setBrush(brush);
	paint->drawPolygon(a); 
	if (g->second.Sticky) {
	  CoordPoint pinpos = translateCoordPoint(g->second.PinPos);
	  int x = tx(pinpos.x);
	  int y = ty(pinpos.y);
	  paint->drawPixmap(x,y,Pin);
	}
      } else paint->drawPolyline(a);

      // draw optional colorpoints
      ColorPoints& cp = g->second.cp;
      if (!cp.v.empty()) {
        //	ofstream prove("prove.dat");
        int size = cp.v.size();
        for (int i  = 0; i < size; i++) {
          float x = cp.v[i][0];
          float y = cp.v[i][1];  
          float z = cp.v[i][2];
          Color c;
          if (z > 0.5) PostscriptPlot::InterpolateColor(c,cp.MaxColor,cp.MiddleColor,z,0.5,1);
          else PostscriptPlot::InterpolateColor(c,cp.MiddleColor,cp.MinColor,z,0,0.5);
          QBrush brush(ToColor(c)); 
          paint->setPen(ToColor(c));	  
          paint->setBrush(brush);
          pair<double,double> tr = translate(x,y);

          float X = tr.first,  Y = tr.second;
          float R = cp.Radius;
          paint->drawEllipse(tx(X-R),ty(Y-R),sx(2*R),sy(2*R));

          //	    paint->drawPoint(tx(tr.first),ty(tr.second));	  
          //	  prove << x << " " << y << endl;
        } 
        // Draw A colored scale for the color points
        float dz = 0.01;
        int w = sx(0.04);
        w = max(w,15);
        int p = tx(0.95);
        p = min(p,Width()-20);

        w = 15;
        p = Width()-20;		
        for (float z =  0; z < 1; z +=dz) {
          float y = 0.05 + z*0.3;
          Color c;
          if (z > 0.5) PostscriptPlot::InterpolateColor(c,cp.MaxColor,cp.MiddleColor,z,0.5,1);
          else PostscriptPlot::InterpolateColor(c,cp.MiddleColor,cp.MinColor,z,0,0.5);
          QBrush brush(ToColor(c)); 
          paint->setPen(ToColor(c));	  
          //	  paint->setBrush(brush);
          int h = sy(0.35*dz);
          h = max(2,h);
          paint->fillRect(p,ty(y),w,h, brush);
        }
        paint->setPen(Qt::black);
        paint->setBrush(Qt::NoBrush);
        paint->drawRect(p,ty(0.05+0.3), w,sy(0.3));
        QString upperZValue = QString::number(cp.Max, 'f', 2);
        QString lowerZValue = QString::number(cp.Min, 'f', 2);
        int zWidth = fontMetrics().width(upperZValue);
        paint->drawText(Width()-zWidth-1, ty(0.05+0.3)-0.2*fontMetrics().height(), upperZValue);
        paint->drawText(Width()-zWidth-1, ty(0.0), lowerZValue);
        // end of drawing a colored scale for the color points
      }
      //      cout << "Pushing back QRegion: " << &current << endl;
      g->second.QRegions.push_back(a);
      // store as QRegion for later detection
      }
    }

  if ( DrawBestFitModel ) {
    float x = BestFitModel.first;
    float y = BestFitModel.second;
    QBrush brush(Qt::red); 
    paint->setPen(Qt::red);	  
    paint->setBrush(brush);

    pair<double,double> tr = translate(x,y);

    float X = tr.first,  Y = tr.second;

    QPainterPath starPath;
    float R = 12e-3;
    starPath.moveTo( tx( X+R ), ty( Y ) );
    for ( int i = 1; i < 5; ++i ) {
      starPath.lineTo(  tx( X ) + sx( R ) * cos( 0.8 *i * M_PI ) ,
                        ty( Y ) + sy( R ) * sin( .8 *i * M_PI ) );
    }
    starPath.closeSubpath();
    paint->drawPath( starPath );
  }

  if ( DrawModelToMark ) {
    float x = ModelToMark.first;
    float y = ModelToMark.second;
    QBrush brush(Qt::blue); 
    paint->setPen(Qt::blue);	  
    paint->setBrush(brush);

    pair<double,double> tr = translate(x,y);
    float X = tr.first,  Y = tr.second;
    float R = 12e-3;
    paint->drawEllipse(tx(X-R),ty(Y-R),sx(2*R),sy(2*R));
  }

  delete paint; paint=0;
  checkPainter();
}

int PlotWidget::Height() {
  if (IsPrinting) {
    return PrintingHeight;
  } else return height(); 
}

int PlotWidget::Width() {
  if (IsPrinting) {
    return PrintingWidth;
  } else return width(); 
}

/*!
  If store is true, then Web is rasterize. Thus, *always* call the drawWebLayer with store = true,
  if it is the first layer to draw

vector< vector<double> >*   PlotWidget::drawWebLayer(double frac, const QColor& color,bool store) {
 double total;
 static  vector< vector<double> >* v;
 if (store) v = Web->rasterize(Width()-xoffset,Height()-yoffset,total,Co.x,Co.X,Co.y,Co.Y);

 // the following line may throw an exception, we catch this in repaint()
 list < pair<int,int> > *lst   = Web->fractionOfRaster(v,frac);
 
 QPen pen(color);
 QBrush brush(color);

 paint->setPen(pen);
 
 for (list< pair<int,int> >::iterator i = lst->begin(); i != lst->end(); i++) {
   int x1=i->first+ xoffset , y1 = Height()-yoffset- i->second;
   i++;
   int x2=i->first+xoffset, y2=Height()-yoffset - i->second;
   if (Printing) {
     paint->fillRect(x1,y2,1,abs(y2-y1),brush);
   }
   else paint->drawLine(x1,y1,x2,y2);
 }
 delete lst;
 return v;
}
*/

 /*! 
    make a 2d-Splineweb the player of the game :-). 
    With this, you give control over to the plotwidget. Do not delete or modify the web yourself !!!  
    Do not specify an anchor for this web, i.e. create without anchor.
  
void PlotWidget::setWeb(SplineWeb* w) {
  if (Web) delete Web;
  Web =w;
  BadDraw = false;
  autoScale();
}
 */

void PlotWidget::setBlocks(list<Block>* b) {
  delete Blocks;
  Blocks =b;
}

/*
void PlotWidget::setRegions(list<ConfidenceRegion>* r) {
  if (Regions) delete Regions;
  Regions =r;
}
*/

void PlotWidget::addRegionGroup(list<ConfidenceRegion>* r) {
  map<int,ConfidenceRegionGroup> replace;	       
  int k =0;
  for (unsigned int i =0; i < RegionGroups.size(); i++) {
    if (RegionGroups[i].Sticky) replace[k++] = RegionGroups[i];
  }
  RegionGroups = replace;
  int size = RegionGroups.size();
  RegionGroups[size] = ConfidenceRegionGroup(*r);
  RegionGroups[size].Clip = Co;
}
void  PlotWidget::attachColorPoints(ColorPoints& cp) {
  int size = RegionGroups.size();
  RegionGroups[size-1].cp = cp;
} 


void PlotWidget::setRegionGroupSticky(int id) {
  RegionGroups[id].Sticky = true;
  RegionGroups[id].PinPos = LastDetectedGroupCollision;
  drawRegions();
  repaint();
}



void PlotWidget::setRectangles(list<Rectangle>* r) {
  if (Rectangles) delete Rectangles;
  Rectangles =r;
}

Color PlotWidget::highlightColor(Color &c,float factor) {
  c.r += 0.5; c.b += 0.5; c.g += 0.5;
  if (c.r > 1.0) c.r =1.0;
  if (c.g > 1.0) c.g =1.0;
  if (c.b > 1.0) c.b =1.0;
  return c;
}


int PlotWidget::idCount(int id) {
  int count =-1;
  for (map<int,bool>::iterator i = validIds->begin(); i != validIds->end(); i++) {
    count++;
    if (i->first == id) return count;
  }
  return 0;
}

void PlotWidget::revertZoom() {
  if (Zoom.empty()) return;
  Co = Zoom.top();
  Zoom.pop();
}


void PlotWidget::displayCoords(double x, double y) {
  QString t = "["+toStr(x,4) + " : " +toStr(y,4)+" ]";
  paint->drawText(width()-fontMetrics().width(t), 14, t);
}

