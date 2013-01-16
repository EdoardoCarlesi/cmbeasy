#include "postscriptplot.h"
#include "global.h"
#include <iostream>
#include <math.h>

//PostscriptPlot::PostscriptPlot(string resourcefile) :  LowLevelPlot(), Resource(resourcefile) {}
PostscriptPlot::PostscriptPlot(const LowLevelPlot& l) : LowLevelPlot(l) {}

void PostscriptPlot::plot(ofstream &out) {
  prepare();
  ifstream commands(Resource.c_str());

  if (commands.bad()) throw Bad_Error("No command file found");
  
  out << "%!PS-Adobe-3.0 EPSF-3.0" << endl; 
  out << "%%BoundingBox: 0 0 " << Width << "  " << Height << endl;
  out << "%%Creator: cmbeasy" << endl;
  
  out << "/scalex{" << 1.0/Height << " mul } def" << endl;
  out << "/scaley{"<< 1.0/Width << " mul } def" << endl;
  out << "/scalexy{" << sqrt( 1.0/(Height*Height) + 1.0/(Width*Width) ) << " mul } def" << endl;

  char ch;
  while (commands.get(ch)) out.put(ch);
  
  
  // =============================================
  //  == Translate and scale coordinate system
  // =============================================
  

  //
  // go to landscape and (0..1) coordinates
  //
  
  out << endl;
  out << Width << " " << Height << " scale" << endl;
  out << "1 0 translate" << endl;
  out << "90 rotate" << endl;

  //
  // move and scale further to viewport
  // 

  //out << XShiftLow << " " << YShiftLow << " translate" << endl;
  out << LeftMargin << "  " << BottomMargin <<  " translate" << endl;
  out << ScaleX2View << "  " << ScaleY2View << " scale"<<endl;
  

  //
  // scale standard line-width
  // 
  // 
  out << StandardLineWidth * CurveLineWidth * 0.01 <<  " setlinewidth" << endl;

  //
  //  make clipping part for everything that is inside the plot
  //
  out << "gsave" << endl;
  out << "newpath" << endl;
  out << "0 0 moveto" << endl;
  out << "1 0 lineto" << endl;
  out << "1 1 lineto" << endl;
  out << "0 1 lineto" << endl;
  out << "closepath clip" << endl;
 
  // ==============================================
  // == PRINT Rectangles
  // ==============================================
  
  if (! Things.Rectangles.empty()) {
    for (list<Rectangle>::iterator i = Things.Rectangles.begin(); i != Things.Rectangles.end(); i++) {
      
      double x = i->C.x, X=i->C.X;
      double y = i->C.y, Y=i->C.Y;

      out << "newpath" << endl;
      out << x << " " << Y << " " << " moveto" << endl;
      out << X << " " << Y << " " << " lineto" << endl;
      out << X << " " << y << " " <<  " lineto" << endl;
      out << x << " " << y << " " <<  " lineto" << endl;
      out <<"closepath" << endl;
      out << i->FillColor.r << " " <<  i->FillColor.g << " " << i->FillColor.b << " " << " setrgbcolor" << endl;
      out <<"fill" << endl;

      out << "newpath" << endl;
      out << x << " " << Y << " " << " moveto" << endl;
      out << X << " " << Y << " " << " lineto" << endl;
      out << X << " " << y << " " <<  " lineto" << endl;
      out << x << " " << y << " " <<  " lineto" << endl;
      out <<"closepath" << endl;
      out << "0 setgray"  << endl;
      out <<"stroke" << endl;
    }
  }
   // ==============================================
  // == PRINT Regions
  // ==============================================

  if (! Things.Regions.empty()) {
    for (list<ConfidenceRegion>::iterator i = Things.Regions.begin(); i != Things.Regions.end(); i++) {
      out << "% %%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      out << "% Region: " <<     i->PenColor.r << " " <<  i->PenColor.g << " " << i->PenColor.b << " " << " as color" << endl;
      out << "% %%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      out << "newpath" << endl;
      //      string what = " moveto";
      float last_x=0, last_y=0;
      for (list<CoordPoint>::iterator k = i->region.begin(); k!= i->region.end(); k++) {
	if (k == i->region.begin()) {
	  out <<  k->x << " " << k->y << " moveto" << endl;
	  last_x = k->x; last_y = k->y;
	} else {
	  out <<  mac(k->x-last_x) << " " << mac(k->y-last_y) << " r"<<endl;
	  last_x = k->x; last_y = k->y;
	  //	  out <<  mac(k->x) << " " << mac(k->y) << " r"<<endl;
	}
      }
      out << "closepath" << endl;
      if (i->Fill) {
	out << "gsave" << endl;
	out << i->FillColor.r << " " <<  i->FillColor.g << " " << i->FillColor.b << " " << " setrgbcolor" << endl;
	out <<"fill grestore" << endl;
      }
      out << "gsave" << endl; 
      out << i->PenColor.r << " " <<  i->PenColor.g << " " << i->PenColor.b << " " << " setrgbcolor" << endl;
      out << StandardLineWidth * CurveLineWidth * 0.01 <<  " setlinewidth" << endl;
      out << "stroke grestore" << endl;
    }
  }

  // ==============================================
  // == PRINT Blocks
  // ==============================================
  
  if (! Things.Blocks.empty()) {
    for (list<Block>::iterator i = Things.Blocks.begin(); i != Things.Blocks.end(); i++) {
      out << i->FillColor.r << " " <<  i->FillColor.g << " " << i->FillColor.b << " " << " setrgbcolor" << endl;
      float mc = 100000; // scaling factor for mac [saves file space]
      for (map<float, list<CoordPoint> >::iterator l = i->lines.begin(); l !=i->lines.end(); l++) {
	float x = l->first;
	for (list<CoordPoint>::iterator j = l->second.begin(); j != l->second.end(); j++) {
	  float y1 = j->x;
	  float y2 = j->y;
	  out << mac(i->Width,mc) << " " << mac(y2 - y1,mc) << " " <<  mac(y1,mc) << " " << mac(x,mc)<< " b" << endl;
	}
      }
    }
  }

  // ==============================================
  // == PRINT ColorPoints
  // ==============================================
  
  if (! Things.CPs.empty()) {
    for (list<ColorPoints>::iterator i = Things.CPs.begin(); i != Things.CPs.end(); i++) {
      //  out << i->FillColor.r << " " <<  i->FillColor.g << " " << i->FillColor.b << " " << " setrgbcolor" << endl;
      ColorPoints &cp = *i;
      int size = cp.v.size();
      for (int i  = 0; i < size; i++) {
	float x = cp.v[i][0];
	float y = cp.v[i][1];  
	float z = cp.v[i][2];
	Color c;
	if (z > 0.5) InterpolateColor(c,cp.MaxColor,cp.MiddleColor,z,0.5,1);
	else InterpolateColor(c,cp.MiddleColor,cp.MinColor,z,0,0.5);
	out << x << " " << y << " " << cp.Radius << "  " << c.r << " " << c.g << " " << c.b << " dc" << endl;
      }
      // draw color scale
      Color c;
      float dz = 0.01;
      float y,x =  cp.LeftLower.x;
      float Y,X = x + cp.Width;
      for (float z =  0; z < 0.99; z +=dz) {
	y = cp.LeftLower.y + z*cp.Height;
	if (z > 0.5) InterpolateColor(c,cp.MaxColor,cp.MiddleColor,z,0.5,1);
	else InterpolateColor(c,cp.MiddleColor,cp.MinColor,z,0,0.5);
	
	Y = y + dz*cp.Height;
	out << "newpath" << endl;
	out << x << " " << Y << " " << " moveto" << endl;
	out << X << " " << Y << " " << " lineto" << endl;
	out << X << " " << y << " " <<  " lineto" << endl;
	out << x << " " << y << " " <<  " lineto" << endl;
	out <<"closepath" << endl;
	out << c.r << " " <<  c.g << " " << c.b << " " << " setrgbcolor" << endl;
	out <<"fill" << endl;
      }   
      // finally, put a black box around the color scale
      y =  cp.LeftLower.y;
      Y = y + cp.Height;
      out << "newpath" << endl;
      out << x << " " << Y << " " << " moveto" << endl;
      out << X << " " << Y << " " << " lineto" << endl;
      out << X << " " << y << " " <<  " lineto" << endl;
      out << x << " " << y << " " <<  " lineto" << endl;
      out <<"closepath" << endl;
      out << " 0 0 0 setrgbcolor" << endl;
      out <<"stroke" << endl;
    }
  }


  // ==============================================
  // == PRINT SpecialPoints (like marker for best fit model)
  // ==============================================

  for (list<SpecialPoint>::iterator i = Things.SpecialPoints.begin(); i != Things.SpecialPoints.end(); ++i) {
    SpecialPoint& p = *i;
    float x = p.position.x;
    float y = p.position.y;
    Color c;
    float R = p.size;
    switch (p.type) {
      case SpecialPoint::Circle:
        c.r=0;c.g=10;c.b=100;
	out << x << " " << y << " " << R << "  " << c.r << " " << c.g << " " << c.b << " dc" << endl;
        break;
      case SpecialPoint::Star:
        c.r=100;c.g=0;c.b=0;
        out << "newpath" << endl;
        out << (x+R) << " " << y << " " << " moveto" << endl;
        for ( int i = 1; i < 5; ++i ) {
          out << ( x +  R  * cos( 0.8 *i * M_PI ) ) << " " <<
            ( y +  R * sin( .8 *i * M_PI ) ) << " lineto" << endl;
        }
        out <<"closepath" << endl;
        out << c.r << " " <<  c.g << " " << c.b << " " << " setrgbcolor" << endl;
        out <<"fill" << endl;
        break;
      default: //shouldn't happen
        throw Bad_Error("PostscriptPlot::plot() - unknown SpecialPoint::type to plot");
    }
  }

  // ==============================================
  // == PRINT Curves
  // ==============================================
  
  if (! Things.Curves.empty()) {
    for (list <Curve>::iterator i = Things.Curves.begin(); i != Things.Curves.end(); i++) {
      out << "newpath" << endl;
      out << "dash" << i->dash % 4 << endl; // 5  dashes  (0..4) in total, so modulo 4 
      out << i->color.r << " " <<  i->color.g << " " << i->color.b << " " << " setrgbcolor" << endl;
      for (list<CoordPoint>::iterator p = i->Points.begin(); p != i->Points.end(); p++) {
	if (p == i->Points.begin()) {
	  out << p->x << " " << p->y << " moveto" << endl;
	} else {
	  out << p->x << " " << p->y << " lineto" << endl;
	}	
      }
      out << "stroke" << endl;
    }
  }


  //
  // Done printing inside the plot so we have to grestore from the clipping gsave
  //
  out << "% Done printing inside the frame" << endl;
  out << "grestore" << endl;


  // ==============================================
  // == PRINT TickMarks
  // ==============================================

  out << StandardLineWidth * TicksLineWidth * 0.01 <<  " setlinewidth" << endl;

  if (! Things.XTickMarks.empty()) {
    for (list <TickMark>::iterator i = Things.XTickMarks.begin(); i != Things.XTickMarks.end(); i++) {
      out << "newpath" << endl;
      out << i->axis << " " << -MajorTickLength/2 << " moveto" << endl;
      out << i->axis << " " << MajorTickLength/2 << " lineto" << endl;
      out << "stroke" << endl;
    }
  }

  if (! Things.YTickMarks.empty()) {
    for (list <TickMark>::iterator i = Things.YTickMarks.begin(); i != Things.YTickMarks.end(); i++) {
      out << "newpath" << endl;
      out <<  -MajorTickLength/2 << " " << i->axis << " moveto" << endl;
      out <<  MajorTickLength/2 << " " << i->axis << " lineto" << endl;
      out << "stroke" << endl;
    }
  }


  // ==============================================
  // == PRINT Labels
  // ==============================================

  float XStauch = (Width/Height) * ScaleY2View / ScaleX2View ; // stauch w.r.t canonical
  XStauch = 1.0/XStauch;
  //  float LabelScale = 0.05;

  if (! Things.XTickLabels.empty()) {
    for (list <TickLabel>::iterator i = Things.XTickLabels.begin(); i != Things.XTickLabels.end(); i++) {
      
      out << "% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      out << "%  label at: " << i->axis << "  deci: " << i->decimal << "  exp: " << i->exponent << endl;
      out << "% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      string grace = i->grace;
      if (grace =="") { 
	grace = i->decimal;
	if (i->exponent != "") {
	  grace += "\\S.\\N10\\S " + i->exponent;
	} 
      }
      out << i->axis << " " << - i->offset << " moveto" << endl;
      out << "getboundin" << endl;
      out << parseGraceString(grace, false,i->size,XStauch);
      out <<  "getboundout" << endl;
      out << "getheightandwidth" << endl;
      out << "-0.5 mul" << endl;
      out << "exch" << endl;
      out << "neg" << endl;
      out << "rmoveto" << endl;
      out <<  parseGraceString(grace, true,i->size,XStauch);
    }
  }

  

  //  ===========================
  //        Y-Labels
  //  ===========================

  if (! Things.YTickLabels.empty()) {
    for (list <TickLabel>::iterator i = Things.YTickLabels.begin(); i != Things.YTickLabels.end(); i++) {
      
      out << "% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      out << "%  label at: " << i->axis << "  deci: " << i->decimal << "  exp: " << i->exponent << endl;
      out << "% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

      string grace = i->grace;
      if (grace =="") { 
	string grace = i->decimal;
	if (i->exponent != "") {
	  grace += "\\S.\\N10\\S " + i->exponent;
	} 
      }
  
      out << - i->offset << "  " << i->axis << " moveto" << endl;
      out << "getboundin" << endl;
      out << parseGraceString(grace, false,i->size,XStauch);
      out <<  "getboundout" << endl;
      out << "getheightandwidth" << endl;
      out << "neg" << endl;
      out << "exch" << endl;
      out << "-0.5 mul" << endl;
      out << "rmoveto" << endl;
      out <<  parseGraceString(grace, true,i->size,XStauch);
    }
  }

  //  ===========================
  //     X -   AXIS -Label
  //  ===========================

  if (Things.XAxisLabel.text != "") {
    out << "% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    out << "%   Axis Label for x-axis " << endl;
    out << "% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    // we divide by ScaleY2View such that an offset of 0.1 really means 10% of the paper
    out << Things.XAxisLabel.axis << " " << - Things.XAxisLabel.offset/ScaleY2View  << " moveto" << endl;
    out << "getboundin" << endl;
    out << parseGraceString(Things.XAxisLabel.text, false,Things.XAxisLabel.size,XStauch,true);
    out <<  "getboundout" << endl;
    out << "getheightandwidth" << endl;
    out << "-0.5 mul" << endl;
    out << "exch" << endl;
    out << "neg" << endl;
    out << "rmoveto" << endl;
    out << parseGraceString(Things.XAxisLabel.text, true,Things.XAxisLabel.size,XStauch,true);
  }

 //  ===========================
  //     Y -   AXIS -Label
  //  ===========================

  if (Things.YAxisLabel.text != "") {
    out << "% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    out << "%   Axis Label for y-axis " << endl;
    out << "% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    // we divide by ScaleX2View such that an offset of 0.1 really means 10% of the paper
    out << - Things.YAxisLabel.offset / ScaleX2View << "  " << Things.YAxisLabel.axis << " moveto" << endl;
    out << "gsave" << endl;
    out << "90 rotate" << endl;
    out << "getboundin" << endl;
    out << parseGraceString(Things.YAxisLabel.text, false,Things.YAxisLabel.size,XStauch,false);
    out <<  "getboundout" << endl;
    out << "getheightandwidth" << endl;
     out << "-0.5 mul" << endl;
    out << "exch" << endl;
    out << "neg" << endl;
    out << "rmoveto" << endl;
    out << parseGraceString(Things.YAxisLabel.text, true,Things.YAxisLabel.size,XStauch,false);
    out << "grestore"<<endl;
  }


 
  // ==============================================
  // == PRINT Frame
  // ==============================================
  out << StandardLineWidth * FrameLineWidth * 0.01 <<  " setlinewidth" << endl;
  out << "newpath" << endl;
  out << "0 0 moveto" << endl;
  out << "1 0 lineto" <<endl;
  out << "1 1 lineto" <<endl;
  out << "0 1 lineto" <<endl;
  out <<"closepath" << endl;
  out << "0 setgray"  << endl;
  out <<"stroke" << endl;

  out << "showpage" << endl;

}


string PostscriptPlot::parseGraceString(string grace, bool show,float LabelSize, float XStauch, bool Landscape) {
  
/*!
  Parse label string which is  in a format resembling  xmgr (xmgrace):
  It uses Escape (\) to indicate a command:
  
  \1: charset 1
  \2: charset 2
  \3: symbol charset
  \s: subscript on
  \S: superscript on
  \N: normal
*/
 

  cout << "parsegrace XSTAUCH: " << XStauch << endl;

  string out;
  float ptscale = 1.0;
  
  // We scale the shift such that a shift of 0.1 really means 10% of the actual paper
  float ShiftScale = 1.0/ScaleY2View;
  if (Landscape) ShiftScale = 1.0/ScaleX2View;  

  //vector<double> AnalyzeThis::parsePrintString(string& p, ostream &out, double ptscale, bool dummywrite) {

  list <GraceLabelHelper> labels;
  //  vector<double> heights(3,0.0);  // 0: total, 1: below base: 2 above upper edge
  //fontbehaviour former = normal, next = normal;
  double scale[]={1.0,0.8,0.8};
  for (int i =0; i < 3; i++) scale[i] *= LabelSize;

  for (int i = 0;i<3;i++) 
  cout << "Scale is : " << scale[i] << endl;

  //  cout << "parsing: " << p << endl;
  
  bool flush=false;
  bool escaped = false;
  bool hassuper=false,hassub=false;

  GraceLabelHelper current,next;
  current.scale = scale[0];
  current.family = "family1";
  current.shift = 0;
  current.stopwagon = false;
  next = current;

  //string text;
  int count =0;
  for (unsigned int i =0; i < grace.size() ; i++) {
    char c = grace.at(i);
    if (c == '\\') escaped = true;
    else { 
      if (escaped) {
	//	cout << "escaped! " << c << "  int: " << (int)c << endl;
	escaped = false;
	if (c == 's') { next.scale = scale[1]; next.shift = 0.4*ptscale*scale[1]; next.stopwagon = true; hassub=true;   }
	if (c == 'S') { next.scale = scale[2]; next.shift = -0.4*ptscale*scale[2]; next.stopwagon = true; hassuper= true;}
	if (c == 'N') { next.scale = scale[0]; next.shift = 0; next.stopwagon=false; }
	if (c == ' ') { 
	  cout << "CAUGHT ESCAPED SPACE"<<endl;
	  cout << "current.text is: " << current.text << endl;
	  if (current.stopwagon) cout << "stopwagon is on" << endl;
	  current.stopwagon = false; next.text += " ";
	}
	if (c == '1') next.family = "family1";
	if (c == '2') next.family = "family2";
	if (c == '3') next.family = "family3";
	flush = true;
      } else {
	bool use=true;
	if (c == '(') { current.text += "\\050";  use =false;}
	if (c == ')') { current.text += "\\051";  use =false;}
	if (use) current.text += c;
      }
    }
    if (flush) { 
      flush = false;
      if (current.text != "") labels.push_back(current);
      current = next;
      next.text = "";
    }   
  }
  
  //  heights[0] = ptscale;
  //  if (hassub) { heights[0] += 0.5*ptscale*scale[1]; heights[1] = 0.7*ptscale*scale[1];}
  //  if (hassuper) { heights[0] += 0.5*ptscale*scale[2];  heights[2] = 0.5*ptscale*scale[2];}
  //  if (dummywrite) return heights;
  // final flush
  if (current.text != "") labels.push_back(current);
  list<GraceLabelHelper>::reverse_iterator k = labels.rbegin();

  float x_stauch = XStauch;
  float y_stauch = 1.0;

  if (! Landscape) {
    y_stauch = x_stauch;
    x_stauch = 1.0;
  }

  while (k != labels.rend()) {
    //    cout << "For text (" << k->text <<")  we assume stopwagon: " << k->stopwagon << endl;
    cout << "For text  (" << k->text <<")  we assume scaling: " << k->scale << "  the values are then 1/scale: " << 1.0/k->scale << "  stauch/scale: " << XStauch / k->scale << "   scale/stauch: " << k->scale / XStauch << endl;
    if (k->stopwagon) out += " 1 "; else out += " 0 ";
    out += " 0  " + double2string(k->shift) + " ";
    out += double2string(x_stauch / k->scale,6) + " " + double2string(y_stauch/k->scale,6);
    out += " (" + k->text + ") ";
    if (show) out += " 1 "; else out += " 0 ";
    out += " ("+ k->family + ") ";
    out += double2string(k->scale/x_stauch,6) + " " + double2string(k->scale/y_stauch,6);
    out += " 0 " + double2string(-k->shift); 
    if (k->stopwagon) out += " 1 "; else out += " 0 ";
    out += "\n";
    count++;
    k++;
  }
  out += " " + int2string(count) += "\n";
  out += "parsetext\n";
  //  return heights;
  return out;

}



/*!
  Parses string label which should use latex syntax and
  produces xmgrace-like format that is actually used for
  the postscript output.
*/
string PostscriptPlot::parseTexLabel(string label) {
  string parsed;
  string memory="";
  label += " ";
  //  list <LabelState> state;
  // state.push_back(state(1,0)); 
  bool escape = false;
  bool returntonormal = false;
  bool appendnormal = false;
  bool addspace=false;
  int bracket =0;
  
  int SaveFont[100];
  SaveFont[0] = 1;
  cout << "parsing: " << label << endl;


  for (unsigned int i =0; i < label.size() ; i++) {
    char c = label.at(i);
  
    //cout << c << "  " << (int) c << "   " << escape << endl;
    if (c == '\\') {  
      escape = true; memory =""; 
      if (appendnormal) {
	appendnormal=false;
	parsed += "\\ \\N";
      }
    }  else {
      bool checkmemory = false;
      string append = "";
      if (c == '^') { escape = false; append = "\\S"; appendnormal=false;checkmemory = true; returntonormal=true; addspace=false; }
      if (c == '_') { escape = false; append = "\\s";  appendnormal=false;checkmemory = true; returntonormal = true; addspace=false;}
      if (addspace) { parsed+="\\ \\N"; addspace=false;}
      if (c == ' ') { 
	escape=false; checkmemory = true; 
	if (bracket==0) append = "\\ "; // signals "no stop wagon"
      }   
      if (c == '~') { 
	escape=false; checkmemory = true; 
	if (bracket==0) append = ""; // signals "no stop wagon"
      }   
      if (c == ')' || c=='(' ||  c=='[' || c==']') {
	escape=false; checkmemory = true; 
	append = "";
	if (appendnormal) { append = "\\ \\N"; appendnormal = false;}
	append += c;
	cout << "Habe runde klammer: append ist: " << append << endl;
	if (returntonormal) cout << "returntonormal" << endl;
      } 
     
      if (c == '{') { 
	escape=false; checkmemory = true; 
	bracket++; 
	SaveFont[bracket] = SaveFont[bracket-1];
      }
      if (c == '}') { 
	escape=false; checkmemory = true; 
	bracket --; 
	cout << "Klammer zu" << endl;
	cout << "parsed: " << parsed << endl;
	cout << "append: " << append << endl;
	cout << "memory: " << memory << endl;
	if (returntonormal) {
	  // this situation arises, if we have something like "log_{10}..." now the question
	  // is: should we insert a \ \N to return to normal and have no stopwagon for grace
	  // or is there maybe something like "log_{10}^2" and we should keep the stopwagon ?
	  // so what we do is use a new flag, called addspace
	  addspace=true;
	  returntonormal = false;
	  //	  appendnormal = true; returntonormal = false;
	} 
	if (bracket < 0) throw Bad_Error("PostscriptPlot::parseTexLabel() one closing curly bracket } too much");
	if (SaveFont[bracket] != SaveFont[bracket+1]) { parsed += "\\" + int2string(SaveFont[bracket]) ; cout << "\\" << SaveFont[bracket] << endl;}
      } 
      if (checkmemory && memory != "") {
	cout << "COMPARING: " << memory << endl;
	if (memory == "Omega") parsed += "\\3W\\1";
	if (memory == "omega") parsed += "\\3w\\1";
	if (memory == "Sigma") parsed += "\\3S\\1";
	if (memory == "sigma") parsed += "\\3s\\1";
	if (memory == "Gamma") parsed += "\\3G\\1";
	if (memory == "gamma") parsed += "\\3g\\1";
	if (memory == "Alpha") parsed += "\\3A\\1";
	if (memory == "alpha") parsed += "\\3a\\1";
	if (memory == "Beta") parsed += "\\3B\\1";
	if (memory == "beta") parsed += "\\3b\\1";
	if (memory == "Delta") parsed += "\\3D\\1";
	if (memory == "delta") parsed += "\\3d\\1";
	if (memory == "Rho") parsed += "\\3R\\1";
	if (memory == "rho") parsed += "\\3r\\1";
       	if (memory == "Pi") parsed += "\\3P\\1";
       	if (memory == "pi") parsed += "\\3p\\1";
	if (memory == "nu") parsed += "\\3n\\1";
	if (memory == "tau") parsed += "\\3t\\1";
	if (memory == "Tau") parsed += "\\3T\\1";
	if (memory == "epsilon") parsed += "\\3e\\1";
	if (memory == "zeta") parsed += "\\3z\\1";
	if (memory == " ") parsed += "\\ ";
	if (memory == "quad") parsed += "\\ \\ \\ \\ \\  ";
	if (memory == "rm") { parsed += "\\2"; SaveFont[bracket] = 2; }
	memory = "";
      }
      if (checkmemory) parsed += append;
      if (appendnormal) { appendnormal = false; parsed += "\\N"; }
      if (escape) memory += c; 
    
      //     cout << "memory: " << memory << endl;
      if (!checkmemory) {
	if (!escape && c != '{' && c != '}') parsed += c;
	if (returntonormal && bracket == 0) {
	  returntonormal = false;
	  appendnormal = true;  // don't append yet, cause maybe there needs to be a space
	  // first to indicate a wagon move to the right (no stopwagon)
	  // parsed += "\\N";
	}
      }
    }
  }
  cout << "Tex parsed: "<< parsed << endl;
  return parsed;
}

/*!  
  Little helper for drawArea. Acronym for Multiply And Cut
*/ 
int PostscriptPlot::mac(float x,float scale) {
  x *= scale;
  return (int)(rint(x));
}

void PostscriptPlot::InterpolateColor(Color &c,Color &up, Color& down, float z, float Min,float Max) {
  float r = (z - Min)/(Max-Min);  // relative value
  c.r = up.r*r  + down.r*(1-r);
  c.g = up.g*r  + down.g*(1-r);
  c.b = up.b*r  + down.b*(1-r);
  if (c.r > 1.0) c.r = 1.0;
  if (c.g > 1.0) c.g = 1.0;
  if (c.b > 1.0) c.b = 1.0;
  //  cout << "Interp: " << z << " r: " << c.r << " g: " << c.g << " b: " << c.b << endl;
}
