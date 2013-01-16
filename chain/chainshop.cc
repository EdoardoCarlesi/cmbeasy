#include "chainshop.h"
//#include "spline.h"
#include "miscmath.h"
//#include "data.h"
//#include "controlpanel.h"
//#include "cl.h"
#include "minmax.h"
#include "anchor.h"
#include "spline.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_errno.h>

#include <fstream>
#include <algorithm>
#include <limits>
#include <numeric>

using namespace std;

#ifdef MACOSX_PANTHER
extern "C" int isnan(double); 
extern "C" int isinf(double); 
#endif 


#define LARGESTMULTIPLE 1000*1000
#define MAGICANALYZE 37551

//silence valgrind (?)
static int intread(istream& i) {
  int s=0;
  i.read((char*) &s,sizeof(int));
  return s;
}

ChainShop::ChainShop() {
//debug only
#ifndef PRERELEASE
  mAveragingType=UseTotal;
  mDistributionType=UseWeight;
#endif
}

ChainShop::~ChainShop() {

}
/*!
  Use the map<float, map<float,float> > m which already gives the likelihood on a grid (x,y)
  to draw 1 and 2 sigma likelihood contours. Output is in postcript to a file specified in the
  Printing structure print.
  
  If returnlines is true, then don't write to file but return the lines of the sigma regions
  in a struct RasterizReturn. This is used by the graphical front-end.

*/
list<Block>* ChainShop::rasterize(Grid& m,int *Progress) {
  list<Block>* blocks = new list<Block>();
  if (Progress) *Progress = 20;
  //  float range_x = (max_x - min_x);
  // float range_y = (max_y - min_y);


  int max_regions=100;
  vector<double> DeltaChi2(max_regions); 
  vector<float> Sigma(max_regions);
  vector<bool> border(max_regions,false); 
    
  int regions =0;
  Miscmath misc;
  for (double sig = 0.2; sig <= 3.001; sig += 0.1) { 
    //cout << "sig: " << sig << endl;
    DeltaChi2[regions] = misc.DeltaChi2CorrespondingToSigma(2,sig);
    Sigma[regions] = sig;
    //    cout << "SIGMA: " << sig << "  chi2: " << DeltaChi2[regions] << endl;
    for (int snap=1; snap <= 3 ; snap++) 
      border[regions] = border[regions] || Miscmath::relError(sig,snap) < 1e-2;  
    regions++;
  }
  double a,b;
  vector< map<float, map<float, bool> > > findlines(regions);
  for (Grid::iterator i = m.begin(); i != m.end(); i++) {
    a = i->first;
    for (map<float,float>::iterator j = i->second.begin(); j != i->second.end(); j++) {
      b = j->first;
      for (int r = regions-1; r >= 0; r --) {
	if (j->second < DeltaChi2[r]) { 
	  findlines[r][a][b] = true;
	}
      }
    }
  }      
  

  
  
  if (Progress) *Progress = 60;
  
  int  Size = findlines.size();
  int count = 0;
  for (unsigned int sigma =0; sigma < findlines.size(); sigma++) {
    Block b;
    if (Progress) *Progress = 60 + (int)(20 * count++ / Size);
    for ( map<float, map<float, bool> >:: iterator i = findlines[sigma].begin(); i != findlines[sigma].end(); i++) {
      float column = i->first;
      map<float, bool>::iterator j = i->second.begin();
      float lastpos = j->first;
      float left = lastpos;
      while (++j != i->second.end()) {
	float pos = j->first;
	if (pos - lastpos > step_y*1.001) {  // ha! here we make a jump
	  b.lines[column].push_back(CoordPoint(left, lastpos));
	  left = pos;
	}
	lastpos = pos;
      }
      b.lines[column].push_back(CoordPoint(left,lastpos));      
    }
    Color OneSigma(1,1,1);
    //    if (MaxSigma < 2) OneSigma = Color(0,0,0); 
    float sig = Sigma[sigma];
    if (sig > 0.999 && sig < 1.001) { b.PenColor = OneSigma; b.UsePen=true;}
    if (sig > 1.999 && sig < 2.001) { b.PenColor = Color(0,0,0); b.UsePen=true; }
    if (sig > 2.999 && sig < 3.001) { b.PenColor = Color(0,0,0); b.UsePen=true; } 
    b.FillColor = shade(sig / 3.0);
    b.Width = step_x;
    blocks->push_back(b);
  }    
  return blocks;
}

static ostream& operator<<(ostream& os, const CoordPoint& p)
{
  os << "(" << p.x << "/" << p.y << ")";
  return os;
}

list<Block>* ChainShop::rasterize_Bayesian2(Grid&  m, const ChainShop::GridInfo& gridInfo, int *Progress) {
  list<Block>* blocks = new list<Block>();
  if (Progress) *Progress = 20;
  //  float range_x = (max_x - min_x);
  // float range_y = (max_y - min_y);


  int max_regions=100;
  vector<double> relArea(max_regions);
  vector<float> Sigma(max_regions);
  vector<bool> border(max_regions,false); 
    
  int regions =0;
  Miscmath misc;
  for (double sig = 0.2; sig <= 3.001; sig += 0.1) { 
    relArea[regions] = misc.fractionCorrespondingToSigma(sig);
    Sigma[regions] = sig;
    for (int snap=1; snap <= 3 ; snap++) 
      border[regions] = border[regions] || Miscmath::relError(sig,snap) < 1e-2;  
    regions++;
  }
  double totalMult=0;
  double maxMult=0;
  map<double, double> areaMap;
  for (map<float, map<float,float> >::iterator i = m.begin(); i != m.end(); i++) {
    //a = i->first;
    for (map<float,float>::iterator j = i->second.begin(); j != i->second.end(); j++) {
      //b = j->first;
      double mult = j->second;
      totalMult += mult;
      areaMap[mult]+=mult;
      maxMult = max(maxMult, mult);
    }
  }

  vector<double> partialSums(areaMap.size());
  map<double, double>::const_iterator it, end = areaMap.end();
  unsigned long c=0;
  for (it=areaMap.begin(); it!=end; ++it) {
    partialSums[c++] = it->second;
  }

  partial_sum(partialSums.rbegin(), partialSums.rend(), partialSums.rbegin());
  map<double, double>::iterator mit, mend = areaMap.end();
  for (mit = areaMap.begin(), c=0; mit!=mend; ++mit) {
    (*mit).second=partialSums[c++];
  }


  double a,b;
  vector< map<float, map<float, bool> > > findlines(regions);
  for (map<float, map<float,float> >::iterator i = m.begin(); i != m.end(); i++) {
    a = i->first;
    for (map<float,float>::iterator j = i->second.begin(); j != i->second.end(); j++) {
      b = j->first;
      for (int r = regions-1; r >= 0; r --) {
        double mult = j->second;
        if ((areaMap[mult]/totalMult)<relArea[r]) {
          findlines[r][a][b] = true;
        }
      }
    }
  }
  cout << "**************total Mult: " << totalMult << endl;

  if (Progress) *Progress = 60;

  int  Size = findlines.size();
  int count = 0;
  for (unsigned int sigma =0; sigma < findlines.size(); sigma++) {
    Block b;
    if (Progress) *Progress = 60 + (int)(20 * count++ / Size);
    for ( map<float, map<float, bool> >:: iterator i = findlines[sigma].begin(); i != findlines[sigma].end(); i++) {
      float column = i->first;
      map<float, bool>::iterator j = i->second.begin();
      float lastpos = j->first;
      float left = lastpos;
      while (++j != i->second.end()) {
        float pos = j->first;
        if (pos - lastpos > step_y*1.001) {  // ha! here we make a jump
          b.lines[column].push_back(CoordPoint(left, lastpos));
          left = pos;
        }
        lastpos = pos;
      }
      b.lines[column].push_back(CoordPoint(left,lastpos));      
    }
    Color OneSigma(1,1,1);
    //    if (MaxSigma < 2) OneSigma = Color(0,0,0); 
    float sig = Sigma[sigma];
    if (sig > 0.999 && sig < 1.001) { b.PenColor = OneSigma; b.UsePen=true;}
    if (sig > 1.999 && sig < 2.001) { b.PenColor = Color(0,0,0); b.UsePen=true; }
    if (sig > 2.999 && sig < 3.001) { b.PenColor = Color(0,0,0); b.UsePen=true; } 
    b.FillColor = shade(sig / 3.0);
    b.Width = step_x;
    blocks->push_back(b);
  }
  return blocks;
}

/*!
  Scan monte carlo file filename and determine minum and maximum values of collumns
  pos_x and pos_y (starting from 0)
*/
vector<float> ChainShop::autoScale(string filename,int pos_x, int pos_y) {
  int line=0;
  float min_x = 1e20, max_x = -1e20, min_y = 1e20, max_y=-1e20,x=0,y=0;
  float remind[100];
  ifstream in(filename.c_str());
  if (!in) throw Bad_Error("ChainShop::autoScale()  file not found");
  if (read<int>(in) != MAGICANALYZE) throw Bad_Error("ChainShop::autoscale..() Not a valid .mcc file");
  int size = read<int>(in); // reading number of columns 
  cout << "expected columns: " << size << endl;
  while (! in.eof()) {
    for (int i = 0; i < size-1; i++) {
      remind[i] = read<float>(in);
      if (in.fail()) break;
      if (i == pos_x) x = remind[i];
      if (i == pos_y) y = remind[i];
    }
    int multiple = read<int>(in);
    if (in.fail()) break;
    line++;
    min_x = min(min_x,x);
    max_x = max(max_x,x);
    min_y = min(min_y,y);
    max_y = max(max_y,y);
     

     if (multiple < 1 || multiple > LARGESTMULTIPLE) {
       cout << "multiple: " << multiple << " wrong in line: " << line <<  endl;
       throw Bad_Error("ChainShop::autoscale..() multiple suspiciously larger or small");    
     }
  }

  cout << "Total line count " << line << endl;

  vector<float> v(4);
  v[0] = min_x;
  v[1] = max_x;
  v[2] = min_y;
  v[3] = max_y;
  return v;


}

 

/*! 
  Read in monte carlo file filename, using columns pos_x and pos_y (starting from 0) as x and y axis.
  min_x, min_y, max_x,max_y determine the range of parameters for x and y. For the 2-dim plane spanned
  by min_x, max_x etc, fill the map s with the (relative) likelihood of parameters (x,y). 
  It uses int steps in x and y directions (i.e. the size of map s is steps^2). For each parameter point in the
  file, it smears out the discrete point using a gaussian window with width range / steps * SmearFactor.
  
  While s[x][y] holds the main result (i.e. the likelihood surface), it also writes the discrete points to a file 
  called "gdh.dat" which you can for instance plot by gnuplot: "> plot gdh.dat" should do the job in gnuplot.
*/
void ChainShop::get2DimMCDistribution(string filename, Grid& s,
                                      int pos_x, int pos_y, int pos_likeli,
                                      GridInfo& gridInfo, int bins, float SmearFactor,
                                      int* Progress, bool GridFile)
{
  if (Progress) *Progress = 0;
  ifstream in(filename.c_str());
  if (!in) throw Bad_Error("get2dimMCDistribution():: file not found");

  float Min_x = 1e50;
  float Max_x = -1e50;
  float Min_y = 1e50;
  float Max_y = -1e50;
  float x=0, y=0;
  float lnLike=0;

  list<vector<float> >  values;
  int line=0;
  bool ignore=false;
  float remind[100];
  //  char c;
  int multiple=0;
  if (read<int>(in) != MAGICANALYZE) throw Bad_Error("ChainShop::get2dim..() Not a valid .mcc file");
  int size = read<int>(in); // reading number of columns 
  //  cout << "expected columns: " << size << endl;
  double debug_totalMult=0;
  while (!in.eof()) {
    for (int i = 0; i < size-1; i++) {
      remind[i] = read<float>(in);
      if (in.fail()) {ignore = true; break;}
      if (i == pos_x) x = remind[i];
      if (i == pos_y) y = remind[i];
      if (i == pos_likeli) lnLike = remind[i];
    }
    line++;
    multiple = intread(in);
    debug_totalMult += multiple;
    //cout << "line: [ " << line << " : " << multiple <<"]  ";
    //for (int k = 0; k < size-2; k++) cout << remind[k] << "  ";
    //cout << multiple << endl; 

    if ((multiple < 1 || multiple > LARGESTMULTIPLE) && !ignore) {
      cout << "multiple: " << multiple << endl;
      throw Bad_Error("ChainShop::get2dim..() multiple suspiciously larger or small");    
    }

    if (!ignore) {
      vector<float> v(4);
      v[0] = x;
      v[1] = y;
      v[2] = lnLike;
      v[3] = multiple;
      values.push_back(v);

      Min_x = min(x,Min_x);
      Min_y = min(y,Min_y);
      Max_x = max(x,Max_x);
      Max_y = max(y,Max_y);
    }
  }
  //cout << "***debug val of total mult: " << debug_totalMult << endl;

  step_x = (Max_x - Min_x)/bins;
  step_y = (Max_y - Min_y)/bins;
  min_x = Min_x + 0.5*step_x;
  min_y = Min_y + 0.5*step_y;
  max_x = min_x + (bins-1)*step_x;
  max_y = min_y + (bins-1)*step_y;

  gridInfo.min_x = min_x;
  gridInfo.min_y = min_y;
  gridInfo.max_x = max_x;
  gridInfo.max_y = max_y;
  gridInfo.Min_x = Min_x;
  gridInfo.Min_y = Min_y;
  gridInfo.Max_x = Max_x;
  gridInfo.Max_y = Max_y;
  gridInfo.step_x = step_x;
  gridInfo.step_y = step_y;


  if (GridFile) {
    //    cout << "WRITING GRID FILE: " << endl;
    unsigned int ValSize = values.size();
    ofstream histo("gdh.dat");
    int cnt=0;
    for (list<vector<float> >::iterator i = values.begin(); i != values.end(); i++,cnt++) {
      if (Progress) *Progress = (10*cnt) / ValSize;
      vector<float> &t = *i;
      histo << t[0] << " " << t[1] << endl;
    }
  }

  DistributionType type = UseWeight;
  AveragingType averagingType = UseTotal;
#ifndef PRERELEASE
  type = mDistributionType;
  averagingType = mAveragingType;
#endif

  // ==============================
  // first initialize the map[x][y] = 0
  // ==============================
  for (int i = 0; i < bins; i++) {
    for (int j = 0; j < bins; j++) {
      float a = min_x + i*step_x;
      float b = min_y + j*step_y;
      s[a][b] = 0;
    }
  }

  int count =0;
  unsigned int Size = values.size();

  float smear_x = SmearFactor * step_x;
  float smear_y = SmearFactor * step_y;
  float sx = -1.0/(2*smear_x*smear_x);  // gaussian width ^2 for smearing
  float sy = -1.0/ (2*smear_y*smear_y); 
  int dp = (int)rint(5*SmearFactor);  // size of square
  int dq = (int)rint(5*SmearFactor);

  vector<vector<float> > smearedOutWeights;
  if (SmearFactor != 0.0) {
    float f;
    float total=0;
    int vecSize1=2*dp+2;
    int vecSize2=2*dq+2;
    smearedOutWeights = vector<vector<float> >(vecSize1, vector<float>(vecSize2, 0));
    for (int pp=-dp; pp<=dp; pp++) {
      for (int qq=-dq; qq<=dq; qq++) {
        float chi2=pp*pp+qq*qq;
        f=exp(-0.5*chi2/SmearFactor/SmearFactor);
        smearedOutWeights[pp+dp][qq+dq]=f;
        total+=f;
      }
    }
    for (int pp = -dp; pp <= dp; pp++) {
      for (int qq = -dq; qq <= dq; qq++) {
        smearedOutWeights[pp+dp][qq+dq]/=total;
      }
    }
  }

  vector<vector<float> > binnedVals(bins, vector<float>(bins));
  vector<vector<float> > binCount(bins, vector<float>(bins));
  // run through the list of (x,y) parameter pairs
  double debug_totalMult2=0;
  for (list<vector<float> >::iterator i = values.begin(); i != values.end(); i++) {
    if (Progress) *Progress = 10 + (90*count) / Size;
    //  cout << "Progress: " << Progress << endl;
    count++;
    x = (*i)[0];
    y = (*i)[1];
    int multiple = (int)rint((*i)[3]);
    debug_totalMult2 += multiple;
    float current_lnLike = (*i)[2];

    float a,b; // = min_x; 

    int p = (int)floor((x - Min_x)/step_x);  // discrete point closest to x 
    int q = (int)floor((y - Min_y)/step_y); // discrete point closest to y

    bool inside=true;
    if (p < 0 || p > bins-1 || q < 0 || q > bins-1) inside=false; // happens only for the very outmost points, i.e for 4 in total

    if (inside) {
      a = min_x + p*step_x;
      b = min_y + q*step_y;

      if (a > max_x + 0.01*step_x   || b  > max_y + 0.01*step_y) {
        cout <<"x: " << x << " y: " << y << endl;
        cout << "a: " << a << " b: " << b << endl;
        cout << "min_x: " << min_x << "  max_x: " << max_x << "  bins: " << bins << " step_x: " << step_x<<endl;
        cout << "min_y: " << min_y << "  max_y: " << max_y << "  bins: " << bins << " step_y: " << step_y<<endl;
        cout << "Min_y: " << Min_y << " Max_y: " << Max_y << endl;
        throw Bad_Error("stop");
      }

#ifndef PRERELEASE
      if (SmearFactor != 0.0 && (type!=UseWeight || averagingType!=UseTotal) )
      {
        throw Bad_Error("ChainShop::get2dimMCDistribution() - if you use smearing, use DistributionType=UseWeight and averagingType=UseTotal.");
      }
#endif

      if (SmearFactor != 0.0) {
        float val = multiple;
        if (type==UseChi2)
          val = current_lnLike;

        int dp = (int)rint(5*SmearFactor);  // size of square under investigation
        int dq = (int)rint(5*SmearFactor);

        if (dp ==0 || dq == 0) throw Bad_Error("smearing is smaller than stepsize !");

        int start_p = max(p-dp,0);
        int stop_p = min(p+dp,bins-1);
        int start_q = max(q-dq,0);
        int stop_q = min(q+dq,bins-1);
        for (int pp = start_p; pp <= stop_p; pp++) {
          for (int qq = start_q; qq <= stop_q; qq++) {
            a = min_x + pp*step_x;
            b = min_y + qq*step_y;
            float weight = smearedOutWeights[p-pp+dp][q-qq+dq];
            binnedVals[pp][qq] += weight*val;
          }
        }
      }  else {
        float val;
        val = multiple;
#ifndef PRERELEASE
        if (type==UseChi2)
          val = current_lnLike;
        if (averagingType==UseMaximum) {
          if (binnedVals[p][q]==0) {
            binnedVals[p][q] = val;
          } else {
            binnedVals[p][q] = max(binnedVals[p][q], val);
          }
        } else
#endif
          binnedVals[p][q] += val;

#ifndef PRERELEASE
        if (averagingType==UseAverage || type==UseChi2) {
          ++binCount[p][q];
        }
#endif
      }
    }
  }

  //cout << "**** totalMul2: " << debug_totalMult2 << endl;

  float MinChi2=1e100;

  for (int i = 0; i < bins; i++) {
    for (int j = 0; j < bins; j++) {
      float val = binnedVals[i][j];
      float a = min_x + i*step_x;
      float b = min_y + j*step_y;

#ifndef PRERELEASE
      if (type==UseChi2) {
        val=exp(-0.5*val);
        if (binCount[i][j]==0)
          val = std::numeric_limits<double>::infinity();
        if (val < MinChi2) {
          MinChi2 = val;
          gridInfo.minCoords = make_pair(a,b);
        }
      }
      if (averagingType==UseAverage && binCount[i][j]!=0) {
        val/=binCount[i][j];
      }
#endif
      if (isinf(val) || isnan(val))
        throw Bad_Error("ChainShop::get2dimMCDistribution() - encountered unexpected nan or inf");
      s[a][b]=val;
    }
  }

  if (type==UseChi2) {
    Grid::iterator it, end = s.end();
    for (it=s.begin(); it!=end; ++it) {
      map<float, float>::iterator it1, end1 = it->second.end();
      for (it1=it->second.begin(); it1!=end1; ++it1) {
        (*it1).second -= MinChi2;
      }
    }
  }
  //cout << "MinChi2: " << MinChi2 << endl;
}

/*!
  Read in the files given in vector<string> infiles (Filename is without .dat) and
  parse the ascii data in these raw montecarlo files.

  There are two versions, one is obsolete (for the public at least), indicated by
  OldVersion = true:

  For each input file, it discards
  the first int discard models. Please note, that the time spent at each parameter 
  point, i.e. the weigh of that parameter point is counted towards discard. For example,
  if the first 12 models investigated have a total multiplicity of 199, they will be discarded,
  the next one (and all following in that input file) is not.

  The new version uses the new raw data file format, OldVersion = false:

  Scan all input files for parameter entry with zero multiplicity. The montecarlo driver uses
  zero multiplicity to indicate the transition from burn in to frozen stepsize.
  For this (new) version, int discard is not used.

  \param infiles Vector of strings with the names of the input files
  \param outfile Name of output file without trailing .mcc and .dat
  \param discard Number of chain points to discard at the beginning. Obsolete if OldVersion=false
  \param Progress Pointer to integer. Can be used to monitor progress
  \param OnlyBinary If no human readable merged .dat file should be produced, just the binary .mcc
  \param OldVersion Switch to an old version. The new version is the one implemented thoughout here

*/


double ChainShop::distillChain(vector<string> infiles, string outfile, int* Progress, int discard, bool OnlyBinary,bool OldVersion) {
  string outfileName=outfile+".dat";
  vector<string>::iterator it = find(infiles.begin(), infiles.end(), outfileName);
  if (it!=infiles.end()) { // outfile would overwrite one of the infiles (if mcc name is equal to a chainfile name modulo .dat/.mcc
    outfileName = outfile+"_combined_for_destillation.dat";
  }
  ofstream out;
  if (! OnlyBinary) out.open(outfileName.c_str());
  ofstream binout((outfile+".mcc").c_str());

  // sneak preview of the input: we find out how many columns there are

  //cout << "distilling chain: " << infiles[0] << endl;

  ifstream sneak((infiles[0]).c_str());
  if (!sneak) throw Bad_Error("distillChain():: file not found: "+infiles[0]);
  char buffer[1000];
  for (int i =0; i < 999; i++) buffer[i]=0;
  //  double preview;
  sneak.getline(buffer,999);
  
  int i = 0;
  int number=0;
  int columns =0;
  //cout << "sneak: " << sneak << endl;
  for (;;) {
    char c = buffer[i++];
    if (c == '\n' || c ==  0) break;
    if (c == ' ')  number=0; else number++;
    cout << "char["<<i <<"] " << (unsigned int) c << "  " << c << " number: " << number << endl; 
    if (number==1) columns++;
    if (i > 900) throw Bad_Error("distillChain():: row exceeds 900 chracacters, row too long");
  } 
  //  if (number) columns++; // count the last one
  cout << "NUMBER OF COLUMNS " << columns << endl;
  

  long total_multiple = 0;
  long lines = 0;

  write<int>(binout,MAGICANALYZE);
  write<int>(binout, columns);
  
  int PerFile = 100/infiles.size();  // for Progress
  int LinesPerFile=40000; // for progress, stores number of lines in first file (initialized here with sensible number)
  int MomentaryLines=0; // for progress
  int UpdateCount=0; // for progress

  for (unsigned int f = 0; f < infiles.size(); f++) {
    MomentaryLines=0;
    ifstream in((infiles[f]).c_str());
    //ofstream out((filename + "_distilled.dat").c_str());
    if (!in) throw Bad_Error("distillChain():: file not found");
    out.setf(ios::scientific);
    vector<double> data(columns), store(columns);
    bool virgin = true;
    int count = 0;
    int CurrentLine=0;
    
    // cout << "opening: " << infiles[f] << endl;
    if (OldVersion) { // OLD DATA FORMAT VERSION
      while (! in.eof()) {
	for (int i =0; i < columns; i++) {
	  store[i] = data[i];
	  in >> data[i];
	}
	int multiple = (int)rint(data[columns-1]);
	if (!virgin) count += multiple;
	if (!virgin && count > discard) {
	  for (int i = 0; i < columns - 1; i++) {
	    //	  float tmp = store[i];
	    if (!OnlyBinary) out << store[i] << " ";
	    write<float>(binout, store[i]);
	  }
	  if (!OnlyBinary) out << multiple+1 << endl;
	  int tmpint = multiple+1;
	  write<int>(binout, tmpint);
	  total_multiple += multiple+1;
	  lines++;
	} 
	virgin = false;
      }
    } else { // NEW DATA FORMAT VERSION
      //cout << "NEW VERSION" << endl;
      while (! in.eof()) {
	CurrentLine++;
	for (int i =0; i < columns-1; i++) {
          in >> data[i];
          //cout << "read " << data[i] << endl;
        }
	float mul;
	int multiple; 
	in >> mul;	//	in >> multiple; (was unsuited if output in decimal format)
        //cout << "read mult: " << mul << endl;
	multiple = (int)mul;
	if (in.fail()) break;
	//	if (virgin) cout << "virgin -- ";
	// cout << "multiple: " << multiple << endl; 
	if (!virgin) {
	  for (int i = 0; i < columns - 1; i++) {
	    if (!OnlyBinary) out << data[i] << " ";
	    write<float>(binout, data[i]);
	  }
	  if (multiple < 1 || multiple > LARGESTMULTIPLE) {
	    cout << "SUSPICION: " << infiles[f] << "  in line: " << CurrentLine << " value: " << multiple << endl;	    
	  }

	  if (!OnlyBinary) out << multiple << endl;
	  write<int>(binout, multiple);
	  total_multiple += multiple;
	  lines++;
	  MomentaryLines++;
	  UpdateCount++;
	} else {
	  lines++; //remove
	  total_multiple+= multiple; //remove
	}
	if (multiple == 0) {
	  if (virgin == false) throw Bad_Error("ChainShop::distill() burn-in ended more than once. Corrupt data file ?");
	  virgin = false;  // from now on, note the lines
	  cout << "total_multiple: "<< total_multiple << " for : " << lines << endl;
	  cout << "BURN in average stay: " << (float)total_multiple / (float)lines << endl; //remove
	  total_multiple = 0;  //remove
	  lines = 0; //remove 
	}
	if (Progress && UpdateCount == 100) {
	  *Progress = f*PerFile + (MomentaryLines*PerFile)/LinesPerFile;
	  UpdateCount = 0;
	}
      }
      if (f==0) LinesPerFile = lines;
      if (virgin) throw Bad_Error("ChainShop::distillChain() New format raw data file, yet no line with zero multiplicity\nindicating start of frozen stepsize");
    }
    cout << "total_multiple: "<< total_multiple << " for : " << lines << endl;
    cout << "AFTER Burn in stay: " << (float)total_multiple / (float)lines  << endl;
  }
  return  (float)total_multiple / (float)lines; // average stay at model
  //cout << "Average stay at a model: " << (float)total_multiple / (float)lines << endl;
}

/*!
  Write 1-Dim distribution to file
*/
void ChainShop::dumpDistribution(string filename, map<double,double>& d) {
  ofstream out(filename.c_str());
  for (map<double,double>::iterator i = d.begin(); i != d.end() ; i++) 
    out << i->first << "  " << i->second << endl;
}

void ChainShop::likeli2chi2(map<double,double>& m) {
  for (map<double,double>::iterator i = m.begin(); i != m.end(); i++) i->second = -2.0*log(i->second);
}
void ChainShop::chi2likeli(map<double,double>& m) {
  for (map<double,double>::iterator i = m.begin(); i != m.end(); i++) i->second = exp(-0.5*i->second);
}

/*
void ChainShop::perfectSpline(map<double,double>& data, map<double,double> &out,int MaxNr ) {

  ofstream perf("perfect.dat");
  double min=1e100, maxi=-1e100;
  for (map<double,double>::iterator i = data.begin(); i != data.end(); i++) {
    double y = i->second;
    if (y > maxi) maxi = y;
    if (y < min)  min = y;
  }
  
  double range = maxi-min;
   
  map<double,double> subsample;

  for (map<double,double>::iterator i = data.begin(); i != data.end(); i++) {
      double x1 = i->first;
      double x2 = i->first;
      map<double,double>::iterator j = i;
      if (++j != data.end()) x2 = j->first;
      
      if (x1 != x2) {
	for (double x3 = x1; x3 < x2; x3 += (x2 - x1)*0.1) subsample[x3] = 0.0;
      } else subsample[x1] = 0.0;
    }

  int siz = data.size();
  double TotalZ = 0.0;
  
  for (int c = 0; c < MaxNr; c++) {
    Spline *probier = new Spline(siz,"gaga");
    int picked = 0;
    double pickfrac = 0.3;
    for (map<double,double>::iterator i = data.begin(); i != data.end(); i++) {
      
      double x = i->first;
      double y = i->second;
      
      bool pick = Miscmath::posRnd() < pickfrac;
      map<double,double>::iterator j = i;
      if (picked >= siz*pickfrac) pick = false;
      if (i==data.begin() || ++j==data.end()) pick = true;
     
      if (pick) { // pick it 
	picked ++;
	double stddev = max(range*0.01, fabs(y*0.05));
	double gauss = Miscmath::gaussRnd(y,stddev);
	probier->set(x,gauss);      
      }
    }
    
    probier->arm();

    double chi2 = 0;
    for (map<double,double>::iterator i = data.begin(); i != data.end(); i++) {
      double x = i->first;
      double y = i->second;
      double interpol = probier->fastY(x);
      double stddev = max(range*0.01, fabs(y*0.05));
      chi2 += (y-interpol)*(y-interpol) / (stddev*stddev);
    }
    chi2 /= siz;
    double like = exp(-0.5*chi2);
    for (map<double,double>::iterator i = subsample.begin(); i != subsample.end(); i++) {
      i->second += probier->fastY(i->first)*like;
    }
    TotalZ += like;
    //cout << "Spline: " << c << "  chi2: " << chi2 << "  likelihood: "<< exp(-0.5*chi2) << endl;
    delete probier;
  }

  for (map<double,double>::iterator i = subsample.begin(); i != subsample.end(); i++) {
    i->second /= TotalZ;
    perf << i->first << " "<< i->second;
    if (i->second > -1) perf << " " << exp(-0.5*i->second) << endl;
  }

  out = subsample;
}
*/  



/*
  Walk through .mcc MonteCarlo file and bin the models in steps bins from min_x to max_x where the x value is in column pos_x of the .mcc file. */
void ChainShop::getMCDistribution(string filename, map<double, vector<double> >& s,int pos_x, double min_x, double max_x, int steps,int *Progress, bool GridFile) {

  if (Progress) *Progress = 0;
  double step_x = (max_x - min_x)/steps; 
  vector<vector<float> > v;  readMCC(filename,v);

  int size = v.size();
  int collumns = v[0].size(); 

  map<int,float> store; 

  for (int i =0; i < steps; i++) store[i] =0;
  for (int line=0; line < size; line++) { 
    if (Progress) *Progress = (int)(100*line/(double)size);
    float x = v[line][pos_x];
    float multiple = v[line][collumns-1];

    int bin = (int)floor( (x - min_x)/step_x);
    store[bin] += multiple;
  }

  float highest = -1e10;

  for (int i =0; i < steps; i++) highest = max(highest,store[i]);

  vector<double> result(2);
  for (int i =0; i < steps; i++) {
    float a = min_x + i*step_x;
    result[0] = step_x;
    result[1] = store[i]/highest;
    if (isnan(result[1]) || isinf(result[1]))
        throw Bad_Error("isnan/isinf in ChainShop::getMCDistribution detected.\n"
                        "Possible cause: a column of all zeroes.");
    s[a] = result;
  }

  if (GridFile) {
    cout << "WRITING GRID FILE: " << endl;
    ofstream histo("gdh1.dat");
    for (int i =0; i < steps; i++) {
      float a = min_x + i*step_x;
      histo << a << "  " << s[a][1] << endl;
    }
  }

}

/*!
  In contrast to the usual getMCDistributionDerived, this one does not take the number count
  per bin, but the maximum likelihood that it finds within each bin. Use with care!
*/
void ChainShop::getMCDistributionDerived(string filename, map<double, vector<double> >& s,int pos_x, int pos_likeli, double min_x, double max_x, int steps,int *Progress, bool GridFile) {

  if (Progress) *Progress = 0;
  double step_x = (max_x - min_x)/steps;

  vector<vector<float> > v; readMCC(filename,v);
  int size = v.size();
  map<int,float> store; 

  

  for (int i =0; i < steps; i++) store[i] =-1e100;
  for (int line=0; line < size; line++) { 
    if (Progress) *Progress = (int)(100*line/(double)size);
    float x = v[line][pos_x];
    float loglikeli = v[line][pos_likeli];

    int bin = (int)floor( (x - min_x)/step_x);
    store[bin] =  max(store[bin],loglikeli);
    // cout << "LIKEPOS: " << pos_likeli << "  loglike: " << loglikeli << endl;
  }

  float highest = -1e10;  
  for (int i =0; i < steps; i++) highest = max(highest,store[i]);

  vector<double> result(2);
  for (int i =0; i < steps; i++) {
    float a = min_x + i*step_x;
    result[0] = step_x;
    cout << "storing: highest: " << highest << " " << store[i] << endl;
    result[1] = exp(store[i]-highest);
    s[a] = result;
  }
 
}




vector<double> ChainShop::peakAndSigma(map<double,double> &s) {
  map<double, int> trans;
  for (map<double,double>::iterator i = s.begin(); i != s.end(); i++) trans[i->first] = (int)rint(i->second);
  vector<double> nix;
  return nix;
}

/*! 
  Scan the montecarlo .mcc file for the best fit model
 */
vector<float> ChainShop::bestFitModel(string FileName, int likeli_pos) {
  vector<vector<float> > v; readMCC(FileName,v);
  if (likeli_pos>numberOfColumns(FileName)) {
    throw Bad_Error("ChainShop::bestFitModel() - likelipos > number of columns");
  }
  unsigned int lines = v.size();
  float bestloglike = -1e50;
  int bestrow = 0;
  for (unsigned int row = 0; row < lines; row ++) {
    if (v[row][likeli_pos] > bestloglike) {
      bestloglike = v[row][likeli_pos];
      bestrow = row;
    }
  }
  return v[bestrow];
}

/*! 
  Go through the .mcc file and return a structure that contains each point in the chain with it's value of 
  column 'column' and a corresponding relative value of 0..1 describing it's value with respect to 
  the minimum and maximum value of this column in the chain
*/
void ChainShop::getThirdDimension(string FileName, int x, int y, int column, ColorPoints &cp,int MaxPoints) {
  vector<vector<float> > v; readMCC(FileName,v);
						
  int size = v.size();
  int columns = v[0].size();

  if (column > columns || column < 0) throw Bad_Error("ChainShop::getThirdDimension() column out of range");
  
  float Min = v[0][column];
  float Max = v[0][column];
  for (int i = 1; i < size; i++) {
    if (v[i][column] < Min) Min = v[i][column];
    if (v[i][column] > Max) Max = v[i][column];
  }
  cp.Min = Min;
  cp.Max = Max;

  float Threshold = MaxPoints;
  Threshold /= size;

  cout << "THRESHOLD: " << Threshold << endl;

  list<vector<float> > tmp;
  for (int i = 0; i < size; i++) {
    if (Miscmath::posRnd(1.0) < Threshold) {
      vector<float> t(3);
      t[0] = v[i][x]; 
      t[1] = v[i][y]; 
      t[2] =  (v[i][column] - Min)/(Max - Min);
      tmp.push_back(t);     
    } 
  }

  size = tmp.size();
  cp.v.resize(size);
  list<vector<float> >::iterator j = tmp.begin();
  for (int i = 0; i < size; i++) {
    cp.v[i] = *j;
    j++;
  }
  cout << "Total ColorPoints: " << size << endl;
}



unsigned int ChainShop::numberOfColumns(const string& FileName)
{
  ifstream in(FileName.c_str());
  if (!in) throw Bad_Error("ChainShop::numberOfColumns() file not found: " + FileName);
  if (read<int>(in) != MAGICANALYZE) throw Bad_Error("ChainShop::numberOfColumns() Not a valid .mcc file: " + FileName);
  return read<int>(in); // reading number of columns
}


/*!
  Read in .mcc file and return vector of vector of float containing every line (and column)\
  of the .mcc file. The last double in each line is in fact the multiple which may have to
  be cast to integer sometimes
*/
void ChainShop::readMCC(string FileName,  vector< vector<float> >& v) {
  list<vector<float> > tmp;
  ifstream in(FileName.c_str());
  if (!in) throw Bad_Error("ChainShop::readMCC()  file not found");
  if (read<int>(in) != MAGICANALYZE) throw Bad_Error("ChainShop::readMCC() Not a valid .mcc file");
  int size = read<int>(in); // reading number of columns 
  vector<float>  keep(size);
  //  cout << "expected columns: " << size << endl;
  while (! in.eof()) {
    for (int i = 0; i < size-1; i++) keep[i] = read<float>(in);
    if (in.fail() ) break;
    int multiple = read<int>(in);    
    keep[size-1] = (float)multiple;
    tmp.push_back(keep);
  }

  int len = tmp.size();
  v.resize(len);
  int i =0;
  for (list<vector<float> >::iterator k = tmp.begin(); k != tmp.end(); k++,i++) {
    v[i] = *k;
  }
}


void ChainShop::revisitChain(string Original, string Revisited, int TotalLikePos, int likepos) {

  int TotalWeight = 30000;

  vector<vector<float> > v; readMCC(Original,v);
  int size = v.size();
  int columns = v[0].size();

  Miscmath::seed();


  float BestTotalLogLike =-1e100;
  float BestLogLike = -1e100;
  for (int i = 0; i < size; i++) {
    //    cout << "scanning: " << v[i][likepos] << "  " << v[i][TotalLikePos] << endl;
    BestLogLike = max(BestLogLike,v[i][likepos]);
    BestTotalLogLike = max(BestTotalLogLike,v[i][TotalLikePos]);
  }
  cout << "BestLogLike: " << BestLogLike << "  :: BestTotalLogLike: " << BestTotalLogLike << endl;
  
  map<int, int> note;

  int cnt =0;
  do { 
    int pos = (int)rint(Miscmath::posRnd()*(size-1));  // pick one line in the chain random 
    //cout << "pos: " << pos<< endl;
    float ratio = exp( v[pos][likepos] - BestLogLike); // relative propability w.r.t to best one
    
    //    ratio = 1.0;
    float ratio2 = exp( v[pos][TotalLikePos] - BestTotalLogLike) ;

    float ratio3 = ratio/ratio2;
    cout << "cnt: " << cnt << " pos: " << pos << " ratio: " << ratio << "  ratio2: " << ratio2 << " ratio3 " << ratio3 << endl;
  
    int add = (int)rint(ratio3);  
    add = min(100,add);
    add *= (int)rint(v[pos][columns-1]);
    cout << "adding: " << add << "  multipl: " << v[pos][columns-1] << endl;
    if (add > 0) {
      cnt++;
      if (note.find(pos) == note.end()) note[pos] = add; else note[pos]+=add;  // initialize and count updwards
    }
  } while (cnt < TotalWeight);

  int minmult = 1000*1000;
  for (map<int,int>::iterator i = note.begin(); i != note.end(); i++)  minmult = min(i->second,minmult);
  
  cout << "Minmult: " << minmult << endl;
  
  ofstream out(Revisited.c_str());
  write<int>(out,MAGICANALYZE);
  write<int>(out, columns);
  for (map<int,int>::iterator i = note.begin(); i != note.end(); i++) {
    int pos = i->first;
    int multiple = i->second / minmult;
    cout << "writing: "<< pos << "  :: " << multiple << endl;
    for (int k = 0; k < columns-1; k++) write<float>(out,v[pos][k]);
    write<int>(out, multiple);
  }
}


 // FITTING TO POLYNOMIAL
 // this procedure employs the gnu scientific library (GSL) and
 // makes heavy use of the example presented in the documentation concerning " Nonlinear least-squares fitting"
 // this sets up a Levenberg-Marquardt-solver and fits a sixth-order polynomial to the histogram. 
 //************************
/*!
  Fit to exponential polynomial.
  Levenberg-Marquardt-solver and fits a sixth-order polynomial to the histogram
*/

int ChainShop::OrderExpPoly(6);
vector<double> ChainShop::FitExpPoly_x_Values;

ExpPoly ChainShop::fitExpPoly(map<double, vector<double> >& binned, int order) {
  OrderExpPoly = order;
  if (OrderExpPoly > 6 || OrderExpPoly <2) throw Bad_Error("ChainShop::fitExpPoly() order out of bound");
  
  double fitParams[7]; // fitting coefficients will be stored here
  map<double,double> fitted;

  
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  size_t iter = 0;
	
  const size_t n = binned.size();
  const size_t p = order+1;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  vector<double> y(n), sigmaS(n);
  gsl_multifit_function_fdf f;

  double x_init[7] = {0.5,1.0,-1.0,0.0,0.0,0.0,0.0 };

  gsl_vector_view x = gsl_vector_view_array (x_init, p);

  const gsl_rng_type * type;
  gsl_rng * r;
  gsl_rng_env_setup();
  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  //  This is the data to be fitted 
  FitExpPoly_x_Values.resize(n);
  // We shift the x-axis, cause for very large x values , say for the fitting the loglikelihood
  // which may peak at x=-800 or so, we need to shift, cause exp(800) is not a good idea
  // numerically
  double shift = binned.begin()->first; // for simplicity, we shift at the very first value
  cout << "shifting: " << shift << endl;
  double invstretch = fabs(binned.rbegin()->first - shift);
  if ( invstretch == 0 )
    invstretch = 1.;
  invstretch = 1./invstretch;
  int j =0;
  for (map<double,vector<double> >::iterator k = binned.begin(); k != binned.end(); k++,j++) {
    y[j]= k->second[1];
    //cout << "looking at j=" << j << " y[j]=" << y[j] << endl;
    FitExpPoly_x_Values[j] = k->first + k->second[0]*0.5 - shift; // middle point
    cout << k->first  << " " <<  k->second[0]  <<  " " <<  shift << " " << invstretch << endl;
    FitExpPoly_x_Values[j] *= invstretch;
    cout << j << " - Entr " <<   FitExpPoly_x_Values[j]  << " at: " << y[j] << endl; 
    sigmaS[j]=0.2;
  }

  struct ExpFitPolyData d = { n, &y[0], &sigmaS[0]};
  f.f = &ChainShop::expb_f;
  f.df = &ChainShop::expb_df;
  f.fdf = &ChainShop::expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  // this is the loop for fitting
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate (s);

    // This is useful for debugging, uncomment if nesseccary
    cout << "status at iteration no. " << iter  << " is " << gsl_strerror (status) << endl;
    //  print_state (iter, s);

    // fit is done
    if (status) break;
    status = gsl_multifit_test_delta (s->dx, s->x, 1e-7, 1e-7);
  } while (status == GSL_CONTINUE && iter < 500);

  gsl_multifit_covar (s->J, 0.0, covar);  
  for(int l =0; l <= 6; l++) fitParams[l]=0;
  for (int l=0; l <= order; l++) fitParams[l]= gsl_vector_get(s->x,l);

  //This would give the parameters
  //cout << sqrt(gsl_matrix_get(covar,i,i));
  double err[7];
  for (int k=0; k <= order ; k++) {
    err[k]= sqrt(gsl_matrix_get(covar,k,k));
    if(!(err[k] == err[k]) || isinf(err[k])) {
      throw Bad_Error("ChainShop::fitExpPoly()  Unable to fit exponential");
    }
  }

  cout << "Fitting coefficients: " << endl; 
  for (int k=0; k <= order ; k++)  cout << "a" << k << " = " << fitParams[k] <<" +/- " << err[k] << endl;


  double chi = gsl_blas_dnrm2(s->f);
  printf("chisq/dof = %g\n",  pow(chi, 2.0)/ (n - p));
  cout << endl;
  
    
  // printf ("status = %s\n", gsl_strerror (status));
  gsl_multifit_fdfsolver_free (s);
  
  cout << "got fit " << endl;


  for (vector<double>::iterator i = FitExpPoly_x_Values.begin(); i != FitExpPoly_x_Values.end(); i++) {
    double x = *i;
    cout << "xing: "<< x << endl;
    double corresponds = x/invstretch +shift;
    cout << "corresponds to : " << corresponds << endl;
    fitted[corresponds] = exp(fitParams[0]+fitParams[1]*x+fitParams[2]*x*x+fitParams[3]*x*x*x  +
		    fitParams[4]*x*x*x*x+fitParams[5]*x*x*x*x*x
		    +fitParams[6]*x*x*x*x*x*x);
  }

  // initialize the communication structure (so it's for you to take it home)

  ExpPoly exppoly;

  exppoly.InvStretch = invstretch;
  exppoly.Shift = shift;
  exppoly.FitParams.resize(7);
  for (int i =0; i < 7; i++) exppoly.FitParams[i] = fitParams[i];
  return exppoly;
}


map<double,double> ChainShop::evaluateExpPoly(ExpPoly e, double x, double X,  int steps) {
  map<double,double> fitted;
  
  double Step = (X-x) / steps;
  double Stop = X + Step*1e-1;
  
  for (double y = x; y <= Stop; y+=Step) {
    double s = (y-e.Shift)*e.InvStretch;  
    cout << "evaluating at: "  << y <<  "  shifted etc is: " << s;
    fitted[y] = exp(e.FitParams[0]+e.FitParams[1]*s+e.FitParams[2]*s*s+e.FitParams[3]*s*s*s  +
		    e.FitParams[4]*s*s*s*s+e.FitParams[5]*s*s*s*s*s
		    +e.FitParams[6]*s*s*s*s*s*s);
    cout << "  result: "  << fitted[y] << endl;
    
  }
  return fitted;
}

int ChainShop::expb_f (const gsl_vector * x, void *params,gsl_vector * f)
{
  size_t n = ((struct ExpFitPolyData*)params)->n;
  double *y = ((struct ExpFitPolyData*)params)->y;
  
  double a3=0,a4=0,a5=0,a6=0;
  double a0= gsl_vector_get(x,0);
  double a1 = gsl_vector_get(x,1);
  double a2= gsl_vector_get(x,2);
  if(OrderExpPoly > 2) a3 = gsl_vector_get(x,3);
  if(OrderExpPoly > 3) a4= gsl_vector_get(x,4);
  if(OrderExpPoly > 4) a5 = gsl_vector_get(x,5); 
  if(OrderExpPoly > 5) a6 = gsl_vector_get(x,6);
  size_t i;
  
  for (i = 0; i < n; i++)
    {
      // Model Yi = Exp(a0+a1 x + a2 x^2 + ...)
      double t = FitExpPoly_x_Values[i];
      double Yi = exp(a0+a1*t+a2*t*t+a3*t*t*t+a4*t*t*t*t+a5*t*t*t*t*t+a6*t*t*t*t*t*t);
      gsl_vector_set (f, i,(Yi- y[i])/10); 
      //   cout <<  " Value: " << (Yi - y[i])  << " y[i] " << y[i] << " Yi " << Yi <<   endl;
    }
  return GSL_SUCCESS;
  
}

int ChainShop::expb_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  size_t n = ((struct ExpFitPolyData*)params)->n;
  double a0= gsl_vector_get(x,0);
  double a1 = gsl_vector_get(x,1);
  double a2= gsl_vector_get(x,2);
  double a3=0,a4=0,a5=0,a6=0;
  if(OrderExpPoly > 2) a3 = gsl_vector_get(x,3);
  if(OrderExpPoly > 3) a4= gsl_vector_get(x,4);
  if(OrderExpPoly > 4) a5 = gsl_vector_get(x,5); 
  if(OrderExpPoly > 5) a6 = gsl_vector_get(x,6);
  
  size_t i;
  
  for (i = 0; i < n; i++)
    {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/10,      */
      /*       Yi = exp(a0+a1*x+...)  */
      /* and the xj are the parameters (a0,a1,...) */
      double t =  FitExpPoly_x_Values[i];
      double e = exp((a0+a1*t+a2*t*t+a3*t*t*t+a4*t*t*t*t+a5*t*t*t*t*t+a6*t*t*t*t*t*t)); 
      gsl_matrix_set (J, i, 0, e/10.0); 
      gsl_matrix_set (J, i, 1, t*e/10.0); 
      gsl_matrix_set (J, i, 2, t*t*e/10.0);
      if (OrderExpPoly> 2) gsl_matrix_set (J, i, 3, t*t*t*e/10.0); 
      if (OrderExpPoly> 3) gsl_matrix_set (J, i, 4, t*t*t*t*e/10.0); 
      if (OrderExpPoly> 4) gsl_matrix_set (J, i, 5, t*t*t*t*t*e/10.0);
      if (OrderExpPoly> 5) gsl_matrix_set (J, i, 6, t*t*t*t*t*t*e/10.0);
    }
  return GSL_SUCCESS;
}

int ChainShop::expb_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J)
{
  expb_f (x, params, f);
  expb_df (x, params, J);
  return GSL_SUCCESS;
}

/*!
  Given a map of (x,y,chi2) values and a map which contains as input the sigma confidence and in 
  the end will be filled with a list of CoordinatePoints such that connecting the points will yield a contour
  
  \param s a map with (x,y,chi2) usually from get2dimMc...
  \param confidence a map where in the beginning, the list<CoordPoint> are all empty and for each sigma confidence
  that one is interested in there is one entry. E.g. confidence[1]=list<CoordPoint>(); confidence[2]=list<CoordPoint>();
  Fractional sigma's are perfectly ok.
*/

void ChainShop::getConfidenceLevels_DeltaChi2(map<float, map<float, float> > &s,map<float, list< CoordPoint >  >& confidence,int order,int quality) { 
  
 if (order > 4 || order <2) throw Bad_Error("Chainshop::getConfidenceLevels_DeltaChi2() order out of bound");
  int ExpansionParameters = 6;  // that's for 2nd order
  if (order == 3) ExpansionParameters = 10;
  if (order == 4) ExpansionParameters = 15;

  const int MAXCHI=30; // cutoff value
  float alpharange,betarange;  
  float alpha0=0, beta0=0,temp=1e100; // incredibly high chi2, this will not last ...

  //the following has been taken from the GSL example and modified
  // (sorry for the nasty C code)

  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *val, *w, *c;
  float x,y;
  //this is the parameter vector, 15 fit parameters
  c = gsl_vector_alloc (ExpansionParameters);

  cov = gsl_matrix_alloc (ExpansionParameters, ExpansionParameters);
  // find out the range of the parameters
  float min_X = s.begin()->first;
  float max_X = s.rbegin()->first;
  float min_Y = s.begin()->second.begin()->first;
  float max_Y = s.begin()->second.rbegin()->first;
  alpharange=max_X-min_X;
  betarange=max_Y - min_Y;

  // a counter 
  int points=0;
  float chi2=0;
  //find the minimum chi2 for normalizing later
  for (map<float, map<float,float> >::iterator i=s.begin(); i!=s.end();i++) {
    for (map<float,float>::iterator j=i->second.begin(); j!=i->second.end(); j++) {
      chi2=j->second;
      if(chi2 < temp){
	alpha0=i->first;
	beta0=j->first;
	temp = chi2;
      }
      if (chi2 < MAXCHI) points++;
    }
  }
  // allocate space
  X = gsl_matrix_alloc (points, ExpansionParameters);
  val = gsl_vector_alloc (points);
  w = gsl_vector_alloc (points);

  points=0;
   //need to find the coefficients of the fit, a polynomial of fourth order in x and y
  float alpha,beta;

  for (map<float, map<float,float> >::iterator i=s.begin(); i!=s.end();i++) {
    for (map<float,float>::iterator j=i->second.begin(); j!=i->second.end(); j++) {
      chi2=j->second;
      if(chi2 < MAXCHI){  
	alpha=i->first;
	beta=j->first; 
	x= (alpha-alpha0)/alpharange;  // we perform a transformation to different variables
	y= (beta-beta0)/betarange;    // this improves convergence behaviour
	chi2=chi2/MAXCHI;           // this is important when using the fit later
	
	gsl_matrix_set (X, points, 0, 1.0);
	gsl_matrix_set (X, points, 1, x);
	gsl_matrix_set (X, points, 2, y);
	gsl_matrix_set (X, points, 3, x*x);
	gsl_matrix_set (X, points, 4, x*y);
	gsl_matrix_set (X, points, 5, y*y);
	if (ExpansionParameters > 6) { 
	  gsl_matrix_set (X, points, 6, x*x*x);
	  gsl_matrix_set (X, points, 7, x*x*y);
	  gsl_matrix_set (X, points, 8, x*y*y);
	  gsl_matrix_set (X, points, 9, y*y*y);
	}
	if (ExpansionParameters > 10) {
	  gsl_matrix_set (X, points, 10, x*x*x*x);
	  gsl_matrix_set (X, points, 11, x*x*x*y);
	  gsl_matrix_set (X, points, 12, x*x*y*y);
	  gsl_matrix_set (X, points, 13, x*y*y*y);
	  gsl_matrix_set (X, points, 14, y*y*y*y);
	}
	gsl_vector_set (val, points, chi2);

	// all points are equally important (this sets the sigmas)
	gsl_vector_set (w, points, 0.1);
	points++;
      }
    }
  }
  gsl_multifit_linear_workspace * work 
    = gsl_multifit_linear_alloc (points, ExpansionParameters);
  gsl_multifit_wlinear (X, w, val, c, cov,&chisq, work);
  gsl_multifit_linear_free (work);
  
  // we would like to put the parameters into this vector for easy access
  vector<float> ps(15,0.0);
  for(int l=0; l < ExpansionParameters ; l++){
    ps[l]=gsl_vector_get(c,(l));
  }


  //
  //  find best chi2 by steepest decent.
  // 

  float x_best = 0; 
  float y_best = 0;
  float chi2_best = 1e10;

  bool improved;
  do  {
    x = x_best;
    y = y_best;
    float dchi2dx= MAXCHI*(ps[1] + 2*ps[3]*x+ps[4]*y + 3*ps[6]*x*x+
			   2*ps[7]*x*y+ ps[8]*y*y + 4*ps[10]*x*x*x+ 3*ps[11]*x*x*y+
			   2*ps[12]*x*y*y+ps[13]*y*y*y);
    float dchi2dy= MAXCHI*(ps[2] +ps[4]*x+ 2*ps[5]*y+
			   ps[7]*x*x+ 2*ps[8]*x*y+ 3*ps[9]*y*y+ps[11]*x*x*x+
			   2*ps[12]*x*x*y+3*ps[13]*x*y*y+ 4*ps[14]*y*y*y);
    
    float norm = sqrt(dchi2dx*dchi2dx + dchi2dy*dchi2dy);
    
    dchi2dx /= -norm;
    dchi2dy /= -norm;
   

    float old_chi2= MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
		  ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
		  ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y);
    
    if (chi2_best == 1e10) chi2_best = old_chi2;

    float best_inloop = 1e100;
    improved = false;
    for (float logstep = -20; logstep < -1; logstep += 0.5) {
      float step = pow((float)10.0,logstep);
      float x1 = x + dchi2dx * step;
      float y1 = y + dchi2dy * step;

      float X = x1*alpharange + alpha0;	
      float Y = y1*betarange + beta0; 
      bool inside=true;
      if (X < min_X || X > max_X) { inside = false;}
      if (Y < min_Y || Y > max_Y) { inside = false;}
      
      if (inside) {
	chi2= MAXCHI*(ps[0]+ps[1]*x1+ps[2]*y1+ps[3]*x1*x1+ps[4]*x1*y1+ps[5]*y1*y1+ps[6]*x1*x1*x1+
		  ps[7]*x1*x1*y1+ps[8]*x1*y1*y1+ ps[9]*y1*y1*y1+ps[10]*x1*x1*x1*x1+ps[11]*x1*x1*x1*y1+
		  ps[12]*x1*x1*y1*y1+ps[13]*x1*y1*y1*y1+ps[14]*y1*y1*y1*y1);
	if (chi2 < old_chi2 && chi2 < best_inloop  && (old_chi2 - chi2) < 2 && chi2 < chi2_best-1e-6) {
	  best_inloop = chi2;
	  chi2_best = chi2;
	  x_best = x1;
	  y_best = y1;
	  improved = true;
	}
      }
    }
  } while (improved);
  // little test
  /*
  float verybest = 1e20;
  for (int k = 0; k <100*100; k++) {
    
    x = x_best + Miscmath::rnd(1e-2);
    y = y_best + Miscmath::rnd(1e-2);

    chi2  = MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
	    ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
	    ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y);


    if (chi2 < chi2_best) cout << "FOUND BETTER: " << x << "  " << y << " :: " << chi2 <<  " still: verybest " << verybest << endl;
    
    if (chi2 < verybest) verybest = chi2;
  }
  */

  //
  // Uncomment the region below if you would like to see the input chi^2 (input) and the fit  to it (berg) 
  //
  /*
  ofstream berg("berg.dat"), input("input.dat");
  for (map<float, map<float,float> >::iterator i=s.begin(); i!=s.end();i++) {
    for (map<float,float>::iterator j=i->second.begin(); j!=i->second.end(); j++) {
      float X = i->first;
      float Y = j->first;
      x =(X-alpha0)/alpharange;
      y = (Y-beta0)/betarange;
      chi2  = MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
		      ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
		      ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y);
      berg << X << " " << Y << " " << chi2 << endl;
      input << X << " " << Y << " " << s[X][Y] << endl;
    }
  }
  */

  // we loop through for each value of sigma...
  Miscmath misc;
  ofstream onesigma("onesigma.dat"), twosigma("twosigma.dat");
  for (map<float, list< CoordPoint > >:: iterator conf = confidence.begin(); conf != confidence.end(); conf++) {
    float VALUE=misc.DeltaChi2CorrespondingToSigma(2.0,conf->first);   
    //    cout << "CORRESPONDING IS: " << conf->first << " -> " << VALUE << endl;
    for (float phi = 0; phi < 2*M_PI; phi += M_PI*0.015/quality) {
      float dx = cos(phi);
      float dy = sin(phi);
      //      cout << "Phi: " << phi << " dx: " << dx << "dy: " << dy << endl;
      float  r = 0;
      float dr = 1e-4;
      float X,Y;
      bool outof =false;
      for (;;) {
	x = x_best + r*dx;
	y = y_best + r*dy;	

	chi2= MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
		      ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
		      ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y)  - chi2_best;
	if (chi2 >  VALUE)  break;

	r+=dr;
	x = x_best + r*dx;
	y = y_best + r*dy;	
	float next_chi2 =  MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
			       ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
			       ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y)  - chi2_best;

	if (r > 1) {
	  outof = true;
	  break;
	}
	if (next_chi2 - chi2 > 0.2) dr *=0.5;
	if (next_chi2 - chi2 < 0.1) dr *= 1.2;
      }     
      if (! outof) {
	float dr = 0.5*r;
	//	cout << "initiating half step with: " << r << endl;
	for (int count=0; count < 12; count++) {
	  x = x_best + r*dx;
	  y = y_best + r*dy;
	  
	  chi2= MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
			ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
			ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y)  - chi2_best;
	  
	  //	  cout << "count: " << count << " x: "<< x << " y: " << y << " chi2: " << chi2 <<  "  r: " << r << endl;
	  if (fabs( (chi2-VALUE) / VALUE) < 1e-4) break;
	  if (chi2 > VALUE) r -= dr; else r += dr;
	  dr *= 0.5;
	}
      } else {
	//	cout << "so I didn't see a crossing" << endl;
	x = x_best + r*dx;
	y = y_best + r*dy;
      }     
      X = x*alpharange + alpha0;	
      Y = y*betarange + beta0;
      conf->second.push_back(CoordPoint(X,Y));
      if (conf->first == 1.0) onesigma << X << " " << Y << endl;
      if (conf->first == 2.0) twosigma << X << " " << Y << endl;
    }
  }
}
// END OF ORIGINAL VERSION


/*!
  Given a map of (x,y,chi2) values and a map which contains as input the sigma confidence and in 
  the end will be filled with a list of CoordinatePoints such that connecting the points will yield a contour
  
  \param s a map with (x,y,chi2) usually from get2dimMc...
  \param confidence a map where in the beginning, the list<CoordPoint> are all empty and for each sigma confidence
  that one is interested in there is one entry. E.g. confidence[1]=list<CoordPoint>(); confidence[2]=list<CoordPoint>();
  Fractional sigma's are perfectly ok.
*/

void ChainShop::getConfidenceLevels_Bayesian(map<float, map<float, float> >& s,
                                             map<float, list< CoordPoint > >& confidence,
                                             int order, int quality, int *Progress) {

  if (order > 4 || order <2) throw Bad_Error("Chainshop::getConfidenceLevels_Bayesian() order out of bound");
  int ExpansionParameters = 6;  // that's for 2nd order
  if (order == 3) ExpansionParameters = 10;
  if (order == 4) ExpansionParameters = 15;

  const int MAXCHI=30; // cutoff value
  float alpharange,betarange;
  float alpha0=0, beta0=0,temp=1e100; // incredibly high chi2, this will not last ...

  //the following has been taken from the GSL example and modified
  // (sorry for the nasty C code)

  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *val, *w, *c;
  float x,y;
  //this is the parameter vector, 15 fit parameters
  c = gsl_vector_alloc (ExpansionParameters);

  cov = gsl_matrix_alloc (ExpansionParameters, ExpansionParameters);
  // find out the range of the parameters
  float min_X = s.begin()->first;
  float max_X = s.rbegin()->first;
  float min_Y = s.begin()->second.begin()->first;
  float max_Y = s.begin()->second.rbegin()->first;
  alpharange=max_X-min_X;
  betarange=max_Y - min_Y;

  // a counter
  int points=0;
  float chi2=0;
  //find the minimum chi2 for normalizing later
  for (map<float, map<float,float> >::iterator i=s.begin(); i!=s.end();i++) {
    for (map<float,float>::iterator j=i->second.begin(); j!=i->second.end(); j++) {
      chi2=-2.*log(j->second);
      j->second=chi2;
      if(chi2 < temp){
        alpha0=i->first;
        beta0=j->first;
        temp = chi2;
      }
      if (chi2 < MAXCHI) points++;
    }
  }
  // allocate space
  X = gsl_matrix_alloc (points, ExpansionParameters);
  val = gsl_vector_alloc (points);
  w = gsl_vector_alloc (points);

  points=0;
   //need to find the coefficients of the fit, a polynomial of fourth order in x and y
  float alpha,beta;

  for (map<float, map<float,float> >::iterator i=s.begin(); i!=s.end();i++) {
    for (map<float,float>::iterator j=i->second.begin(); j!=i->second.end(); j++) {
      chi2=j->second;
      if(chi2 < MAXCHI){
        alpha=i->first;
        beta=j->first;
        x= (alpha-alpha0)/alpharange;  // we perform a transformation to different variables
        y= (beta-beta0)/betarange;    // this improves convergence behaviour
        chi2=chi2/MAXCHI;           // this is important when using the fit later

        gsl_matrix_set (X, points, 0, 1.0);
        gsl_matrix_set (X, points, 1, x);
        gsl_matrix_set (X, points, 2, y);
        gsl_matrix_set (X, points, 3, x*x);
        gsl_matrix_set (X, points, 4, x*y);
        gsl_matrix_set (X, points, 5, y*y);
        if (ExpansionParameters > 6) {
          gsl_matrix_set (X, points, 6, x*x*x);
          gsl_matrix_set (X, points, 7, x*x*y);
          gsl_matrix_set (X, points, 8, x*y*y);
          gsl_matrix_set (X, points, 9, y*y*y);
        }
        if (ExpansionParameters >10) {
          gsl_matrix_set (X, points, 10, x*x*x*x);
          gsl_matrix_set (X, points, 11, x*x*x*y);
          gsl_matrix_set (X, points, 12, x*x*y*y);
          gsl_matrix_set (X, points, 13, x*y*y*y);
          gsl_matrix_set (X, points, 14, y*y*y*y);
        }

        gsl_vector_set (val, points, chi2);

        // all points are equally important (this sets the sigmas)
        gsl_vector_set (w, points, 0.1);
        points++;
      }
    }
  }
  gsl_multifit_linear_workspace * work
    = gsl_multifit_linear_alloc (points, ExpansionParameters);
  gsl_multifit_wlinear (X, w, val, c, cov,&chisq, work);
  gsl_multifit_linear_free (work);

  // we would like to put the parameters into this vector for easy access
  vector<float> ps(15,0.0);  // we keep 15 and set the ExpansionParameters ... 14 range to zero!
  for(int l=0; l < ExpansionParameters ; l++){
    ps[l]=gsl_vector_get(c,(l));
  }

  // if the fit is really weird, then the minimal chi2 of the fit might not be in the vicinity
  // of the data's minimum chi2 so we roughly scan the entire fit for better values 
  // later, we will do a steepest decent from the our best estimate of the location
  // of the minimum to really find it
  float x_best = 0;
  float y_best = 0;
  float chi2_proposedbest =  MAXCHI*ps[0]; // we did expand around the proposed value
  cout << "proposed best: "<< chi2_proposedbest << endl;
  for (x = -1; x <= 1.001; x+=1e-2) {
    for (y = -1; y<=1.001; y+=1e-2) {
      float X = x*alpharange + alpha0;
      float Y = y*betarange + beta0;
      bool inside=true;
      if (X < min_X || X > max_X) { inside = false;}
      if (Y < min_Y || Y > max_Y) { inside = false;}
      if (inside) {
        chi2  = MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
                        ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
                        ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y);

        if (chi2 < chi2_proposedbest) {
          cout << "FOUND BETTER at: " << x << " , " << y << " with chi2: " << chi2 << " old: " << chi2_proposedbest << endl;
          x_best = x;
          y_best = y;
          chi2_proposedbest = chi2;
        }
      }
    }
  }

  //
  //  find best chi2 by steepest decent.
  //


  float chi2_best = 1e10;

  bool improved;
  do  {
    x = x_best;
    y = y_best;
    float dchi2dx= MAXCHI*(ps[1] + 2*ps[3]*x+ps[4]*y + 3*ps[6]*x*x+
                          2*ps[7]*x*y+ ps[8]*y*y + 4*ps[10]*x*x*x+ 3*ps[11]*x*x*y+
                          2*ps[12]*x*y*y+ps[13]*y*y*y);
    float dchi2dy= MAXCHI*(ps[2] +ps[4]*x+ 2*ps[5]*y+
                           ps[7]*x*x+ 2*ps[8]*x*y+ 3*ps[9]*y*y+ps[11]*x*x*x+
                           2*ps[12]*x*x*y+3*ps[13]*x*y*y+ 4*ps[14]*y*y*y);

    float norm = sqrt(dchi2dx*dchi2dx + dchi2dy*dchi2dy);

    dchi2dx /= -norm;
    dchi2dy /= -norm;


    float old_chi2= MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
                            ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
                            ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y);

    if (chi2_best == 1e10) chi2_best = old_chi2;

    float best_inloop = 1e100;
    improved = false;
    for (float logstep = -20; logstep < -1; logstep += 0.5) {
      float step = pow((float)10.0,logstep);
      float x1 = x + dchi2dx * step;
      float y1 = y + dchi2dy * step;

      float X = x1*alpharange + alpha0;
      float Y = y1*betarange + beta0;
      bool inside=true;
      if (X < min_X || X > max_X) { inside = false;}
      if (Y < min_Y || Y > max_Y) { inside = false;}

      if (inside) {
        chi2= MAXCHI*(ps[0]+ps[1]*x1+ps[2]*y1+ps[3]*x1*x1+ps[4]*x1*y1+ps[5]*y1*y1+ps[6]*x1*x1*x1+
                      ps[7]*x1*x1*y1+ps[8]*x1*y1*y1+ ps[9]*y1*y1*y1+ps[10]*x1*x1*x1*x1+ps[11]*x1*x1*x1*y1+
                      ps[12]*x1*x1*y1*y1+ps[13]*x1*y1*y1*y1+ps[14]*y1*y1*y1*y1);
        if (chi2 < old_chi2 && chi2 < best_inloop  && (old_chi2 - chi2) < 2 && chi2 < chi2_best-1e-6) {
          best_inloop = chi2;
          chi2_best = chi2;
          x_best = x1;
          y_best = y1;
          improved = true;
        }
      }
    }
  } while (improved);
  // little test
  /*
  float verybest = 1e20;
  for (int k = 0; k <100*100; k++) {

    x = x_best + Miscmath::rnd(1e-2);
    y = y_best + Miscmath::rnd(1e-2);

    chi2  = MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
                    ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
                    ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y);


    if (chi2 < chi2_best) cout << "FOUND BETTER: " << x << "  " << y << " :: " << chi2 <<  " still: verybest " << verybest << endl;

    if (chi2 < verybest) verybest = chi2;
  }
  */

  //
  // Uncomment the region below if you would like to see the input chi^2 (input) and the fit  to it (berg) 
  //

  /*
  ofstream berg("berg.dat"), input("input.dat");
  for (map<float, map<float,float> >::iterator i=s.begin(); i!=s.end();i++) {
    for (map<float,float>::iterator j=i->second.begin(); j!=i->second.end(); j++) {
      float X = i->first;
      float Y = j->first;
      x =(X-alpha0)/alpharange;
      y = (Y-beta0)/betarange;
      chi2  = MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
                      ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
                      ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y);
      berg << X << " " << Y << " " << chi2 << endl;
      input << X << " " << Y << " " << s[X][Y] << endl;
    }
  }
  cout << "best chi2: " << chi2_best << endl;
  */

  // first, get the total likelihood

  double total = 0;
  double dphi = 0.015*M_PI/quality; // advance the angle by dphi
  cout << "dphi: " << dphi << endl;
  for (float phi = 0; phi < 2*M_PI; phi += dphi) {
    float dx = cos(phi);
    float dy = sin(phi);
    double r =1e-5;  // first new r 
    double dr = 1e-5;  // stepsize
    double lastchi2 = -1e10;
    double last_likeli = 1.0; // last likelihood (we start at the best chi2=0, so likeli = 1.0).
    double last_r = 0; // r at which the last likeli was (we start at r =0);
    float X,Y; //,x,y;
    for (;;) {
      x = x_best + r*dx;
      y = y_best + r*dy;

      X = x*alpharange + alpha0;	
      Y = y*betarange + beta0;

      bool inside=true;
      if (X < min_X || X > max_X) { inside = false;}
      if (Y < min_Y || Y > max_Y) { inside = false;}
      if (!inside) break;

      chi2= MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
                    ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
                    ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y)  - chi2_best;

      if (chi2 < lastchi2 && false) {
        cout << "breaking: " << chi2 << " lastchi2: " << lastchi2 << endl;
        break;
      }
      lastchi2 = chi2;

      float likeli = exp(-0.5*chi2);
      //total += likeli * dr * r*dphi;
      float f1 = 0.5*last_likeli*last_r;
      float f2 = 0.5*likeli*r;
      total += (f1 + f2) * dr *dphi;
      //total += likeli * dr * r*dphi;

      //cout << "phi: " << phi << " r: " << r << " lastlike"  << last_likeli << " like: " << likeli << " dr: " << dr << endl;
      if (likeli/last_likeli > 0.97) dr *= 1.7; else dr *=0.5;  // 0.99

      if (likeli < 1e-5) break;  // 1e-5
      /*
      float c=0.1;
      double fprime = fabs(likeli - last_likeli)/dr;
      double f2fprime = 0.5*(last_likeli + likeli)/(fprime+1e-12);
      dr = min(c*f2fprime  , 3*dr);
      dr = max(1e-5,dr);
      */



      last_likeli = likeli;
      last_r = r;
      r += dr;

    }
    //    cout << "end of slice: " << phi << " at: " << r << endl;
  }
  //  cout << "TOTAL LIKELIHOOD IS: " << total << endl;
  //  throw Bad_Error("stop");

  // we loop through for each value of sigma...
  Miscmath misc;
  ofstream onesigma("onesigma.dat"), twosigma("twosigma.dat");
  // determine up to which area content we need it:

  float MaxRelativeArea =misc.fractionCorrespondingToSigma(confidence.rbegin()->first);   
  //  cout << "Max Relative Area: " << MaxRelativeArea << endl;

  // A map of splines that will store for each phi the radius(area)
  map<float,Spline*> Radius;
  map<float,float> R,R_new;
  Anchor anchor;
  for (float phi = 0; phi < 2*M_PI; phi += dphi) {
    R[phi] = 0; // start at 0
    Radius[phi] = new Spline(100,"radius",&anchor);
  }
  float delta_chi2=1e-3; // iso deltachi2
  double area=0;
  for (;;) {
    //
    // first, find the iso line that corresponds to some deltachi2
    //
    for (float phi = 0; phi < 2*M_PI; phi += dphi) {
      float dx = cos(phi);
      float dy = sin(phi);
      //X cout << "Phi: " << phi << " dx: " << dx << "dy: " << dy << endl;
      float  r = R[phi];
      float dr = 1e-6;
      //float X,Y;
      bool outof =false;
      for (;;) {
        x = x_best + r*dx;
        y = y_best + r*dy;

        chi2= MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
            ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
            ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y)  - chi2_best;
        if (chi2 >  delta_chi2)  break;

        r+=dr;
        x = x_best + r*dx;
        y = y_best + r*dy;
        float next_chi2 =  MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
            ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
            ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y)  - chi2_best;


        if (r>1.5) {
          outof = true;
          break;
        }
        if (next_chi2 - chi2 > 0.2) dr *=0.5;
        if (next_chi2 - chi2 < 0.1) dr *= 1.2;
      }
      if (! outof) {
        float dr = 0.5*r;
        //cout << "initiating half step with: " << r << endl;
        for (int count=0; count < 12; count++) {
          x = x_best + r*dx;
          y = y_best + r*dy;

          chi2= MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
              ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
              ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y)  - chi2_best;

          //   cout << "count: " << count << " x: "<< x << " y: " << y << " chi2: " << chi2 <<  "  r: " << r << endl;
          if (fabs( (chi2-delta_chi2) / delta_chi2) < 1e-4) break;
          if (chi2 > delta_chi2) r -= dr; else r += dr;
          dr *= 0.5;
        }
      } else {
        //X cout << "so I didn't see a crossing" << endl;
        // x = x_best + r*dx;
        //y = y_best + r*dy;
      }
      R_new[phi] = r;
      //X cout << "delta_chi2: " << delta_chi2 << " phi: " << phi << " R: " << R[phi] << " " << R_new[phi] << endl;
    }

    //
    // compute the area in the shell between R[phi] and R_new[phi]
    //
    //    cout << "olde area: " << area;
    for (float phi = 0; phi < 2*M_PI; phi += dphi) {
      float dx = cos(phi);
      float dy = sin(phi);
      float r1 = R[phi];
      float r2 = R_new[phi];
      int steps = 20;
      steps = 200;
      float dr = (r2-r1)/(steps);

      x = x_best + r1*dx;
      y = y_best + r1*dy;
      chi2 = MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
                     ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
                     ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y)  - chi2_best;
      float last_likeli = exp(-0.5*chi2);

      for (int i = 1; i <= steps; i ++) {
        float r =  r1 + i*dr;
        x = x_best + r*dx;
        y = y_best + r*dy;

        float X = x*alpharange + alpha0;
        float Y = y*betarange + beta0;
        bool inside=true;
        if (X < min_X || X > max_X) { inside = false;}
        if (Y < min_Y || Y > max_Y) { inside = false;}
        if (!inside) break;
        chi2= MAXCHI*(ps[0]+ps[1]*x+ps[2]*y+ps[3]*x*x+ps[4]*x*y+ps[5]*y*y+ps[6]*x*x*x+
                      ps[7]*x*x*y+ps[8]*x*y*y+ ps[9]*y*y*y+ps[10]*x*x*x*x+ps[11]*x*x*x*y+
                      ps[12]*x*x*y*y+ps[13]*x*y*y*y+ps[14]*y*y*y*y)  - chi2_best;

        float likeli = exp(-0.5*chi2);
        area += 0.5*(likeli*r + last_likeli*(r-dr) )*dr *dphi;
        last_likeli = likeli;
      }
      R[phi] = R_new[phi];
    }
    cout << "    new area: " << area << "  in percent: " << area/total << endl;
    // store the radius as a function of relative area for each angle
    for (float phi = 0; phi < 2*M_PI; phi += dphi) Radius[phi]->set(area/total,R[phi]);
    cout << "waiting for area to be bigger than max: " << MaxRelativeArea << endl;
    if (area/total > MaxRelativeArea) break;
    delta_chi2 += 1e-1;
  }
  for (float phi = 0; phi < 2*M_PI; phi += dphi) Radius[phi]->arm();

  //
  //  So now we have for each direction a Spline that given the relative integrated likelihood area/total
  //  yields the radius such
  //  all we have to do now is to infer for all requested sigmas the radius :-)
  //

  //  Radius[0]->dump("radius",false);

  for (map<float, list< CoordPoint > >:: iterator conf = confidence.begin(); conf != confidence.end(); conf++) {
    float RelativeArea =misc.fractionCorrespondingToSigma(conf->first);
    cout << "Looking for " << RelativeArea << endl;
    for (float phi = 0; phi < 2*M_PI; phi += dphi) {
      float dx = cos(phi);
      float dy = sin(phi);
      //      Radius[phi]->dump("radius",false);
      float r = Radius[phi]->fastY(RelativeArea);
      //      cout << RelativeArea << " phi: " << phi << " r: " << r << endl;
      x = x_best + r*dx;
      y = y_best + r*dy;
      float X = x*alpharange + alpha0;
      float Y = y*betarange + beta0;
      conf->second.push_back(CoordPoint(X,Y));
      if (conf->first == 1.0) onesigma << X << " " << Y << endl;
      if (conf->first == 2.0) twosigma << X << " " << Y << endl;
    }
  }
}



/*!
  Given the map (x,y,chi2) s, return ConfidenceRegions to plot ranging up to MaxSigma and having PerSigma steps per one additional sigma
  confidence. If Fill is true, allow color shading of the regions.
*/
list<ConfidenceRegion>* ChainShop::getConfidenceRegions(map<float, map<float, float> >& s,
                                                        float MaxSigma, int PerSigma, bool Fill,
                                                        ConfidenceInference infer,int order,
                                                        int quality, int *Progress) {
  int regions = 0;
  list<ConfidenceRegion>* C = new list<ConfidenceRegion>();
  map<float, list< CoordPoint > > confidence;
  float step = 1.0/PerSigma;
  for (double sig = step; sig <= MaxSigma+0.001; sig += step) {
    confidence[sig] = list<CoordPoint>();
    regions++;
  }

  switch (infer) {
  case DeltaChi2:
    getConfidenceLevels_DeltaChi2(s,confidence,order,quality);
    break;
  default:
    getConfidenceLevels_Bayesian(s,confidence,order,quality,Progress);
  }

  int k =regions-1;
  for (map<float, list<CoordPoint> >::reverse_iterator i=confidence.rbegin(); i != confidence.rend(); i++,k--) {
    ConfidenceRegion R;
    float sig = i->first;
    float  x = sig / 3.0;

    // Store color values for gui
    R.FillColor = shade(x);
    R.PenColor = R.FillColor;

    if (Fill) {
      Color OneSigma(1,1,1);
      if (MaxSigma < 2) OneSigma = Color(0,0,0);
      if (sig > 0.999 && sig < 1.001) { R.PenColor = OneSigma; }
      if (sig > 1.999 && sig < 2.001) { R.PenColor = Color(0,0,0); }
      if (sig > 2.999 && sig < 3.001) { R.PenColor = Color(0,0,0); }
    } else {
      R.PenColor.multiply(0.5);
    }
    R.region = i->second;
    R.Fill = Fill;
    C->push_back(R);
  }
  return C;
}


/*!
  Assymetric distribution for coloring the likelhood
*/
float ChainShop::asymmetric(float x, float peak, float peak_x, float up, float down) {
  float  y=0;
  if (x >= peak_x) {
    y=  peak + down*(peak_x-x);
  } else {
    y = peak + up*(x - peak_x);
  }
  if (y < 0) y = 0;
  if (y > 1) y = 1;
  return y;
}

/*!
  Heurisitc color shading scheme. Nasty but works
*/
Color ChainShop::shade(float x) {
  float green,blue,red;
    blue = asymmetric(x, 1.0 , 0.25, 1.5, 2.5);
    green =  asymmetric(x, 1 , 0.4 , 4, 3);
    red = asymmetric(x, 1, 0.7, 3,0);

    if (red > 1) red =1;
    if (blue > 1 ) blue = 1;
    if (green > 1) green = 1;
    green += asymmetric(x,0.6,0.5,5,1.5);

    if (red > 1) red =1;
    if (blue > 1 ) blue = 1;
    if (green > 1) green = 1;
    green += asymmetric(x,0.4,0.6,5,1);

    blue += asymmetric(x,0.98,1.0,2.3,2);
    green += asymmetric(x,0.93,1.0,2.2,2);

    if (red > 1) red =1;
    if (blue > 1 ) blue = 1;
    if (green > 1) green = 1;

    if (red < 0) red =0;
    if (blue < 0) blue = 0;
    if (green < 0) green = 0;

    return Color(red,green,blue);
}

float ChainShop::CrossCorrelation(string FileName, int x, int y) {
  vector< vector<float> > v;  readMCC(FileName,v);
  float cross = Correlation(v,x,y);
  float auto_x = Correlation(v,x,x);
  float auto_y = Correlation(v,y,y);

  cout << "cross: " << cross << "  auto_x: " << auto_x << "  auto_y: " << auto_y << endl;
  cout << "normalized: " <<  cross/sqrt(auto_x*auto_y) << endl;

  return cross/sqrt(auto_x*auto_y);
}

float ChainShop::Correlation(vector< vector<float> > &v, int x, int y)
{
  double sum_x=0;
  double sum_y=0;
  const double n = v.size();
  if (n<=0)
    throw Bad_Error("ChainShop::Correlation() - no entries in chain file.");
  const int multPos = v[0].size()-1;
  double totalMult = 0;
  for (unsigned int i = 0; i < n; ++i) {
    const double mult = v[i][multPos];
    totalMult += mult;
    sum_x += mult*v[i][x];
    sum_y += mult*v[i][y];
  }
  const double av_x=sum_x/totalMult;
  const double av_y=sum_y/totalMult;

  cout << "av_x: " << av_x << ", av_y: " << av_y << endl;

  double corr=0;
  for (unsigned int i = 0; i < n; ++i) {
    const double mult = v[i][multPos];
    corr +=mult*(v[i][x]-av_x)*(v[i][y]-av_y);
  }
  corr /= totalMult;
  return corr;
}

