#ifndef CHAINSHOP_H
#define CHAINSHOP_H
#include "../mathobject/mathobject.h"
#include <string>
#include <vector>
#include <list>
#include <map>
#include <fstream>
#include "lowlevelplot.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define tarrsiz 900
#define TE_ARRSIZE 512

using namespace std;
class Spline;
//class Model;
class Data;
class CL;
class Anchor;
class Cosmos;
class ControlPanel;

/*!
  Low level (internal) struct for fitting an exponential in
  a polynomial to the 1 dim marginalized data
*/
struct ExpFitPolyData {
  size_t n;
  double * y;
  double * sigmaS;
};

/*!
  High level (i.e. you can use it) struct for communicating
  and storing the exp(polynom) data to re-use etc...
  */
struct ExpPoly {
  int Order;
  float InvStretch;
  float Shift;
  vector<float> FitParams;
};

/*!
  Class for handling of Monte Carlo output. Mostly, the 
  routines are for "distilling" (i.e. converting to binary format and
  discarding the burn-in phase), marginalizing and likelihood
  distribution and confidence level extraction.

  These routines are used by the graphical user interface.
*/
class ChainShop
{
 private:
  float min_x, min_y, max_x,max_y; //!< used by rasterize; obsolete - use GridInfo
  float step_x,step_y; //!< used by rasterize and set by get2dim; obsolete - use GridInfo

  static vector<double> FitExpPoly_x_Values;  //!< For the 1-dim exp(polynom) fit to distribution
  static int OrderExpPoly; //!< Order of ExpPoly fitting
 public:
  ChainShop();
  ~ChainShop();

  enum ConfidenceInference { Bayesian, DeltaChi2 };
  
  /* ***************************************************************************************
  **
  **  The following functions are for marginalization, likelihood contour finding,
  **  plotting, rasterizing etc. 
  **
  ***************************************************************************************** */

  //! read binary variable from instream
  template<class T> static T read(istream& i) {
    T s;
    i.read((char*) &s,sizeof(T));
    return s;
  }

  //! write binary variable to outstream
  template<class T> static void write(ostream& o,const T& s) {
    o.write((const char*) &s,sizeof(T));
  }

  //! convert several chain-output data files to one binary .mcc file
  static double distillChain(vector<string> infiles, string outfile, int* Progress=0, int discard=0,bool OnlyBinary=false,bool OldVersion=false);
  //! Read .mmc file and return each line as a float vector
  static void readMCC(string FileName, vector< vector<float> > &v);

  //! Return the number of columns in a the MCC file FileName
  static unsigned int numberOfColumns(const string& FileName);

  static void revisitChain(string Original, string Revisited, int totallikepos, int likepos);

  static void likeli2chi2(map<double,double>&);  //!< convenience function
  static void chi2likeli(map<double,double>&); //!< convenience function

  static vector<float> bestFitModel(string FileName, int likeli_pos);

  //  void perfectSpline(map<double,double> &, map<double,double>& s,int MaxNr=1000); 

  //! Scan .mcc file, and build map of x vs. likelihood of 1-dim marginalized likelihood of column pos_x
  /*
  void getMCDistribution(string filename, map<double,double>&,int pos_x, int pos_likeli,double smear, double min_x, double max_x, int steps,int *Progress=0,bool GridFile=true);
  */

  //!< Scan .mcc file, and build map of x vs. likelihood of 1-dim marginalized likelihood of column pos_x. Uses binned number count
  void getMCDistribution(string filename, map<double, vector<double> >&,int pos_x, double min_x, double max_x, int steps,int *Progress=0,bool GridFile=true);
   //!< Scan .mcc file, and build map of x vs. likelihood of 1-dim marginalized likelihood of column pos_x. Uses  maximum likelihood in bin (use with care)
  void getMCDistributionDerived(string filename, map<double, vector<double> >&,int pos_x, int pos_likeli, double min_x, double max_x, int steps,int *Progress=0,bool GridFile=true);

  
  
  static int expb_f (const gsl_vector * x, void *params,gsl_vector * f); //!< 
  static int expb_df (const gsl_vector * x, void *params, gsl_matrix * J);
  static int expb_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);
  static ExpPoly fitExpPoly(map<double, vector<double> > & binned, int order=6);
  
  //! Use the exp(poly) fit ExpPoly e and return a map of points (x,likelihood) with x=x...X and steps number of steps
  static map<double,double> evaluateExpPoly(ExpPoly e, double x, double X,  int steps); 

  //! small struct to pass information about the grid for rastered 2-d likelihood plots around
  struct GridInfo
  {
    float Min_x, min_x, Min_y, min_y, Max_x, max_x, Max_y, max_y, step_x, step_y, bins;
    std::pair<float, float> minCoords;

    float bin2x(unsigned int i) const { return min_x + i*step_x; }
    float bin2y(unsigned int j) const { return min_y + j*step_y; }
    unsigned int x2bin(float x) const { return (unsigned int)((x-min_x)/step_x); }
    unsigned int y2bin(float y) const { return (unsigned int)((y-min_y)/step_y); }
  };

  typedef map<float, map<float, float> > Grid;

  enum DistributionType { UseChi2, UseWeight };
  enum AveragingType { UseMaximum, UseAverage, UseTotal };
// debugging
#ifndef PRERELEASE
  void setDistributionType(DistributionType t) {mDistributionType=t;}
  void setAveragingType(AveragingType t) { mAveragingType=t; }
  private:
  DistributionType mDistributionType;
  AveragingType    mAveragingType;
  public:
#endif

  //! Scan .mcc file and build map of (x,y) vs likelihood of 2-dim marginalized likelihood of columns pos_x,pos_y  
  void get2DimMCDistribution(string filename, Grid& s,int pos_x, int pos_y, int pos_likeli,
                             GridInfo& gridInfo, int bins=500, float SmearFactor=3,int* Progress=0,bool GridFile=true);
  //! Given a (x,y) vs likelihood map, return a line-wise rastered image. 
  list<Block>* rasterize(map<float, map<float, float> >& ,int* =0);
  list<Block>* rasterize_Bayesian2(Grid& , const GridInfo&, int* =0);
  //! Little helper for color shading 
  static float asymmetric(float x, float peak, float peak_x, float up, float down);    
  //! Core routine to find 2-d confidence levels given a map with (x,y,likelihood)
  static void getConfidenceLevels_DeltaChi2(Grid &s, map<float, list<CoordPoint > >&,int order=4,int quality=1);
  static void getConfidenceLevels_Bayesian(Grid &s, map<float, list<CoordPoint > >&,int order=4,int quality=1,int *Progress=0);
  //! High level confidence region routine. Will return Confidence regions nicely shaded
  static list<ConfidenceRegion>* getConfidenceRegions(Grid& s, float MaxSigma=3.0, int PerSigma = 1,bool Fill=true,ConfidenceInference=Bayesian,int order=4,int quality=1,int *Progress=0);
    //! Little shading helper
  static Color shade(float x); //!< shade from x=0 -> dark blue to x=1 -> light red
  //! return third dimension information for some third dimension column, allows to color encode third dimension in 2-d plot
  static void getThirdDimension( string FileName, int x, int y, int column, ColorPoints& result,int MaxPoints=500);

   
  static vector<double> peakAndSigma(map<double,double>&);  // return peak and +- sigma regions
   
   //! Write 1-dim distribution to file
   static void dumpDistribution(string filename, map<double,double>&);
   
  //! Return parameter range covered by parameter in column pos_x and pos_y
   static vector<float> autoScale(string filename,int pos_x, int pos_y);
  
   static float Correlation(vector< vector<float> > &v, int x, int y);
   static float CrossCorrelation(string filename, int x, int y);


};

#endif
