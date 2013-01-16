#ifndef ANALYZETHIS_H
#define ANALYZETHIS_H
#include "../mathobject/mathobject.h"
#include "../mathobject/snedata.h"
#include "newdatreader.h"
#include <string>
#include <vector>
#include <list>
#include <map>
#include <fstream>
//#include "lowlevelplot.h"
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
//#include "chainshop.h"

#define tarrsiz 900
#define TE_ARRSIZE 512

using namespace std;
class Spline;
//class Model;
class Data;
class CL;
class Anchor;
class baseCosmos;
class Cosmos;
class ControlPanel;

/*!
  Small struct controlling some parameters for postscript output of
  2-dimensional marginalized likelihood plots

struct Printing  {
  float xtick, ytick;
  float StartTick_x, StartTick_y;
  string XLabel, YLabel;
  int labelsize, ticksize;
  float xlabeloffset,ylabeloffset;
  float scaleX,scaleY; // paper scaling to 0..1, 0..1
  
  LowLevelPlot::LabelStyle XTickLabelStyle, YTickLabelStyle;

  // cmbeasy gui stuff
  bool ProtectedSize;
  bool AutomaticTick;
  Printing(int ls=32,int ts=32,float xo=0.05, float yo=0.05, float sx=593.3, float sy=841.9) : labelsize(ls), ticksize(ts), xlabeloffset(xo), ylabeloffset(yo),scaleX(sx), scaleY(sy) , XTickLabelStyle(LowLevelPlot::decimal), YTickLabelStyle(LowLevelPlot::decimal), ProtectedSize(true) , AutomaticTick(true) {};
  };*/

/*!
  Small struct used to return a normalization factor and the chi^2 for WMAP
*/
struct WMAPNorm {
  double norm;
  double chi2_tt, chi2_te;
};

/*
  Small struct (with a lot of content) that holds the likelihood regions in a rasterized
  format. For as long as you do not plan to enhance the graphical user interface, there
  is no need to try and understand my (potentially) weird implementation of this task :-)

struct RasterizeReturn {
  RasterizeReturn() : lines(0) {};
  ~RasterizeReturn() { if (lines) delete lines; }
  int regions;
  float min_x,min_y,max_y,range_x,range_y,step_x;
  pair<float,float> best;
  vector<bool> border;
  vector< vector< float > > color;
  vector< map<float, list< pair<float, float> > > > *lines;
};
*/
/*!
  Low level (internal) struct for fitting an exponential in
  a polynomial to the 1 dim marginalized data

struct ExpFitPolyData {
  size_t n;
  double * y;
  double * sigmaS;
};


  High level (i.e. you can use it) struct for communicting
   and storing the exp(polynom) data to re-use etc... 

struct ExpPoly {
  float InvStretch;
  float Shift;
  vector<float> FitParams;
  }; */

/*!
  Main data analysis class of cmbeasy. Provides comparison of
  numerical results to CMB, SNe Ia and Large Scale Structure measurements.

*/
class AnalyzeThis : public Mathobject {
 private:
  float range_x, range_y; //!< used by rasterize 
  float min_x, min_y, max_x,max_y; //!< used by rasterize
  float scaleX, scaleY; //!< used by rasterize
  float step_x,step_y;
  bool printit;
  vector < vector<float> > mz;
  vector < vector<float> > filter;
  vector <int> dimensionality;
  vector< map<float, int> > values;
  map <long , float* > HashMap;

  //SNe Data sets
  //Knop03
  SNeData Knop03Extended, Knop03ExtendedLowe, Knop03ExtendedLoweStrict;
  SNeData Knop03HST, Knop03HSTLowe;
  //Riess04
  SNeData Riess04all, Riess04gold;
  //Tonry03
  SNeData Tonry03full, Tonry03_172, Tonry03_130;
  //Barris04
  SNeData Barris04;
  //Astier05
  SNeData Astier05;

  Data *BinnedWMapTT, *BinnedWMapTE,*DiagMapTE;
  Data *QSO;
  Spline *tmpSig8;

  double cl_data[tarrsiz+1], neff[tarrsiz+1], fskyeff[tarrsiz+1];
  double te_data[TE_ARRSIZE+1], cltt_data[TE_ARRSIZE+1], ntt[TE_ARRSIZE+1], nee[TE_ARRSIZE+1];
  vector< vector< double> > off_diag, r_off_diag;
  vector< vector< double> > te_off_diag;
  
  double sig8Integrator(const double k) const; //!< used by sigma8()


  //Routines for reading in SN data Conley
  bool knop03read; //!< Have the Knop03 SNe been read in?
  bool riess04read; //!< Have the Riess04 SNe been read in?
  bool tonry03read; //!< Have the Tonry03 SNe been read in?
  bool barris04read; //!< Have the Barris04 SNe been read in?
  bool astier05read; //!< Have the Astier05 SNe been read in
  int readKnop03Data(); //!< Reads in Knop03 data
  int readRiess04Data(); //!< Reads in Riess04 data
  int readTonry03Data(); //!< Reads in Tonry03 data
  int readBarris04Data(); //!< Reads in Barris04 data
  int readAstier05Data(); //!< Reads in Astier05 data

  //Old routines for sne Ia
  Data *Riess;
  Data *Riess06_gold, *Riess04_gold, *Riess04_all;
  Data *Tonry, *Tonry_subsample172;
  Data *Tonry_subsample130, *IfADS ;
  Data *HSTSn1a, *HSTLowZ;
 

 public:
  AnalyzeThis();
  ~AnalyzeThis();
  
  enum Tonry03Sample { TonryFull, Tonry172, Tonry130 };
  enum Knop03Sample { K03Primary, K03Lowe, K03Lowestrict, K03HST, K03HSTLowe };

  double sigma8(Spline* s,double h); //!< given powerspectrum s, return sigma_8
    
  /* WMAP  1yr */
  double WMAP_TT(Spline& ts); //!< return chi^2 w.r.t WMAP TT
  double WMAPDiag_TT(Spline &ts); //!< return chi^2 w.r.t only the diagonal part of the covariance matrix of WMAP TT
  double WMAP_TE(Spline& ts, Spline &cs, Spline &es); //!< return chi^2 w.r.t WMAP TE
  double WMAPDiag_TE(Spline &ts, Spline &cs, Spline &es); //!< return chi^2 w.r.t only the diagonal part of the covariance matrix of WMAP TE

  
  void scaleCls(CL&, int n, double normalize); //!< multiply spectra by a factor
  void scaleScalarCls(CL&, int n, double normalize); //!< multiply scalar spectra by a factor
  void scaleTensorCls(CL&, int n, double normalize); //!< multiply tensor spectra by a factor

  void quickWMAPNormalize(Cosmos& cosmos,const ControlPanel& control,CL& cl); //!< Normalize to WMAP diagonal TT data. This is fast and a replacement for the old cobenormalize(). However, cl's must be in muK^2
  
  bool isSane(CL&,int,double*,double*);
  vector<WMAPNorm> WMAPNormalize(Cosmos &cosmos,ControlPanel& control, CL &cl,bool fast=false); //!< Get the best fit to the total WMAP likelihood. Implies quickWMAPNormalize(). Rather slow. Cl's must be in muK^2 

  pair<double,double> bestFittingBinnedWMAP_TT(Spline& ts); //!< return the normalization factor of the Temperature spectrum ts w.r.t the binned WMAP TT data
  pair<double,double> bestFitting(Spline& ts, Data*, int ScaleDof,bool FindFactor=true); //!< compute the best fitting normalization factor w.r.t to some data 
  

  void initWMAPCommon(int *Progress=0); //!< wrapper, calls initMapCommon for TT and TE 
  void initWMAPCommonTT(int *Progress=0); //!< read in WMAP TT covariance matrix 
  void initWMAPCommonTE(int *Progress=0); //!< read in WMAP TE covariance matrix
  bool WMAPNotYetInitialized; //!<  read in covariance of WMAP yet ?

  /**********************************************************
   *
   *                     WMAP 7-year likelihood
   *
   **********************************************************/
  struct WMAP7Likelihood
  {
    double TTlike;           //!<  master tttt
    double TTlowllike;       //!<  low tttt chisq
    double TTlowldet;        //!<  low tttt determinant
    double Beamlike;         //!<  beam/pt source correction to tttt
    double TElike;           //!<  master tete
    double TEdet;            //!<  master tete determinant
    double Lowllike;         //!<  TE/EE/BB lowl
    double Lowldet;          //!<  lowl det
    double TBlike;           //!<  master tbtb chisq flag
    double TBdet;            //!<  master tbtb determinant flag
    double total;            //!<  total -2Ln(L)
  };

#ifdef WMAPTHREEFORTRAN   //when using the Fortran90 code directly
  void WMAP3F90SetMaxTT( int l ); //!< for the options, see the README of the WMAP distribution; if these routines are not called, the defaults in the ..._options.f90 file are used
  void WMAP3F90SetMinTT( int l );
  void WMAP3F90SetMaxTE( int l );
  void WMAP3F90SetMinTE( int l );
  void WMAP3F90SetLowlMax( int l );
  void WMAP3F90SetUseLowlTT( bool b );
  void WMAP3F90SetUseLowlPol( bool b );
  void WMAP3F90SetUseTT( bool b );
  void WMAP3F90SetUseTE( bool b );
  void WMAP3F90SetSZAmplitude( double amp );
  bool WMAP3F90Error();
  void WMAP3F90ErrorReport();


  WMAP7Likelihood WMAP3F90computeLikelihood( CL& cl );
#endif

#ifndef NOWMAP7
  /*! \brief compute likelihood with respect to WMAP 7-year data
   * 
   *  returns a WMAP7Likelihood struct that contains the likelihoods
   *  for the given Cls:
   *  
   *  \code
   *  WMAP_OPTIONS::self()->setSz_amp(  sz_amp );
   *  AnalyzeThis::WMAP5Likelihood like = ai.WMAP7computeLikelihood( cl );
   *  cout << "Total likelihood: " << like.total;
   *  task[ LOGLIKEPOS+2 ] = -0.5*like.TTlike;
   *  \endcode
   */

  WMAP7Likelihood WMAP7computeLikelihood( CL& cl ); //!< return WMAP7 likelihood for the given cl; for options see the file WMAP_3yr_options.h

#endif

 private:
  void ConvertWMAP2BinaryTT(int *Progress=0); //!< convert WAMP TT data files to binary for faster startup
  void ConvertWMAP2BinaryTE(int *Progress=0); //!< convert WMAP TE data files to binary for faster startup
 public:
  /*! Apply the initial amplitudes A_s (for scalars) and A_t (for tensors) to the cl's and sigma8 etc.;
   *  when force=true, don't warn if normalization has already been set
   */
  void rescaleSpectra(Cosmos& cosmos,const ControlPanel& control,CL& cl, vector<double>& A_s, vector<double>& A_t, bool force=false);
  //! same as above, applies the same A_s and A_t to all spectrac cl::xx[i]
  void rescaleSpectra(Cosmos& cosmos, const ControlPanel& control, CL& cl, double A_s, double A_t, bool force=false);
  vector<WMAPNorm> WMAPLikelihood(Cosmos& cosmos,ControlPanel& control,CL &cl);
  void  fiducialAmplitudes(Cosmos& cosmos,vector<double> &A_s,vector<double> &A_t);
  void applyInflationaryTensorRatio(Cosmos &cosmos,vector<double> &A_s, vector<double> &A_t);

  /*
  **
  ** Supernovae stuff author: A. Conley
  **
  */
  double SNIaCore(const baseCosmos&, const SNeData&, double,
		  double, double, double, double,
		  bool, double, double) const; //!< SN core chi^2 routine (A. Conley)
  double estimateScriptM(const baseCosmos& cosmos,const SNeData& data,
			 double alpha, double beta, double intrinsicdisp, 
			 double pecz, double widthmean) const; //!< Estimate scirptM given other params  (A. Conley)

  double SNIaKnop03(const baseCosmos& cosmos, double alpha,
		    Knop03Sample sample = K03Primary, double Rb = 4.1,
		    double intrinsicdisp = 0.11, double pecz = 0.001);//!< chi^2 routine for Knop 03 data  (A. Conley)
  double estimateScriptMKnop03(const baseCosmos& cosmos, double alpha,
			       Knop03Sample sample = K03Primary, double Rb=4.1,
			       double intrinsicdisp=0.11, 
			       double pecz=0.001); //!< Estimate scriptM for Knop03 sample (A. Conley)

  double SNIaRiess04(const baseCosmos& cosmos, bool goldsample=true); //!< chi^2 routine for Riess04 data  (A. Conley)
  double SNIaRiess06(const baseCosmos& cosmos, bool goldsample=true); //!< chi^2 routine for Riess06 data  (A. Conley)
  double estimateScriptMRiess04(const baseCosmos& cosmos, bool goldsample=true); //!< Estimate scriptM for Riess04 sample  (A. Conley)

  double SNIaTonry03(const baseCosmos& cosmos, Tonry03Sample sample = Tonry172); //!< chi^2 routine for Tonry03 data  (A. Conley)
  double estimateScriptMTonry03(const baseCosmos& cosmos, Tonry03Sample sample = Tonry172); //!< Estimate scriptM for Tonry03 sample

  double SNIaBarris04(const baseCosmos& cosmos); //!< chi^2 routine for Barris04 data  (A. Conley)
  double estimateScriptMBarris04(const baseCosmos& cosmos); //!< Estimate scriptM for Barris04 sample  (A. Conley)

  double SNIaAstier05(const baseCosmos& cosmos, double alpha, double beta,
		      double intrinsicdisp=0.131, double pecz=0.001); //!< chi^2 routine for Astier05 data  (A. Conley)
  double estimateScriptMAstier05(const baseCosmos& cosmos, double alpha, 
				 double beta, double intrinsicdisp=0.131,
				 double pecz=0.001); //!< Estimate scriptM for Astier05 sample  (A. Conley)


  //
  // From here on, the old supernovae routines
  //
  double Sn1aCore(Spline& luminositydistance, Data &,double=5.0); //!< **Legacy stuff.  Better use routines by A. Conley  ** [Core routine for Tonry and HST]
  double Sn1aRiess04(Spline& luminosityDistance, bool goldsample=true); //!< **Legacy stuff. Better use routines by A. Conley  **  [chi^2 w.r.t Riess04 data]
  double Sn1aRiess06(Spline& luminosityDistance, bool goldsample=true); //!< **Legacy stuff. Better use routines by A. Conley  **  [chi^2 w.r.t Riess04 data]
  double Sn1aRiess( Spline& luminosityDistance); //!< **Legacy stuff. Better use routines by A. Conley  **  [chi^2 w.r.t to Riess data]
  double Sn1aHST(Spline &luminosityDistance,bool subsample=true); //!< **Legacy stuff. Better use routines by A. Conley  **  [chi^2 w.r.t HST data]

  /*
  **
  ** BOOMERANG from F.PIACENTINI (astro-ph/0507507, TE), T.E. Montroy (astro-ph/0507514, EE & BB), W.C. Jones (astro-ph/0507494,TT)
  **
  */
 public: 
  map<int,bool> BOOMERANGBandDeselect; //!< All integers appearing here will be considered de-selected from the BOOMERANG likelihood. Counting goes from 1 to 47
  /*! Deselect entire spectra with this function
    \param TT if true, deselect all TT points
    \param TE if true, deselect all TE points
    \param EE if true, deselect all EE points
    \param BB if true, deselect all BB points
  */
  void BOOMERANGDeselect(bool TT, bool TE, bool EE, bool BB); //!< deselect entire spectra
  void BOOMERANGDeselect(int band); //!< deselect individual bands
  void BOOMERANGInit( CL &cl ); // all C'ls in muK^2
  //! return chi2 w.r.t BOOMERANG05, call BOOMERANGInit first
  double BOOMERANGChi2(double beamErrorAdjust=0, double calibrationFactor=1);
  //! return marginalized chi2 w.r.t BOOMERANG05, call BOOMERANGInit first
  double BOOMERANGChi2WithCalibration(); 
 private:
  double BOOMERANGWindowConv(int windowNumber, Spline &ClTT, Spline &ClTE, Spline &ClEE, Spline &ClBB);
  double BOOMERANGthband[47];
  vector<NewdatReader::Band> BOOMERANGexpData;
  NewdatReader::Matrix BOOMERANGFisherMatrix;
    /*
  **
  ** VSA
  **
  */
 public: 
  map<int,bool> VSABandDeselect; //!< All integers appearing here will be considered de-selected from the VSA likelihood. Counting goes from 1 to 16.
  double VSAChi2(Spline &Cl); //!< return chi2 w.r.t VSA2004 all C'ls in muK^2
  double VSAChi2WithCalibration(Spline &Cl);   //!< all Cl's in muK^2
 private:
  double VSAWindowConv(int windowNumber,Spline &Cl, double cut_l);
    /*
  **
  **CBI 2000+2001 (from A. C. S. Readhead et al. (2004, astro-ph/0402359).
  ** very much like the VSA computation
  **
  */
 public: map<int,bool> cbiBandDeselect; //!< All integers appearing here will be considered de-selected from the CBI2000+01 likelihood. Counting goes from 1 to 13.

  double CBIChi2(Spline &Cl); //!< return chi2 w.r.t CBI200+2001. Cl's in muK^2
  double CBIChi2WithCalibration(Spline &Cl); //!< Cl's in muK^2 
 private:
  double CBIWindowConv(int windowNumber,Spline &Cl, double cut_l);

  /*
  **
  ** CBI and ACBAR
  **
  */
 public:
  map<int,bool> cbiMosaicBandDeselect; //!< All integers appearing here will be considered de-selected from the CBI Mosaic likelihood. Counting goes from 1 to 14.
  map<int,bool> acbarBandDeselect; //!< All integers appearing here will be considered de-selected from the ACBAR likelihood. Counting goes from 1 to 14.
 
  double CBIDeepChi2(Spline &); //!< return chi2 w.r.t CBI-DEEP
  double CBIMosaicChi2(Spline &); //!< return chi2 w.r.t  CBI-MOSAIC 
  double ACBARChi2(Spline &); //!< return chi2 w.r.t ACBAR
  double CBIMosaicChi2WithCalibration(Spline&); 
  double ACBARChi2WithCalibration(Spline&);
 private:
  double WindowConv(string windowFileName, Spline &Cl, int  binnumber, int binwidth, double cut_l);
  double WindowConvACBAR(string windowFileName, Spline &Cl, double cut_l);
 public:
  /*
  **
  ** SDSS
  **
  */
  double SDSS_chiSquared(Spline *powerSpline, double bias=1.0, double Kbreak=0.15);
  double SDSS_bestBias(Spline *powerSpline, double Kbreak=0.15);
  double SDSS_Pk(float k, Spline *powerSpline);
  
  /*
  **
  ** 2dFGRS
  **
  */
  double TwoDF_bestBias(Spline*);
  double TwoDF_convolutePowerSpectrum( Spline*);
  float TwoDF_Pk(float k, Spline *spline); 
  float TwoDF_PkconvW(float rk, Spline *powerSpline);
  float TwoDF_Wksq(float k);
  double TwoDF_chiSquared(vector<double> convTheoreticalData, vector<double> k_theor);

  /*
  **
  ** Baryon Acoustic Peak
  **
  **/
  double SDSS_BAP_A_chiSquared(baseCosmos&, double ns, double Ode=0) const; //!< Returns chi^2 based on the A parameter of Eistenstein '05

  double Angular_BAP_A_chiSquared(baseCosmos&, double ns, double Ode=0) const; //!< Returns chi^2 based on the angular distance ratio of several BAO datasets

  /*
  **
  ** Varying fine structure constant
  **
  */
  double qso(Spline &delalpha); //!< return likelihood given relative delta alpha(z) / alpha(0)
  double oklo(Cosmos&,Spline &alpha);  //!the zero -result
  double okloLamoreaux(Cosmos&,Spline &delalpha); //! The new non-zero result by Lamoreaux and Torgerson
  double etaChi2(double eta, double log10q); //!< equivalence principle constraint

  /*
  ** Lyman Alpha constraints
  **
  */ 
  double lymanAlphaPatMcDonaldChi2(Cosmos& cosmos); //!< Compute Lyman Alpha constraints using pivot scale k=0.009 s/km. Routines from Pat McDonald

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

};

#endif
