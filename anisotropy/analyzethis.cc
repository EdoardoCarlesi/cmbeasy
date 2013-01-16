#include "analyzethis.h"

#include "newdatreader.h"
#include "spline.h"
// #include "model.h"
#include "miscmath.h"
#include "data.h"
#include "snedata.h"
#include "controlpanel.h"
#include "cl.h"
#include "cosmos.h"
#include "minmax.h"

#include  "LyaFchi2Interpolator.h"


#ifndef NOWMAP7
#include "WMAP_7yr_likelihood.h"
#include "WMAP_7yr_options.h"
#endif

//#include "gsl/gsl_matrix.h"
//#include "gsl/gsl_linalg.h"

#include <numeric>
#include <limits>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm>

#ifdef MACOSX_PANTHER
extern "C" int isnan(double); 
extern "C" int isinf(double); 
#endif 


AnalyzeThis::AnalyzeThis() :  BinnedWMapTT(0), WMAPNotYetInitialized(true) {
  //  double DivideOutMuK = pow(2.725e6,-2.0);

  //Don't read in supernova data yet
  knop03read = riess04read = tonry03read = barris04read = astier05read = false;

  BinnedWMapTT = new Data(ControlPanel::cmbeasyDir("/resources/wmap/wmap_binned_tt_spectrum_7yr_v4p1.txt").c_str(),"binned map tt",Data::Yerror);
  BinnedWMapTE = new Data(ControlPanel::cmbeasyDir("/resources/wmap/wmap_binned_te_spectrum_7yr_v4p1.txt").c_str(),"binned map te",Data::Yerror);
  QSO =  new Data(ControlPanel::cmbeasyDir("/resources/qso/cps_bin_z.dat").c_str(),"qso",Data::Yerror);
  printit=false;
  // post process binned wmap te data (only used for the quick tests and consitency checks)
  for (list<DataEntry>::iterator i = BinnedWMapTE->points.begin(); i != BinnedWMapTE->points.end();i++) {
    double l = i->x();
    (*i) *= l*(l+1)/(2*3.141);
  }

  // Old Supernovae routines
  Riess04_gold = new Data(ControlPanel::cmbeasyDir("/resources/sne1a/riess04_gold.dat").c_str(),"riess04_gold");
  Riess06_gold = new Data(ControlPanel::cmbeasyDir("/resources/sne1a/riess06_gold.dat").c_str(),"riess06_gold");
  Riess04_all = new Data(ControlPanel::cmbeasyDir("/resources/sne1a/riess04_all.dat").c_str(),"riess04_all");
  Riess = new Data(ControlPanel::cmbeasyDir("/resources/sne1a/riess.dat").c_str(),"riess",Data::Yerror);
  Tonry = new Data(ControlPanel::cmbeasyDir("/resources/sne1a/tonry.dat").c_str(),"tonry");
  Tonry_subsample172 = new 
    Data(ControlPanel::cmbeasyDir("/resources/sne1a/tonry_subsample172.dat").c_str(),"tonry_subsample172");
  Tonry_subsample130 = new 
    Data(ControlPanel::cmbeasyDir("/resources/sne1a/tonry_subsample130.dat").c_str(),"tonry_subsample130");
  IfADS = new Data(ControlPanel::cmbeasyDir("/resources/sne1a/IfADS.dat").c_str(),"IfADS");
  HSTSn1a = new Data(ControlPanel::cmbeasyDir("/resources/sne1a/hstsn1a.dat").c_str(),"hstsn1a");
  HSTLowZ = new Data(ControlPanel::cmbeasyDir("/resources/sne1a/hst_low_z.dat").c_str(),"hst_low_z");
  // sort Ries04 (for plotting reasons only, it's nicer to have z re-arranged internally
  Riess04_gold->sort();
  Riess06_gold->sort();
}

AnalyzeThis::~AnalyzeThis() {
  if (BinnedWMapTT) delete BinnedWMapTT;
  if (BinnedWMapTE) delete BinnedWMapTE;
  if (QSO) delete QSO;
}

/* ************************************************************************************************************************************
**
**
** WMAP 5-year stuff
**
**
**
************************************************************************************************************************************** */

#ifdef WMAPTHREEFORTRAN
extern "C"
{
  void __wmap_pass2__pass2_compute_likelihood(  double*, double*, double*, double*, double* );
  //void wmap_options_mp_set_WMAP_data_dir_( char* ); //hehe, passing char* around is a tiny bit more tricky :)
  void __wmap_options__set_ttmax(  int* );
  void __wmap_options__set_ttmin(  int* );
  void __wmap_options__set_temax(  int* );
  void __wmap_options__set_temin( int* );
  void __wmap_options__set_lowl_max( int* );
  void __wmap_options__set_use_lowl_tt( int* );
  void __wmap_options__set_use_lowl_pol( int* );
  void __wmap_options__set_use_tt( int* );
  void __wmap_options__set_use_te( int* );
  void __wmap_options__set_sz_amp( double* );
  void __wmap_util__wmap_likelihood_error_occured( int* );
  void __wmap_util__wmap_likelihood_error_report();
}

bool AnalyzeThis::WMAP3F90Error()
{
  int tmp;
  __wmap_util__wmap_likelihood_error_occured( &tmp );
  return !tmp; // tmp is wmap_likelihood_ok, so return !tmp
}

void AnalyzeThis::WMAP3F90ErrorReport()
{
  __wmap_util__wmap_likelihood_error_report();
}

void AnalyzeThis::WMAP3F90SetMaxTT( int l )
{
  int tmp;
  tmp = l;
  __wmap_options__set_ttmax( &tmp );
}

void AnalyzeThis::WMAP3F90SetMinTT( int l )
{
  int tmp;
  tmp = l;
  __wmap_options__set_ttmin( &tmp );
}

void AnalyzeThis::WMAP3F90SetMaxTE( int l )
{
  int tmp;
  tmp = l;
  __wmap_options__set_temax( &tmp );
}

void AnalyzeThis::WMAP3F90SetMinTE( int l )
{
  int tmp;
  tmp = l;
  __wmap_options__set_temin( &tmp );
}

void AnalyzeThis::WMAP3F90SetLowlMax( int l )
{
  int tmp;
  tmp = l;
  __wmap_options__set_lowl_max( &tmp );
}

void AnalyzeThis::WMAP3F90SetUseLowlTT( bool b )
{
  int tmp; tmp = b;
  __wmap_options__set_use_lowl_tt( &tmp );
}

void AnalyzeThis::WMAP3F90SetUseLowlPol( bool b )
{
  int tmp; tmp = b;
  __wmap_options__set_use_lowl_pol( &tmp );
}

void AnalyzeThis::WMAP3F90SetUseTT( bool b )
{
  int tmp; tmp = b;
  __wmap_options__set_use_tt( &tmp );
}

void AnalyzeThis::WMAP3F90SetUseTE( bool b )
{
  int tmp; tmp = b;
  __wmap_options__set_use_te( &tmp );
}

void AnalyzeThis::WMAP3F90SetSZAmplitude( double amp )
{
  double a; a = amp;
  __wmap_options__set_sz_amp( &a );
}

AnalyzeThis::WMAP7Likelihood AnalyzeThis::WMAP3F90computeLikelihood( CL& cl )
{
  double tt[ 2001 ], te[ 2001 ], ee[ 20001 ],bb[ 20001 ];
  double like[ 9 ];

  Spline TotalT( cl.ts[ 0 ],"total cl" );
  for ( int i = 0; i < cl.ts[ 0 ]->size(); i++ )
    TotalT.set( cl.ts[ 0 ]->y( i ) + cl.tt[ 0 ]->y( i ) ); // scalar + tensor TT C_l's

  TotalT.arm();
  cl.cs[0]->arm();
  cl.es[0]->arm();
  cl.bs[0]->arm();

  //indexing for the fortran arrays starts at 2
  for ( int i =0; i<= 2000-2; ++i )
  {
    tt[ i ] = TotalT.fastY(i+2.);
    te[ i ] = cl.cs[ 0 ]->fastY(i+2.);
    ee[ i ] = cl.es[ 0 ]->fastY(i+2.);
    bb[ i ] = cl.bs[ 0 ]->fastY(i+2.);
  }

  TotalT.disarm();
  cl.cs[0]->disarm();
  cl.es[0]->disarm();
  cl.bs[0]->disarm();

/* Debug/Test:
  ifstream fiducialCls( "/home/tep1/robbers/wmap/data/wmap_likelihood/data/fiducial_test_cls.dat" );
  for ( int i =0; i<= 2000-2; ++i )
  {
    char coma; int num;
    fiducialCls >> num >> coma >> tt[ i ] >> coma >> te[ i ] >> coma >> ee[ i ] >> coma >> bb[ i ];
  }
*/

  __wmap_pass2__pass2_compute_likelihood(&tt[ 0 ], &te[ 0 ], &ee[ 0 ], &bb[ 0 ], &like[ 0 ] );

  WMAP3Likelihood lh;
  lh.TTlike     =  2.*like[-1+1];   //    ! master tttt chisq
  lh.TTlowllike =  2.*like[-1+2];   //    ! low tttt chisq
  lh.TTlowldet  =  2.*like[-1+3];   //    ! low tttt determinant
  lh.Beamlike   =  2.*like[-1+4];   //    ! beam/pt source correction to tttt chisq
  lh.TElike     =  2.*like[-1+5];   //    ! master tete chisq
  lh.TEdet      =  2.*like[-1+6];   //    ! master tete determinant
  lh.Lowllike   =  2.*like[-1+7];   //    ! TE/EE/BB lowl chisq
  lh.Lowldet    =  2.*like[-1+8];   //    ! TE/EE/BB lowl determinant
  lh.total = 2.*std::accumulate( like, like + 8, 0. );

//X   double total = 0;
//X   for ( int i = 0; i < 9; ++i )
//X   {
//X     cout << i << " is " << ( 2.*like[ i ] ) <<  " original: " << like[ i ] << endl;
//X     total += like[ i ];
//X   }
//X   cout << "and total: " << ( 2*total ) << endl;
  cout << "ai.total: " << lh.total << endl;

  return lh;
}
#endif


#ifndef NOWMAP7
AnalyzeThis::WMAP7Likelihood AnalyzeThis::WMAP7computeLikelihood( CL& cl )
{
  unsigned int len = (unsigned int)cl.ts[ 0 ]->stop();
  real8_1d tt(len), te(len), ee(len), bb(len), like(num_WMAP+1);

  Anchor a;
  Spline *ttTot, *teTot, *eeTot, *bbTot;

  ttTot = cl.createTotalTTSpline(0, "totalTT", &a);
  teTot = cl.createTotalTESpline(0, "totalTE", &a);
  eeTot = cl.createTotalEESpline(0, "totalEE", &a);
  bbTot = cl.createTotalBBSpline(0, "totalBB", &a);

  ttTot->arm();
  teTot->arm();
  eeTot->arm();
  bbTot->arm();

  for ( unsigned int i =2; i<len; ++i )
  {
    tt(i)= ttTot->fastY(i);
    te(i)= teTot->fastY(i);
    ee(i)= eeTot->fastY(i);
    bb(i)= bbTot->fastY(i);
  }

//X /* Debug/Test: */
//X   ifstream fiducialCls( "data/fiducial_test_cls.dat" );
//X   do
//X   {
//X     char coma; int i;
//X     fiducialCls >> i; fiducialCls >> coma >> tt( i ) >> coma >> te( i ) >> coma >> ee( i ) >> coma >> bb( i );
//X   } while ( fiducialCls );

  wmap_likelihood_7yr::wmap_likelihood_compute(tt, te, ee, bb, like );

#ifdef USE_HIGHELL_TB
#error AnalyzeThis::WMAP7computeLikelihood: not adapted to be used with TB...
#endif
#ifdef USE_LOWELL_TBEB
#error AnalyzeThis::WMAP7computeLikelihood: not adapted to be used with TBEB...
#endif

  like(0) = 0.;
  WMAP7Likelihood lh;
  lh.TTlike     = 2.*like(1);   //    ! master tttt chisq
  lh.TTlowllike =  2.*like(2);   //    ! low tttt chisq
  lh.TTlowldet  =  2.*like(3);   //    ! low tttt determinant
  lh.Beamlike   =  2.*like(4);   //    ! beam/pt source correction to tttt chisq
  lh.TElike     =  2.*like(5);   //    ! master tete chisq
  lh.TEdet      =  2.*like(6);   //    ! master tete determinant
  lh.Lowllike   =  2.*like(7);   //    ! TE/EE/BB lowl chisq
  lh.Lowldet    =  2.*like(8);   //    ! TE/EE/BB lowl determinant
#if defined USE_HIGHELL_TB
  lh.TBlike     =  2.*like(9);
  lh.TBdet      =  2.*like(10);
#elif defined LOWELL_TBEB
  lh.TBlike     =  2.*like(9);
  lh.TBdet      =  2.*like(10);
#else
  lh.TBlike     =  0.;           // master tbtb chisq flag
  lh.TBdet      =  0.;           // master tbtb determinant flag
#endif
  lh.total = 2*like.sum();

//X   cout << lh.TTlike   << endl; 
//X   cout << lh.TTlowllike<< endl;
//X   cout << lh.TTlowldet << endl;
//X   cout << lh.Beamlike  << endl;
//X   cout << lh.TElike    << endl;
//X   cout << lh.TEdet     << endl;
//X   cout << lh.Lowllike  << endl;
//X   cout << lh.Lowldet   << endl;
//X   cout << "wmp3 year (c++) total: " << lh.total << endl;

  return lh;
}
#endif


/*!
  Scale scalar spectra by a factor 
  (with spectral index position n). 
  \param cl: The Cl's spectra in a struct
  \param n: Number of spectral index
  \param normalize: Factor multiplying the spectra
*/

void AnalyzeThis::scaleScalarCls(CL &cl, int n, double normalize) {
  (*cl.ts[n]) *= normalize;
  (*cl.es[n]) *= normalize;
  (*cl.cs[n]) *= normalize;
}

/*!
  Scale Tensor spectra by a factor 
  (with spectral index position n). 
  \param cl: The Cl's spectra in a struct
  \param n: Number of spectral index
  \param normalize: Factor multiplying the spectra
*/

void AnalyzeThis::scaleTensorCls(CL &cl, int n, double normalize) {
  (*cl.tt[n]) *= normalize;
  (*cl.et[n]) *= normalize;
  (*cl.bt[n]) *= normalize;
  (*cl.ct[n]) *= normalize;
}


/*!
  Scale all spectra by a factor (with scalar index position n). 
  \param cl: The Cl's spectra in a struct
  \param n: Number of spectral index
  \param normalize: Factor multiplying the spectra
*/
void AnalyzeThis::scaleCls(CL &cl, int n, double normalize) {
  (*cl.ts[n]) *= normalize;
  (*cl.es[n]) *= normalize;
  (*cl.cs[n]) *= normalize;
  (*cl.tt[n]) *= normalize;
  (*cl.et[n]) *= normalize;
  (*cl.bt[n]) *= normalize;
  (*cl.ct[n]) *= normalize;
  (*cl.kk[n]) *= normalize;
  (*cl.tk[n]) *= normalize;
}

double AnalyzeThis::qso(Spline &delalpha) {
  //  ofstream monitor("alpha.dat");
  //  double today = alpha(0); // alpha today
  double chi2 = 0;
  for (list<DataEntry>::iterator i=QSO->points.begin(); i!= QSO->points.end(); i++) {
    double z= i->x();
    double meassured = i->y();
    double Delta = i->dy();
    double theory = 1e5*delalpha(z) ;
    chi2 += (meassured-theory)*(meassured-theory)/(Delta*Delta);
    //    monitor << z << " " << theory*1e-5 << " " << meassured*1e-5 << " " << Delta*1e-5 << endl;
  }
  return chi2;
}


/*!
  Return chi2 of new Oklo results from Lamoreaux. 
  
  \param Cosmos a cosmos
  \param DelApha Spline with relative change of alpha
*/
double AnalyzeThis::okloLamoreaux(Cosmos& cosmos, Spline &DelAlpha) {
  double InMpc = cosmos.t_0() -  cosmos.s2mpc(1.8e9*365*24*3600);
  double z = cosmos.t2z(InMpc); // oklo redshift (something like 0.12-0.14 typically)
  
  cout << "OKLO redshift: " << z << endl;
  double theory = DelAlpha(z);  
  double sigma = 7e-9;

  if (theory > 45e-9) sigma = 15e-9;
  cout << " del at this redshift: " << theory << " and sigma: " << sigma << endl;
  
  return (theory - 45e-9)*(theory - 45e-9)/(sigma*sigma);
  
}

double AnalyzeThis::oklo(Cosmos& cosmos, Spline &DelAlpha) {
  double InMpc = cosmos.t_0() -  cosmos.s2mpc(1.8e9*365*24*3600);
  double z = cosmos.t2z(InMpc); // oklo redshift (something like 0.12-0.14 typically)
  
  // cout << "OKLO redshift: " << z << endl;
  double Delta = 1e-7;
  double theory = DelAlpha(z);
  return theory*theory/(Delta*Delta);
}

/*! Return bound from equivalence principle. Eta is the theoretical value obtained
  without "theoretical uncertanties". This eta is allowed to vary by mulitplying it 
  with  10^(-log10q) to 10^(log10q). The best value is then compared to
  the experimental result (which currently is 0).

  In other words: If your theoretical result is eta but there is an uncertainty such
  that eta could be 

  eta_really = eta_theory*(1 + Q) 

  then the eta_really which lies between

  eta_theory*10^(-log10q) .. eta_theory*10^(log10q)

  and fits best the data will be taken for the chi^2.
  
  Currently, of course, the best value is always the one with the lowest
  absolut value...
*/
double AnalyzeThis::etaChi2(double eta, double log10q) {
  double mineta = eta*pow(10.0,-log10q);
  return mineta*mineta * pow(3e-13,-2);
}

/*!
  Replacement for the old CobeNormalize of CMFAST. Uses the binned 
  TT-C_l's from WMAP to normalize the C_l spectrum to this experiment.
  Data is quoted in muK^2

  It is mandatory that rescaleSpectra has been called before this routine,
  because rescaleSpectra() plugs in the necessary factors of pi, h, etc. 
  and should be used to ensure the correct tensor to scalar initial Amplitudes,
  i.e. A_s and A_t
  
*/
void AnalyzeThis::quickWMAPNormalize(Cosmos& cosmos,const ControlPanel& control,CL& cl) {
  if (cosmos.PowerNormalization[0] == 1.0)
    throw Bad_Error("AnalyzeThis::quickWMAPNormalize() Call rescaleSpectra() first.\nSee driver.cc or xdriver.cc for an example of correct setup");

  cosmos.sigma8.resize(cosmos.InitialPower.size());
  for (unsigned int n = 0; n < cosmos.InitialPower.size(); ++n) {
    Spline TotalT(cl.ts[n],"total cl");
    for (int i = 0; i < cl.ts[n]->size(); i++) 
      TotalT.set(cl.ts[n]->y(i) + cl.tt[n]->y(i)); // scaler + tensor TT C_l's 

    TotalT.arm(); // use this 

    // pair<double,double> fit = bestFitting(TotalT,&Binned,899); 
    pair<double,double> fit = bestFitting(TotalT,BinnedWMapTT,899); 
    //    cout << "THIS IS THE TT-FIDUCIAL FIT: " << fit.second << endl;


    double normalize = fit.first;
    scaleCls(cl,n,normalize);
    cosmos.PowerNormalization[n] *= normalize;

    Anchor anchor;  // convenience
    if (control.power_cdm && control.cmb) {
       // Spline* s = cosmos.createPower(n, "quickMapNormalize",cosmos.power_cdm(), &anchor);
      //cosmos.sigma8[n] = sigma8(s,cosmos.h()); // get sigma8 
      cosmos.sigma8[n] *= sqrt(normalize);
    }
  }
}


/*!
  Return factor and a fiducial chi^2 (don't use!) which gives the best normalization
  of the Cl-TT spectrum w.r.t. the binned WMAP data
*/
pair<double,double> AnalyzeThis::bestFittingBinnedWMAP_TT(Spline &ts) {
  return bestFitting(ts,BinnedWMapTT,899); 
}
/*!
  Find factor that, when multiplying Spline ts, leads to the best fit
  to the Data data. Returns this factor and the chi^2 times the
  factor ScaleDof

  If FindFactor is false, then we won't shift the curve but calculate the
  chi^2 without multiplying.
 */
pair<double,double> AnalyzeThis::bestFitting(Spline& ts,Data* data, int ScaleDof,bool FindFactor) {
  double f = 1.0; // just in case that we are not supposed to find the best factor
  if (FindFactor) {
    double Sum1 = 0, Sum2=0;
    for (list<DataEntry>::iterator i=data->points.begin(); i!= data->points.end(); i++) {
      double l = i->x();
      double meassured = i->y();
      double Delta = i->dy();
      Sum1 += ts(l)*meassured / (Delta*Delta);
      Sum2 += ts(l)*ts(l) / (Delta*Delta);
    }
    f = Sum1 / Sum2;
  }

  double chi2 = 0;
  int dof=0;
  for (list<DataEntry>::iterator i=data->points.begin(); i!= data->points.end(); i++) {
    double l = i->x();
    double meassured = i->y();
    double Delta = i->dy();
    //bf << l << " " << meassured << " " << f*ts(l) << " " << f << endl;
    //    cout << l << "  contributing: " <<  pow(fabs(meassured-f*ts(l)), 2)/(Delta*Delta) << " m: " << meassured << " the: " << f*ts(l) << " Del: " << Delta << endl;
    chi2 += pow(fabs(meassured-f*ts(l)), 2)/(Delta*Delta);
    dof++;
  }
  return pair<double,double>(f, chi2*ScaleDof/dof);
}

/*!
  Make a quick estimate based on binned data of the chi^2 (this estimate is more than crude 
  and has not too much to do with later chi^2 values). It will, however filter out models that
  are incredibly bad fitting. For these, even the diagonal likelihood from WMAP may go 
  beserk. 

  Assumes quickNormalized() splines in CL.
*/
bool AnalyzeThis::isSane(CL &cl, int n, double *tt_chi2, double *te_chi2) {
  Spline tt(*cl.ts[n],"isSane_TT");
  Spline te(*cl.cs[n],"isSAne_TE");
  tt.arm();
  te.arm();

  *tt_chi2 = bestFitting(tt,BinnedWMapTT,899,false).second;
  *te_chi2 = bestFitting(te,BinnedWMapTE,449,false).second;

  cout << "isSane checking: " <<  *tt_chi2  << "  :: " << *te_chi2  << endl;
   
  if (isnan(*tt_chi2) || isnan(*te_chi2)) {
    *tt_chi2 = 1e100;
    *te_chi2 = 1e100;
    return false;
  }

  if (*te_chi2 > 2e3) return false;
  if (*tt_chi2 > 4e3) return false;
  return true;
}

/*!
  Normalizes the spectrum with respect to WMAP TT and TE data. First it calls 
  quickWMAPNormalize(). From this first (and very good guess), it checks if 
  the full likelihood will return a sensible value (bad fitting models can 
  have (wrong) and very good  likelihoods in the full covariance routines). 
  If the full routine cannot be used, it returns the diagonal covariance 
  version (which is quite a good estimate).
  If, however, the full routine can be used, it will find the best fit to 
  the full likelihood and return the chi^2. The factor needed to scale the 
  spectrum *after* quickWMAPNormalize() to the best full likelihood fit will 
  be returned. This factor is naturally very close to 1.0 (as the 
  quickWMAPNormalize() is very good already).

  Call this function with un-armed Cl's. It will multiply the normalization factor into   the Cl's rightaway. 

  \param fast If true, then only quickWMAPNormalize() will be used to normalize and the (diagonal only) chi^2 for this normalization will be returned
  
*/
vector<WMAPNorm> AnalyzeThis::WMAPNormalize(Cosmos& cosmos,ControlPanel& control,CL &cl,bool fast ) {
  vector<WMAPNorm> result(cosmos.InitialPower.size());  // our return vector
  // if the spectrum has not yet been quickNormalized(), call quickNormalize
  quickWMAPNormalize(cosmos,control,cl);
  for (unsigned int n = 0; n < cosmos.InitialPower.size(); n++) {
    result[n].norm  = 1.0;  // initialize
    if (isSane(cl,n,&result[n].chi2_tt,&result[n].chi2_te)) {  // in case that it is insane, we already have the values due to this call ...
      Spline like(10, "like");  // so the model is sane...
      Spline cs(*cl.cs[n],"mapn_cs");
      Spline es(*cl.es[n] ,"mapn_es");
      Spline TotalT(cl.ts[n],"total cl");
      for (int i = 0; i < cl.ts[n]->size(); i++) 
	TotalT.set(cl.ts[n]->y(i) + cl.tt[n]->y(i)); // scaler + tensor TT C_l's 

      TotalT.arm(); // use this 
      cs.arm(); es.arm();  

      double diag_tt = WMAPDiag_TT(TotalT) ;
      double diag_te = WMAPDiag_TE(TotalT,cs,es);
      double diag_chi2 = diag_tt + diag_te;
 
      result[n].chi2_tt = diag_tt;
      result[n].chi2_te = diag_te;
      // if chi2 > 1500, i.e. not so well  fitting model, 100% exact likelihood not sooo important
      // and, of course full likelihood is unstable for bad fitting models,hence stick with diagonal only !
      // therefore, the block below is only called for nicely fitting models
      if (diag_chi2 < 1500 && (!fast)) {  // nicely fitting and not fast requested
	// So Diagonal Likelihood is not too bad, let's get the true one around the 
	// factor that did fit best in the diagonal version
	TotalT.disarm(); cs.disarm(); es.disarm();
	double fine_step = 0.05; 
	bool no_nan=true;
	for (double f = 0.95; f <= 1.0501;f += fine_step) {
	  cout << "Full WMAP Likeli at f: " << f << endl;
	  TotalT *= f; cs*=f; es*=f;
	  TotalT.arm(); cs.arm(); es.arm();
	  double tt = WMAP_TT(TotalT);
	  double te = WMAP_TE(TotalT,cs,es);
	  if (isnan(tt) || isnan(te)) {
	    no_nan=false;
	    break;
	  }
	  like.set(f,tt+te);
	  TotalT.disarm(); cs.disarm(); es.disarm();
	  TotalT *= 1.0/f;  es *= 1.0/f; cs *= 1.0/f;
	}
	if (no_nan) {
	  like.arm(); 
	  double x = like.minimum();  // this is the normalization that fits best for temperature
	  
	  TotalT *= x;  cs *=x; es*=x; // rescale ts
	  TotalT.arm(); cs.arm(); es.arm();    
	  result[n].chi2_tt =  WMAP_TT(TotalT);
	  result[n].chi2_te = WMAP_TE(TotalT,cs,es);	
	  result[n].norm = x;
	  scaleCls(cl,n,x);
	  cosmos.PowerNormalization[n] *= x;
	  cosmos.sigma8[n] *= sqrt(x);
	}
      }
    } 
  }
  return result;
}



vector<WMAPNorm> AnalyzeThis::WMAPLikelihood(Cosmos& cosmos,ControlPanel& control,CL &cl) {
  vector<WMAPNorm> result(cosmos.InitialPower.size());  // our return vector
  for (unsigned int n = 0; n < cosmos.InitialPower.size(); n++) {
    bool isok=true;
    result[n].norm  = 1.0;  // initialize

    Spline cs(*cl.cs[n],"mapn_cs");
    Spline es(*cl.es[n] ,"mapn_es");
    Spline TotalT(cl.ts[n],"total cl");
    for (int i = 0; i < cl.ts[n]->size(); i++) 
      TotalT.set(cl.ts[n]->y(i) + cl.tt[n]->y(i)); // scalar + tensor TT C_l's 

    cs.arm(); es.arm();   TotalT.arm();
    double diag_tt = 1.5*bestFitting(TotalT,BinnedWMapTT,899,false).second;  
    double diag_te = 1.5*bestFitting(cs,BinnedWMapTE,449,false).second;  //50 % penalty
    
    //    cout << "bestfit: " << diag_tt/1.5 << " diag_te: " << diag_te/1.5 << endl;
    result[n].chi2_tt = 1e60;
    result[n].chi2_te = 1e60; // initialize

    double diag_chi2 = diag_tt + diag_te;
    if (isnan(diag_chi2)) isok = false;
    
    // next step, if naive counting (plus penalty) is indicating a 
    // model for which diag likelihood will return sensible values,
    // we do this
    if (diag_tt < 6000 && diag_te < 2000 && isok) {
      diag_tt = 1.1*WMAPDiag_TT(TotalT) ;  // 10 % penalty
      diag_te = 1.1*WMAPDiag_TE(TotalT,cs,es);   // 10% penalty
      diag_chi2 = diag_tt + diag_te;
      // cout << "diag: " << diag_tt/1.1 << " diag_te: " << diag_te/1.1 << endl;
    }
    if (isnan(diag_chi2)) isok = false;
   
    if (isok) {
      result[n].chi2_tt = diag_tt;
      result[n].chi2_te = diag_te;
    }
    // if chi2 (including 10% penalty for diagonal estimate)  is good enough,
    // we go on to calculate the true wmap full likelihood 
    // in fact, for a montecarlo chain, this will be the practically always used
    // for models in the chain after initial descend to good parameter region
    if (diag_tt < 1300 && diag_te < 600 && isok) {  // nicely fitting and not fast requested
      // So Diagonal Likelihood is not too bad, let's get the true one 
      double tt = WMAP_TT(TotalT);
      double te = WMAP_TE(TotalT,cs,es);
      //      cout << "full: " << tt << " diag_te: " << te << endl;
      if (! isnan(tt) ) result[n].chi2_tt = tt;
      if (! isnan(te) ) result[n].chi2_te = te;      
    }
  }
  return result;
}

void AnalyzeThis::fiducialAmplitudes(Cosmos& cosmos,vector<double> &A_s,vector<double> &A_t) {
  int N =  cosmos.InitialPower.size();
  A_s.resize(N);
  A_t.resize(N);
  for (int n = 0; n < N; n++) { 
    A_s[n] = 1.0;
    A_t[n]= 1.0;
  }
}

/*!
  Apply Inflationary relation for A_t in terms of A_s and n_T. 
  Make sure that the A_s and A_t vectors are of the right size (i.e. # of spectral indices to compute) 
  and A_s contains the Amplitude you want. If you want to WMAP normalize later, A_s = 1 is
  efectly fine. This you can get conveniently from fiducialAmplitdues
*/

void AnalyzeThis::applyInflationaryTensorRatio(Cosmos &cosmos,vector<double> &A_s, vector<double> &A_t) {
  if (A_s.size() != cosmos.InitialPower.size()) throw Bad_Error("AnalyzeThis::TensorRatio A_s vector size not sufficient for all requested spectral indices");
  if (A_t.size() != cosmos.InitialPower.size()) throw Bad_Error("AnalyzeThis::TensorRatio A_t vector size not sufficient for all requested spectral indices");
  for (unsigned int n = 0; n < cosmos.InitialPower.size(); ++n) {
      A_t[n] = -8*A_s[n]*cosmos.InitialTensorPower[n];
      //      cout << "n_T : " << cosmos.InitialTensorPower[n] << endl;
  }
}

/*!
  Rescale the Cl's and Powerspectra by a constant factor A_s (i.e. the initial 
  scalar power amplitude). Due to conventions of cmbfast, cmbeasy and camb,
  there are some factors of pi etc needed. 
  
  NOTE: If you simply want to rescale the C'ls, use scaleCls() 

  In our (and camb's) convention, 

  sigma_8 = 1/(2pi^2) \int dlnk W(x)^2 k^3 P_cdm(k) 

  P_cdm is related to P_chi (our initial spectrum) with chi=-1 at initial times for each mode

  P_cdm(k) = 2pi^2 P_chi(k)/k^3 delta_cdm^2 

  where it is understood that delta_cdm_intial is scaled such that chi=-1 for each mode
  The Cl's follow from 

  C_l  = 2/pi \int dk k^2 P_cmbeasy * Delta_l^2 

  But P_cmbeasy in this formula is P_cmbeasy defined as P_cdm(k), i.e. an additional factor of 2pi^2
  comes in for the Cl's also

  Finally, the output is in terms of l*(l+1)C_l/(2*pi). So yet another factor...

  Please note that sigma_8 is also computed for convenience. 
*/

void AnalyzeThis::rescaleSpectra(Cosmos& cosmos, const ControlPanel& control, CL& cl,
                                 vector<double>& A_s, vector<double>& A_t, bool force)
{
  if (A_s.size() != cosmos.InitialPower.size()) {
    throw Bad_Error("AnalyzeThis::RescaleSpectrum() A_s vector size not sufficient for all requested spectral indices");
  }
  if (A_t.size() != cosmos.InitialPower.size()) {
    throw Bad_Error("AnalyzeThis::rescaleSpectrum()  A_t vector size not sufficient for all requested spectral indices");
  }
  cosmos.sigma8.resize(cosmos.InitialPower.size());
  for (unsigned int n = 0; n < cosmos.InitialPower.size(); ++n) {
    // we factor in some more usumal factors, such that
    // you may just call Cosmos::dumpPower() and get a publishable
    // spectrum :-)

    // last minute sanity check: has this normalization been carried out before ?
    // if so, we would introduce wrong factors, so throw error
    if (!force && cosmos.PowerNormalization[n] != 1.0)
      throw Bad_Error("AnalyzeThis::rescaleSpectra() Power spectrum has been normalized already.");

    double pi=M_PI;
    //  A factor of h^3 is factored in, in order to plot P(k) in units of Mpc h^-3
    // In addition, to make agree with the definition of P(k) in terms of the intial spectrum P_chi(k),
    // we have to multiply by a factor of 2*pi^2
    cosmos.PowerNormalization[n] =  2*pi*pi* A_s[n] * pow(cosmos.h(),3);

    // For the Cl's,  we have some other factors

    if (control.cmb) {
      double normalizeCl = 1/(4*pi*4*pi);  // undo cmbfast's 4*pi, i.e. do only the integral
      normalizeCl *= 2/pi; // apply the Cl measure of cmbeasy (following my lecture)
      normalizeCl *= 2*pi*pi; // conversion to camb's P_chi - power spectrum definition
      normalizeCl /= 2*pi; // plot is in terms of  l*(l+1)/2*pi

      if (control.scalar && control.cmb) scaleScalarCls(cl,n,normalizeCl*A_s[n]);   // scale scalars propto A_s

      normalizeCl /= 16; // CMBFASTs P_h is 1/16th of CAMB'S & CMBEASY's P_h
      if (control.tensor && control.cmb) scaleTensorCls(cl,n,normalizeCl*A_t[n]);  // tensors propto A_t
    }

    if (control.power_cdm) {
      Anchor anchor;  // convenience
      Spline* s = cosmos.createPower(n, "quickMapNormalize",cosmos.power_cdm(), &anchor);
      cosmos.sigma8[n] = sigma8(s,cosmos.h()); // get sigma8 
    }
  }
}

//! same as above, applies the same A_s and A_t to all spectrac cl::xx[i]
void AnalyzeThis::rescaleSpectra(Cosmos& cosmos, const ControlPanel& control, CL& cl,
                                 double A_s, double A_t, bool force)
{
  vector<double> A_sVec(cosmos.InitialPower.size(), A_s);
  vector<double> A_tVec(cosmos.InitialPower.size(), A_t);
  rescaleSpectra(cosmos, control, cl, A_sVec, A_tVec, force);
}

/*!
  chi2 of cluster abundance constraint, see Wang & Steinhardt

double AnalyzeThis::clusterAbundance(Model &m, const double weff) {
  
  double theta = m.n -1.0 + m.h -0.65;
  double gamma_ = 0.21 - 0.22*weff + 0.33*(m.o_cdm + m.o_b) + 0.25*theta;
  
  double s8o =  m.sigma8 * pow((double)(m.o_cdm + m.o_b), gamma_);

  double delta = 0.1;

  return (s8o - 0.5)*(s8o - 0.5) / (delta*delta);
}


double AnalyzeThis::sigma8(Model &m) {
  Anchor anchor;
  return sigma8(m.power_cdm(&anchor), m.h);
}
*/

/*!
  Given power spectrum spline s as a function of  k / h, 
  and in units of Mpc h^-3, (i.e. the one from createPower() in cosmos),
  return sigma8, i.e. the density contrast on 8 Mpc/h scales
*/
double AnalyzeThis::sigma8(Spline* s,double h) {
  tmpSig8 = s;
  if (!s->isArmed()) s->arm();
  //  s->dump("sigma8.dat");

  return sqrt(Miscmath::rombint((moSingle)&AnalyzeThis::sig8Integrator, this, s->start(), s->stop(),1e-6));
}


double AnalyzeThis::sig8Integrator(const double koverh) const  {  
  static double tpi =1.0/(2*M_PI*M_PI); // 1/(2*pi^2)    
  double r=8;
  double x = koverh*r;
  double win = 3*( sin(x) - x*cos(x) ) / (x*x*x);

  // PowerNormalization  is propto h^3, hence we need to correct for this here
  // in addition, our integral was over k/h and not k and we pick up an additional
  // factor of h^3. The two cancel, so we don't have to do anything.
  double power = tmpSig8->fastY(koverh); 
  return tpi*power*pow(koverh*win,2);
}



//
//
//  Sophisticated Supernovae Routines 
//  Thanks to Alex Conley
//

/*!
  For Knop03 there are three subsets
  - Primary
  - Low-extinction primary (E(B-V) > 0.1 at 2\sigma removed)
  - Low-extinction strict: Low-extinction + most stringent spectral type cut
  - In addition, we consider the Primary with only the HST high-z SNe
    and the Low-extinction with only the HST high-z SNe
*/
int AnalyzeThis::readKnop03Data() {

  //Get basic sample
  Knop03Extended.readData(ControlPanel::cmbeasyDir("/resources/sne1a/conley/knop03_extended.dat"));
  Knop03Extended.setName("Knop03Extended");

  //Sort because it's nicer to look at them that way
  Knop03Extended.zcmbsort();

  //Now form low-extinction sample
  const int n_lowe_remove = 4;
  string lowe_remove[n_lowe_remove] = { "1992ag","1993ag","1998as","1998ax" };
  vector<string> lowe_remove_v(lowe_remove,lowe_remove+n_lowe_remove);
  Knop03ExtendedLowe = Knop03Extended.copy_remove(lowe_remove_v,true);
  Knop03ExtendedLowe.setName("Knop03ExtendedLowExtinction");

  //From that, get strict-Ia sample
  const int n_lowest_remove = 6;
  string lowest_remove[n_lowest_remove] = 
    { "1995as", "1996cf", "1996cg", "1996cm", "1998ay", "1998be" };
  vector<string> lowest_remove_v(lowest_remove,lowest_remove+n_lowest_remove);
  Knop03ExtendedLoweStrict = 
    Knop03ExtendedLowe.copy_remove(lowest_remove_v, true);
  Knop03ExtendedLoweStrict.setName("Knop03ExtnededLowExtinctionStrictIa");

  //Then make only HST sample
  const int n_hst_remove = 27;
  string hst_remove[n_hst_remove] =
    { "1995ar", "1995as", "1995aw", "1995ax", "1995ay", "1995ax", "1995ay",
      "1995az", "1995ba", "1996cf", "1996cg", "1996ci", "1996cl", "1996cm",
      "1997f", "1997h", "1997i", "1997n", "1997p", "1997q", "1997r",
      "1997ac", "1997af", "1997ai", "1997aj", "1997am", "1997ap" };
  vector<string> hst_remove_v(hst_remove,hst_remove+n_hst_remove);
  Knop03HST = Knop03Extended.copy_remove(hst_remove_v,true);
  Knop03HST.setName("Knop03HST");

  //And HST Lowe sample
  Knop03HSTLowe = Knop03HST.copy_remove(lowe_remove_v, false);
  Knop03HSTLowe.setName("Knop03HSTLowExtinction");

  knop03read = true;

  return 0;
}

/*!
  For Riess04 there are two subsets -- gold, and all (gold+silver)
*/
int AnalyzeThis::readRiess04Data() {

  //Read in full set
  Riess04all.readData(ControlPanel::cmbeasyDir("/resources/sne1a/conley/riess04.dat"));
  Riess04all.setName("Riess04GoldAndSilver");

  //Sort because it's nicer to look at them that way
  Riess04all.zcmbsort();

  //Then remove silver ones
  const int n_silver = 29;
  string silver[n_silver]={"1994B","1994C","1995E","1995M","1995ap","1995ao",
			   "1995ae","1995as","1995ar","1996R","1996T","1996V",
			   "1996cg","1996cm","1996cf","1997ck","1998ay",
			   "1998be","1999da","1999fh", "2000ce","2000ea",
			   "2001jb","2001kd","2002P","2002ab", "2002ad",
			   "2002fx","2002kc" };
  vector<string> silver_v(silver,silver+n_silver);
  Riess04gold = Riess04all.copy_remove(silver_v,true);
  Riess04gold.setName("Riess04Gold");

  riess04read = true;
  return 0;
}

/*!
  For Tonry03 there are three subsets -- everything, 172 (with z < 0.01
   and A_V <= 0.5), and 130 (like 172 but with Perlmutter99 SNe removed)
 */
int AnalyzeThis::readTonry03Data() {

  //Read in base sample
  Tonry03full.readData(ControlPanel::cmbeasyDir("/resources/sne1a/conley/tonry03.dat"));
  Tonry03full.setName("Tonry03AllSNe");

  //Sort because it's nicer to look at them that way
  Tonry03full.zcmbsort();

  //Then remove the lowz and heavily extinguished ones
  const int n172 = 58;
  string tonry172[n172]={"sn72E","sn80N","sn81B","sn81D","sn86G","sn89B",
			 "sn90N","sn90Y","sn91M","sn91T","sn91bg","sn92A",
			 "sn92G","sn92K","sn94D","sn94U","sn94ae","sn95D",
			 "sn95E","sn95ak","sn95al","sn95bd","sn96X","sn96Z",
			 "sn96ai","sn96ao","sn96bk","sn96bo","sn96bv","sn97bd",
			 "sn97bp","sn97bq","sn97br","sn97bz","sn97cw","sn97dt",
			 "sn97fb","sn98ab","sn98aq","sn98bu","sn98de","sn98dh",
			 "sn98dk","sn98dm","sn98ea","sn98ec","sn98es","sn99ac",
			 "sn99by","sn99cl","sn99da","sn99ek","sn99fh","sn99gd",
			 "sn99gh","sn00ce","sn00cx","sn00ea"};
  vector<string> tonry172_v(tonry172,tonry172+n172);
  Tonry03_172 = Tonry03full.copy_remove( tonry172_v, true );
  Tonry03_172.setName("Tonry03:172SNe");

  //And remove Perlmutter 99 SNe
  const int np99 = 42;
  string p99[np99]={"sn92bi","sn94F","sn94G","sn94H","sn94al","sn94am",
		    "sn94an","sn95aq","sn95ar","sn95as","sn95at","sn95aw",
		    "sn95ax","sn95ay","sn95az","sn95ba","sn96cf","sn96cg",
		    "sn96ci","sn96ck","sn96cl","sn96cm","sn96cn","sn97F",
		    "sn97G","sn97H","sn97I","sn97J","sn97K","sn97L","sn97N",
		    "sn97O","sn97P","sn97Q","sn97R","sn97S","sn97ac","sn97af",
		    "sn97ai","sn97aj","sn97am","sn97ap"};
  vector<string> p99_v(p99,p99+np99);
  Tonry03_130 = Tonry03_172.copy_remove( p99_v, true );
  Tonry03_130.setName("Tonry03:130SNe");

  tonry03read = true;

  return 0;
}

/*!
  For Barris04 there are no subsamples.  A 'strict-Ia' subsample is
   discussed in the paper, but the data is then combined with the Tonry03_172
   sample, for which no such cut has been performed.  It would also be
   possible to limit the Barris sample to only those from the Barris paper
   in addition to some of the low-z SNe from Tonry03, but this has not
   been done here.

   This function will also read in the Tonry03 SNe if they haven't already
   been read in.
*/
int AnalyzeThis::readBarris04Data() {

  Barris04.readData(ControlPanel::cmbeasyDir("/resources/sne1a/conley/barris04.dat"));
  if (tonry03read == false) readTonry03Data(); //We need the Tonry SNe
  Barris04 += Tonry03_172;

  //Sort because it's nicer to look at them that way
  Barris04.zcmbsort();

  barris04read = true;
  return 0;
}

/*!
  For Astier 05 there are no subsamples.
*/
int AnalyzeThis::readAstier05Data() {
  Astier05.readData(ControlPanel::cmbeasyDir("/resources/sne1a/conley/astier05.dat"));

  //Sort because it's nicer to look at them that way
  Astier05.zcmbsort();

  astier05read = true;
  return 0;
}


/*!
  Extended core routine for computing SNe Ia likelihoods.
  This version supports full fitting of the stretch-luminosity
  and colour-luminosity relations.
  While, for many purposes, SNe Ia don't care about the Hubble
  constant, for convenience this uses the value set in cosmos to work with
  the absolute magnitude of SN Ia.  Returns the chisquare of
  the fit.

  Uses the technique of Goliath et al., 2001 A&A 380, p. 6-18 to 
  analytically marginalize over the magnitude offset (sometimes called
  script-M in supernova circles).  Note that this is only
  possible because the errors do not depend on this quantity.

  \param cosmos A basecosmos derived class that will provide distances 
         and derivates
  \param data The Supernovae redshift - magnitude data
  \param alpha The slope of the stretch-luminosity relation
  \param beta The slope of the colour-luminosity relation
  \param intrinsicdisp The Assumed intrinsic 'dispersion' of SN1a in magnitudes
  \param pecz Peculiar velocity in redshift units -- for example, 
    if the peculiar velocity is 300 km s^{-1}, this should be about 0.001
  \param widthmean Assumed mean of width parameter 
       ( stretch: 1 delta m15: 1.1 ). Not really necessary, but convenient
  \param fixerr If this is set to true, then fixerr_alpha and fixerr_beta
          are used for error propogation purposes, while alpha and beta
          are used to correct the magnitudes.  This approach is sometimes
          known as a 'stepwise chisquare' fit.  This can be useful if the
          data does not strongly constrain these parameters, because they
          can have a tendency to run away to large values since raising
          their values causes the errors on the individual data points to
          become larger.  The idea is then to do multiple runs of the 
          cosmological fitter until alpha and beta converge.
          Note that there is currently no way to fix
          only one of these parameters.
  \param fixerr_alpha Fixed value of alpha used for error propagation.
          Has no effect unless fixerr is set.
  \param fixerr_beta  Fixed value of beta used for error propagation.
          Has no effect unless fixerr is set.
  \returns The chisquare of the fit
*/
double AnalyzeThis::SNIaCore(const baseCosmos& cosmos,
			     const SNeData& data, double alpha,
			     double beta, double intrinsicdisp,
			     double pecz, double widthmean,
			     bool fixerr=false, double fixerr_alpha=0.0,
			     double fixerr_beta=0.0) const {
  const double prefac = 5 / log(10.0);
  const double twopi = 2*3.14159265358979323846264338327950288419716939937510582;

  //Avoid the complications of rebounding universes
  if ( cosmos.isUniverseRebounding() )
    return numeric_limits<double>::infinity();

  if ( ! cosmos.validHistory() ) 
    throw Bad_Error("AnalyzeThis::SNIaCore() History is not valid");

  int nsn = data.size();
  double *dl = new double[ nsn ];
  for (int i = 0; i < nsn; ++i)
    dl[i] = cosmos.luminosityDistance( data[i].zhel, data[i].zcmb );

  double isq = intrinsicdisp*intrinsicdisp;

  //One complication of the Goliath marginalization equations is that,
  // as presented in the paper, they rely on subtracting large
  // numbers from each other, which is generally something worth
  // avoiding if possible.
  //In order to avoid this issue, we form an approximate guess for
  // scriptm and evaluate around that instead of zero, which
  // only causes slight modifications in the formulae.
  //The rough estimate is formed from the weighted mean without
  // properly calculated errors
  double scriptm0, wtval, roughinvvar, deltaM;
  wtval = scriptm0 = 0.0;
  for (int i = 0; i < nsn; ++i) {
    deltaM = data[i].dmag;
    roughinvvar = 1.0 / (isq + deltaM*deltaM);
    wtval += roughinvvar;
    scriptm0 += (data[i].mag - 5*log10(dl[i])) * roughinvvar;
  }
  scriptm0 /= wtval;

  double peczsq = pecz*pecz;

  double GoliathA, GoliathB, GoliathC, invvar;
  double zcmb, diffmag, deltaZ, deltaZmag, mTheory;
  double corrMag, corrVar;
  GoliathA = GoliathB = GoliathC = 0.0;
  for (int i = 0; i < nsn; ++i) {

    zcmb = data[i].zcmb;
    deltaZ = data[i].dz;
    deltaZ = sqrt( deltaZ*deltaZ + peczsq ); //Add in peculiar velocity

    deltaM = data[i].dmag;

    //Check to make sure dl is positive.  How could it be negative?
    // For certain types of closed Universes, there is enough time
    // for light to travel half way across the Universe by the current
    // time.  This leads to infinite magnification at the antipode.
    // Redshifts at larger distances have formally negative luminosity
    // distances.  See Zel'dovich and Novikov, Relativisitic Astrophysics,
    //  volume 2, 1983, University of Chicago Press.
    //Since we are taking a log, this is a problem.
    //In this code, this situation is handled by returning an infinite chisq

    if (dl[i] < 0) return numeric_limits<double>::infinity();

    mTheory = scriptm0 + 5 * log10(dl[i]);

    //In some ways it would be preferrable to correct the theoretical
    // magnitude rather than the observed magnitude.  However, this
    // is a purely pedological distinction in this case, as the
    // mathmatics is identical either way.  And the organization of
    // the code is a little cleaner this way.
    //So corrMag is the magnitude of the SN corrected for stretch and 
    // colour (whether you think colour is due to extinction or not)
    corrMag = data[i].getCorrMag( alpha, beta, widthmean );
    diffmag = corrMag - mTheory;

    // the delta in redshift translates into one on the predicted magnitude...
    deltaZmag = prefac * cosmos.luminosityDistanceDeriv(zcmb) * 
      deltaZ/dl[i];

    if (fixerr) {
      corrVar = data[i].getCorrVar( fixerr_alpha, fixerr_beta );
    } else {
      corrVar = data[i].getCorrVar( alpha, beta );
    }

    //Add in extra pieces for redshift errors and intrinsic 'dispersion'
    double DeltaSquare = corrVar + deltaZmag*deltaZmag + isq;

    invvar = 1.0 / DeltaSquare;

    GoliathA += diffmag * diffmag * invvar;
    GoliathB += diffmag * invvar;
    GoliathC += invvar;

    //Debugging output
    /*
    double corrmag = data[i].mag + alpha*(data[i].widthpar - widthmean) -
      beta * data[i].colourpar;
    printf("SN: %-6s zcmb: %5.3f mag: %6.2f dmag: %4.2f  magcorr: %6.2f err: %4.2f\n",
	   data[i].name.c_str(),zcmb,data[i].mag,deltaM,corrmag,
	   sqrt(DeltaSquare));
    */
  }

  delete[] dl;

  return GoliathA + log(GoliathC/twopi) - GoliathB*GoliathB / GoliathC;

}

/*!
  \param cosmos A basecosmos derived class that will provide distances and 
      derivates
  \param data The Supernovae redshift - magnitude data
  \param alpha The slope of the stretch-luminosity relation
  \param beta The slope of the colour-luminosity relation
  \param intrinsicdisp The Assumed intrinsic dispersion of SN1a in magnitudes
  \param pecz Peculiar velocity in redshift units -- for example, 
    if the peculiar velocity is 300 km s^{-1}, this should be about 0.001
  \param widthmean Assumed mean of width parameter 
       ( stretch: 1 delta m15: 1.1 ). Not really necessary, but convenient
  \returns The estimated value of the scriptM nuisance parameter

  Note that the definition of ScriptM here is slightly different
  from that used in most papers.
*/
double AnalyzeThis::estimateScriptM(const baseCosmos& cosmos,
				    const SNeData& data,
				    double alpha, double beta,
				    double intrinsicdisp, double pecz,
				    double widthmean) const {

  const double prefac = 5 / log(10.0);

  double dl, wt, wtmean, mTheory, varcurr, deltaZ, deltaZmag, corrMag;
  double isq = intrinsicdisp * intrinsicdisp;

  wt = wtmean = 0.0;
  for (SNeData::const_iterator sn = data.begin(); sn != data.end(); ++sn) {
    dl = cosmos.luminosityDistance(sn->zhel, sn->zcmb);
    mTheory = 5*log10( dl );
    corrMag = sn->getCorrMag( alpha, beta );

    deltaZ = sn->dz;
    deltaZ = sqrt( deltaZ*deltaZ + pecz*pecz ); //Add in peculiar velocity
    deltaZmag = prefac * cosmos.luminosityDistanceDeriv(sn->zcmb) * 
      deltaZ/dl;

    varcurr = sn->getCorrVar( alpha, beta ) + isq + deltaZmag*deltaZmag;

    wt += 1.0 / varcurr;
    wtmean += (corrMag - mTheory) / varcurr;

  }

  return wtmean / wt;

}

/*! This function computes the chi^2 statistic of the SN Ia Supernova
  dataset of Knop et al. 2003.  Unlike Sn1aHST, this version allows
  M_B, alpha and R_B to be varied, and allows the user to specify the
  peculiar velocity and intrinsic dispersion used.

  Right now this does not use the covariance matrix between different
  SNe.

  \param cosmos A basecosmos derived class that will provide distances and 
      derivates
  \param data The Supernovae redshift - magnitude data
  \param alpha The slope of the stretch-luminosity relation
  \param Rb The dust R_B. 4.1 corresponds to normal galactic dust
  \param intrinsicdisp The Assumed intrinsic dispersion of SN1a in magnitudes
  \param pecz Peculiar velocity in redshift units -- for example, if 
      the peculiar velocity is 300 km s^{-1}, this should be about 0.001
  \returns The chi^2 of the fit.
*/
double AnalyzeThis::SNIaKnop03(const baseCosmos& cosmos,double alpha, 
			       Knop03Sample sample, double Rb, 
			       double intrinsicdisp, const double pecz) {

  if (knop03read == false) readKnop03Data();

  SNeData const* data = 0;
  switch (sample) {
  case K03Primary :
    data = &Knop03Extended;
    break;
  case K03Lowe :
    data = &Knop03ExtendedLowe;
    break;
  case K03Lowestrict :
    data = &Knop03ExtendedLoweStrict;
    break;
  case K03HST :
    data = &Knop03HST;
    break;
  case K03HSTLowe :
    data = &Knop03HSTLowe;
    break;
  }

  //cout << "Knop03 used " << data->size() << " supernovae" << endl;

  //Knop03 assumes that colour is due to dust.  The mean stretch
  // value is around 1
  return SNIaCore(cosmos,*data,alpha,Rb,intrinsicdisp,pecz,1.0);

}

/*!
  \param cosmos A basecosmos derived class that will provide distances and 
      derivates
  \param data The Supernovae redshift - magnitude data
  \param alpha The slope of the stretch-luminosity relation
  \param Rb The dust R_B. 4.1 corresponds to normal galactic dust
  \param intrinsicdisp The Assumed intrinsic dispersion of SN1a in magnitudes
  \param pecz Peculiar velocity in redshift units -- for example, if the peculiar velocity is 300 km s^{-1}, this should be about 0.001
  \returns The estimate for scriptM
*/
double AnalyzeThis::estimateScriptMKnop03(const baseCosmos& cosmos,
					  double alpha, Knop03Sample sample,
					  double Rb, double intrinsicdisp,
					  const double pecz) {
  if (knop03read == false) readKnop03Data();

  SNeData const* data = 0;
  switch (sample) {
  case K03Primary :
    data = &Knop03Extended;
    break;
  case K03Lowe :
    data = &Knop03ExtendedLowe;
    break;
  case K03Lowestrict :
    data = &Knop03ExtendedLoweStrict;
    break;
  case K03HST :
    data = &Knop03HST;
    break;
  case K03HSTLowe :
    data = &Knop03HSTLowe;
    break;
  }
  return estimateScriptM(cosmos,*data,alpha,Rb,intrinsicdisp,
			 pecz,1.0);
}


/*! This function computes the chi^2 statistic of the SN Ia Supernova
  dataset of Astier et al. 2005.  

  Unlike most of the other fits available, the user must specify
  a value for both the stretch-magnitude slope (alpha) and the
  colour-magnitude slope (beta).  In the Astier05 paper these were fit,
  with best fit values of 1.52466 and 1.56900, respectively.  However,
  to properly replicate the parameters of this paper you must fit
  for these yourself (or your errors will be underestimated), so
  these values are not provided as defaults.

  Right now this does not use the covariance matrix between different
  SNe.
  
  To use this in a monte carlo, leave alpha and beta as monte carlo parameters.
  Say they are on position 8 and 9 in the task-array. Then a call would look like 
  this:
  	double astier = ai.SNIaAstier05(cosmos,task[8],task[9]);
	task[LOGLIKEPOS+8] =  -0.5 *astier;

  \param cosmos A basecosmos derived class that will provide distances and 
      derivates
  \param data The Supernovae redshift - magnitude data
  \param alpha The slope of the stretch-luminosity relation
  \param beta The slope of the colour-magnitude relation.
  \param intrinsicdisp The Assumed intrinsic dispersion of SN1a in magnitudes.
    The default of 0.131, as used in the paper, is provided
  \param pecz Peculiar velocity in redshift units -- for example, if 
      the peculiar velocity is 300 km s^{-1}, this should be about 0.001.
      The default is 0.001
  \returns The chi^2 of the fit.

*/
double AnalyzeThis::SNIaAstier05(const baseCosmos& cosmos,double alpha, 
				 double beta, double intrinsicdisp, 
				 const double pecz) {
  if (astier05read == false) readAstier05Data();

  //cout << "Astier05 used " << Astier05.size() << " supernovae" << endl;

  //The mean stretch value is around 1
  return SNIaCore(cosmos,Astier05,alpha,beta,intrinsicdisp,pecz,1.0);

}

/*!
  \param cosmos A basecosmos derived class that will provide distances and 
      derivates
  \param data The Supernovae redshift - magnitude data
  \param alpha The slope of the stretch-luminosity relation
  \param Rb The dust R_B. 4.1 corresponds to normal galactic dust
  \param intrinsicdisp The Assumed intrinsic dispersion of SN1a in magnitudes
  \param pecz Peculiar velocity in redshift units -- for example, if the peculiar velocity is 300 km s^{-1}, this should be about 0.001
  \returns The estimate for scriptM
*/
double AnalyzeThis::estimateScriptMAstier05(const baseCosmos& cosmos,
					    double alpha, double beta, 
					    double intrinsicdisp,
					    const double pecz) {
  if (astier05read == false) readAstier05Data();

  return estimateScriptM(cosmos,Astier05,alpha,beta,intrinsicdisp,
			 pecz,1.0);
}

/*! This function computes the chi^2-statistic of the SN Ia Supernovae 
 dataset of Riess et al. 2004.  Note that some of the SN are also
 contained in the other routines.  Furthermore, this fit is somewhat
 limited compared to those possible with some of the other data sets because
 relatively little information about the lightcurve fits has been
 published -- i.e., there is no correlation information, and the
 errors on the extinction estimates are not available, making it difficult
 to evaluate any changes in the extinction law.  Many of the SN in this 
 sample do not possess extinction estimates.

\param cosmos A basecosmos derived class that will provided distances and
   derivatives
\param  goldsample If true, then only the data from the gold sample will 
  be used, otherwise silver objects will be inlcluded
\returns The chi^2 of the fit
*/
double AnalyzeThis::SNIaRiess04(const baseCosmos& cosmos, 
				bool goldsample) {

  if (riess04read == false) readRiess04Data();

  //The intrinsic dispersion is already included in the Riess table
  //This data assumes that dust in external galaxies is identical to
  // ours, and that there is no additional colour parameter independent
  // of the lightcurve shape.  Therefore, Riess has converted colour
  // measurements into A_B, so the beta coefficient is defined to be one 
  // for this sample.  The MLCS delta parameter is not available,
  // and has been pre-applied to the data, so alpha is zero

  //Note that the A_B values have been 'unapplied' in the read in
  // data table in order to allow ambitious fitters to try different
  // values of R_B, etc., although the lack of error information makes
  // this a dicey proposition.

  if ( goldsample ) {
    //cout << "Riess 04 used " << Riess04gold.size() << " supernovae" << endl;
    return SNIaCore(cosmos,Riess04gold,0.0,1.0,0.0,0.0,0.0);
  } else {
    //cout << "Riess 04 used " << Riess04all.size() << " supernovae" << endl;
    return SNIaCore(cosmos,Riess04all,0.0,1.0,0.0,0.0,0.0);
  }
}

/*! This function estimates the scriptM parameter for the Riess04
  data set.
  \param cosmos A basecosmos derived class that will provided distances and
   derivatives
  \param  goldsample If true, then only the data from the gold sample will 
    be used, otherwise silver objects will be inlcluded
  \returns The estimated value of scriptM
*/
double AnalyzeThis::estimateScriptMRiess04(const baseCosmos& cosmos,
					   bool goldsample) {

  if (riess04read == false) readRiess04Data();
  if (goldsample) {
    return estimateScriptM(cosmos,Riess04gold,0.0,1.0,0.0,0.0,0.0);
  } else return estimateScriptM(cosmos,Riess04all,0.0,1.0,0.0,0.0,0.0);
}


/*! This function computes the chi^2-statistic of the SN Ia Supernovae 
 dataset of Tonry et al. 2003  Note that some of the SN are also
 contained in the other routines.  Furthermore, this fit is somewhat
 limited compared to those possible with some of the other data sets because
 relatively little information about the lightcurve fits has been
 published (correlation, errors on extinction).  Many of the SN in this 
 sample do not possess extinction estimates.

\param cosmos A basecosmos derived class that will provided distances and
   derivatives
\param sample Determines which subsample to use.  TonryFull uses all SNe from
   Table 8 of Tonry, Tonry172 excludes all SNe with z < 0.01 or A_V > 0.5,
  and Tonry130 also removes all of the SN from Perlmutter et al. 1999.
\returns The chi^2 of the fit
*/
double AnalyzeThis::SNIaTonry03(const baseCosmos& cosmos, 
				Tonry03Sample sample) {
  const double pecz = 500.0 / 299792.458; //500 km/sec / c

  if (tonry03read == false) readTonry03Data();

  //The intrinsic 'dispersion' has already been applied to the
  // magnitudes in this data set.
  //This data assumes that dust in external galaxies is identical to
  // ours, and that there is no additional colour parameter independent
  // of the lightcurve shape.  Therefore, Tonry has converted colour
  // measurements into A_V, so the beta coefficient is defined to be R_B/R_V
  // for this sample.  The MLCS delta and/or dm15 parameters are not available,
  // and have been pre-applied to the data, so alpha is zero
  
  //Note that the A_V values have been 'unapplied' in the read in
  // data table in order to allow ambitious fitters to try different
  // values of R_V, etc., although the lack of error information makes
  // this a dicey proposition.

  switch (sample) {
  case TonryFull :
    //cout << "Tonry 03 used " << Tonry03full.size() << " supernovae" << endl;
    return SNIaCore(cosmos,Tonry03full,0.0,4.1/3.1,0.0,pecz,0.0);
    break;
  case Tonry172 :
    //cout << "Tonry 03 used " << Tonry03_172.size() << " supernovae" << endl;
    return SNIaCore(cosmos,Tonry03_172,0.0,4.1/3.1,0.0,pecz,0.0);
    break;
  case Tonry130 :
    //cout << "Tonry 03 used " << Tonry03_130.size() << " supernovae" << endl;
    return SNIaCore(cosmos,Tonry03_130,0.0,4.1/3.1,0.0,pecz,0.0);
    break;
  default :
    throw Bad_Error("AnalyzeThis::SNIaTonry03() Unknown sample selector");
    break;
  }
}

/*! This function estimates the scriptM parameter for the Tonry 03
  data set.  Again, errors on the individual components that comprise
  the corrected magnitude are not available.
  \param cosmos A basecosmos derived class that will provided distances and
   derivatives
  \param sample Determines which subsample to use.  TonryFull uses all SNe from
   Table 8 of Tonry, Tonry172 excludes all SNe with z < 0.01 or A_V > 0.5,
  and Tonry130 also removes all of the SN from Perlmutter et al. 1999.
  \returns The estimated value of scriptM
*/
double AnalyzeThis::estimateScriptMTonry03(const baseCosmos& cosmos,
					   Tonry03Sample sample) {
  const double pecz = 500.0 / 299792.458; //500 km/sec / c

  if (tonry03read == false) readTonry03Data();

  switch (sample) {
  case TonryFull :
    return estimateScriptM(cosmos,Tonry03full,0.0,4.1/3.1,0.0,pecz,0.0);
    break;
  case Tonry172 :
    return estimateScriptM(cosmos,Tonry03_172,0.0,4.1/3.1,0.0,pecz,0.0);
    break;
  case Tonry130 :
    return estimateScriptM(cosmos,Tonry03_130,0.0,4.1/3.1,0.0,pecz,0.0);
    break;
  default :
    throw Bad_Error("AnalyzeThis::estimateScriptMTonry03() Unknown sample selector");
  }
}

/*! This function computes the chi^2-statistic of the SN Ia Supernovae 
 dataset of Barris et al. 2004  Note that some of the SN are also
 contained in the other routines.  Furthermore, this fit is somewhat
 limited compared to those possible with the some of the other data sets 
 because information about the individual error components of the corrected
 magnitudes is not available.
 Many of the SN in this sample do not possess extinction estimates.
 This data sample does not support subsetting, and is combined with the
 Tonry03_172 sample.  SN 2001jn has been removed from the sample because
 of a large amount of extinction (this SN was also removed from the analysis
 in the paper).

\param cosmos A basecosmos derived class that will provided distances and
   derivatives
\returns The chi^2 of the fit
*/
double AnalyzeThis::SNIaBarris04(const baseCosmos& cosmos) {
  const double pecz = 500.0 / 299792.458; //500 km/sec / c

  //Essentially the same comments apply to this data sample as
  // to Tonry 03

  if (barris04read == false) readBarris04Data();
  //cout << "Barris04 used " << Barris04.size() << " supernovae" << endl;
  return SNIaCore(cosmos,Barris04,0.0,4.1/3.1,0.0,pecz,0.0);
}

/*! This function estimates the scriptM parameter for the Barris04
  data set.
  \param cosmos A basecosmos derived class that will provided distances and
   derivatives
  \returns The estimated value of scriptM
*/
double AnalyzeThis::estimateScriptMBarris04(const baseCosmos& cosmos) {
  const double pecz = 500.0 / 299792.458; //500 km/sec / c
  if (barris04read == false) readBarris04Data();
  return estimateScriptM(cosmos,Barris04,0.0,4.1/3.1,0.0,pecz,0.0);
}


//
//
//   Old Supernovae Routines 
//
//
  
/*!
  Core routine for computing Sne Ia likelihoods. It first finds the best
  fit to the data by shifting the effective magnitude up and down (this is
  the reason why there's no need in providing the Hubble parameter: taking
  the log() would just give an addition constant)
  Then, it calculates the chi^2 of the prediction to the data
  \param luminosityDistance A Spline holding the luminosity distance 
  \param data The Supernovae redshift - magnitude data
  \param f factor multiplying log10(d_l). Usually this is 5, yet Tonry needs 1
*/

double AnalyzeThis::Sn1aCore(Spline& luminosityDistance, Data &data,double f) {
  // derivative will be needed to estimate errors due to dz 
  Spline Ddldz(luminosityDistance,"DLdz"); 
  luminosityDistance.derive(Ddldz);
  Ddldz.arm();

  double S =0, Sy =0,chi2=0,FMSy=0;  // Summation variables for finding the best fit
  // 
  //  first step: find the best fitting M_b. I call it scriptM here
  // 
  for (list<DataEntry>::iterator i=data.points.begin(); i != data.points.end(); i++) {
    double z = i->x();
    double m = i->y();
    double deltaZ = i->dx();
    double deltaM = i->dy();
   
    double dl = luminosityDistance(z);  // as log is acting, no need for pre-factors
    double sm = m - f*log10(dl);

    // the delta in redshift translates into one on the predicted magnitude...
    double deltaTheory = f/(dl*log(10.0)) * Ddldz(z)*deltaZ;  
    // total delta^2
    double DeltaSquare = deltaM*deltaM + deltaTheory*deltaTheory;

    S += 1/DeltaSquare;
    Sy += sm/DeltaSquare;

    // For fiducial model output (just monitoring / plotting stuff)
    // as we also have shift freedom, as said: it's for
    // plotting and not data analyzation, so disregard this if you are not
    // using the monitor.dat outputfile below
    double FMdl = (1+z)*log10(1+z);
    double FMsm = m-f*log10(FMdl);
    FMSy += FMsm / DeltaSquare;
  }
  // scriptM is not the real M, because we left H_0 out of the picture
  // however, this is the constant that yields the best fit to the supernovae
  double scriptM = Sy / S; 

  //
  // Last Step: Now that we have fixed the overall normalization, let's calculate the chi^2
  //
  for (list<DataEntry>::iterator i=data.points.begin(); i!= data.points.end(); i++) {
    double z = i->x();
    double m = i->y();
    double deltaZ = i->dx();
    double deltaM = i->dy();
    
    double dl = luminosityDistance(z); 
    double  mt=  scriptM + f*log10(dl); // m_theory
    // the delta in redshift translates into one on the predicted magnitude...
    double deltaTheory = f/(dl*log(10.0)) * Ddldz(z)*deltaZ;  
    // total delta^2
    double DeltaSquare = deltaM*deltaM + deltaTheory*deltaTheory;
    // For fiducial model output (just monitoring / plotting stuff)
    // as we still have shift freedom, we shift, but very crudly, as said: it's for
    // plotting and not data analyzation
    //double FM = f*log10((1+z)*log10(1+z)) + FMscriptM;    
    //    monitor << z << " " << mt << " " << m << "  " << mt - FM << "  " << m - FM << " " << deltaM << endl;
    chi2 += (mt-m)*(mt-m) / DeltaSquare;
  }		
  return chi2;
}



/*! This function computes the chi^2-statistic of the Sn1a Supernovae dataset of 
Riess et. al. 2004, astro-ph/0402512.  Note that some of the supernovae are also contained in the
other routines, hence the sets are not independent. 
We use the data given in Table 5 of Riess et. al. (columns 2,3 and 4).
 The algorithm shifts the 
data up and down to find the minimum likelihood. 
The subsample selector gives only the gold set (see Riess et. al. 2004) if set to true, all if set to false.

\param luminosityDistance  a Spline of a computed luminosity distance for a model universe
\param  goldsample if true, then only the data from the gold sample will be used, otherwise all
*/
double AnalyzeThis::Sn1aRiess04(Spline &luminosityDistance, bool goldsample) {
  Data data("sn1ariess04");

  if (goldsample) data = *Riess04_gold;
  else data = *Riess04_all;
  cout << "Riess04 used " << data.points.size() << " supernovae" << endl;
  return Sn1aCore(luminosityDistance, data);
}

double AnalyzeThis::Sn1aRiess06(Spline &luminosityDistance, bool goldsample) {
  if (!goldsample)
    throw Bad_Error("not implemented");

  Data data("sn1ariess06");

  data = *Riess06_gold;
  //cout << "Riess06 used " << data.points.size() << " supernovae" << endl;
  return Sn1aCore(luminosityDistance, data);
}

/*! This function computes the chi^2-statistic of the Sn1a Supernovae dataset of 
Knop et. al. 2003, astro-ph/0309368. This includes the low-z Supernovae R99 and H96 and the Perlmutter
et.al. Supernovae. Note that some of the supernovaes are also present in the Torny et. al. dataset, 
so the two sets are not independent. We use the data with strech/luminosity slope and host galaxy
extinction correction applied (column 6 in the tables of Knop et. al.). The algorithm shifts the 
data up and down to find the minimum likelihood. The subsample selector can be used to drop a
few supernovae from the set. If subsample=true then this algorithm computes the likelihood of set 2
of Knop et. al. 
\param luminosityDistance  a Spline of a computed luminosity distance for a model universe
\param subsample wheter or not to exclude the supernovae of subsample 2 in Knop et.al. */
double AnalyzeThis::Sn1aHST(Spline& luminosityDistance, bool subsample){
  Data data("sn1acombined");
  int c =0;
  for (list<DataEntry>::iterator i=HSTLowZ->points.begin(); i!= HSTLowZ->points.end(); i++,c++) {
    if(!subsample || (c !=4 && c !=17)) data.points.push_back(DataEntry(i->x(),i->y(),i->dx(),i->dy()));
  }
  c =0;
  for (list<DataEntry>::iterator i=HSTSn1a->points.begin(); i!= HSTSn1a->points.end(); i++,c++) {
    if(!subsample || (c !=3 && c!=4) ) data.points.push_back(DataEntry(i->x(),i->y(),i->dx(),i->dy()));
  }
  cout << "HST USED: " << data.points.size() << " supernovae" << endl;
  return Sn1aCore(luminosityDistance, data);
}
    
/*!
  Likelihood w.r.t to data from Riess et. al. The data is given relative to a 
  certain cosmological model (I momentarily forgot which one, to be honest).
  Nevertheless, this model is correctly subtracted here.
  \param luminosityDistance A Spline containing the luminosity distance
*/
double AnalyzeThis::Sn1aRiess(Spline&  luminosityDistance) {
  double err=0,diff=0;
  for (list<DataEntry>::iterator i=Riess->points.begin(); i!= Riess->points.end(); i++) {
    double z = i->x();
    double m = i->y();
    double deltaM = i->dy();
    err += 1/(deltaM*deltaM);
    double dl = luminosityDistance(z);
    m  += 5*log(z*(z+2))/log(10.0);
    diff += (m -  5*log(dl)/log(10.0)) / (deltaM*deltaM);
  }
  double scriptM = diff / err;  // best constant, this is not scriptM, for I've discarded some other constants
  double chi2 = 0;
  for (list<DataEntry>::iterator i=Riess->points.begin(); i!= Riess->points.end(); i++) {
    double z = i->x();
    double m = i->y();
    double deltaM = i->dy();
    m  += 5*log(z*(z+2))/log(10.0);
    double dl = luminosityDistance(z);
    double theory =  scriptM + 5*log(dl)/log(10.0);
    chi2 += (theory - m)*(theory-m)/(deltaM*deltaM);
  }
  return chi2;
}



void AnalyzeThis::initWMAPCommon(int *Progress) {
  if (WMAPNotYetInitialized) {
    initWMAPCommonTT(Progress);
    initWMAPCommonTE(Progress);
    WMAPNotYetInitialized=false;
  }
}


//   "common_t_likelihood" fills the common block used by
//   "compute_mapt_likelihood"---this is the initialization subroutine.
//
//   Inputs:
//      clFile  - The name of the file containing the cl data.
//      offDiag - The name of the file containing the off-diagonal terms.
//
//   Outputs:
//      stat - A status code:  0=success.
//
//   Common blocks:
//        tt_data - This common block is used to pass/store information
//                  that does not change between calls.
//
//   Written by Licia Verde and Hiranya Peiris, Princeton University,
//      December 2002.

void AnalyzeThis::initWMAPCommonTT(int *Progress) {
  //cout << "COMMONLIKELI" << endl;
  ifstream Diag(ControlPanel::cmbeasyDir("/resources/wmap/tt_diag.dat").c_str());
  ifstream Offdiag(ControlPanel::cmbeasyDir("/resources/wmap/tt_offdiag.dat").c_str());

  if (! Diag || !Offdiag) throw Bad_Error("AnalyzeThis::mapCommonTLikeli(): Data file not found");

  ifstream DiagBin(ControlPanel::cmbeasyDir("/resources/wmap/tt_diag_binary.dat").c_str());
  ifstream OffdiagBin(ControlPanel::cmbeasyDir("/resources/wmap/tt_offdiag_binary.dat").c_str());
 
  ifstream Check(ControlPanel::cmbeasyDir("/resources/wmap/tt_conversion_completed.txt").c_str());

  if (!DiagBin || ! OffdiagBin || !Check) {
    //    cout << "TT: Binary does not exist" << endl;
    ConvertWMAP2BinaryTT(Progress); 
    initWMAPCommonTT(Progress);
    return;
  }    // else  cout << "InitWMAPCOMMONTT::: Binary exists !!!!" << endl;

  //  int dummy;
  //  for (int l = 2; l <= tarrsiz; l++) Diag >> dummy >> cl_data[l] >> neff[l] >> fskyeff[l];
  for (int l = 2; l <= tarrsiz; l++) {
    cl_data[l] = read<double>(DiagBin);
    neff[l] = read<double>(DiagBin);
    fskyeff[l] = read<double>(DiagBin);
  }  
  
  //
  //		Read in off diag terms:		      
  //		In the *covariance* matrix there are 2 type of terms 
  //		one that does scale with the clth (due to the mask) and 
  //		one that does not (due to beam and point sources marginalization).
  // 		In the curvature matrix (what we use here) all the off diagonal 
  //		terms end up scaling with the clth but the 2 contributions scale
  //		in different ways see paper for details
  // 		thus here we read them in separately
  
  // Dynamically resize the arrays.
  off_diag.resize(tarrsiz+1);
  r_off_diag.resize(tarrsiz+1);
  for (int l = 2; l <= tarrsiz; l++) {
    off_diag[l].resize(tarrsiz+1);
    r_off_diag[l].resize(tarrsiz+1);
  }
  int zero=0,nzero=0;
  for (int l = 2; l <= tarrsiz; l++) {
    if (Progress) *Progress = (int) ( (100* l) / (double) tarrsiz ); 
    for (int ll=l+1; ll <= tarrsiz; ll++) {
      //      Offdiag >> dummy >> dummy >> off_diag[l][ll] >> r_off_diag[l][ll];
      off_diag[l][ll] = read<double>(OffdiagBin);
      r_off_diag[l][ll] =   read<double>(OffdiagBin);
      if (off_diag[l][ll] == 0.0 && r_off_diag[l][ll] == 0.0) zero++; else nzero++;
      r_off_diag[ll][l] = r_off_diag[l][ll];
      off_diag[ll][l] = off_diag[l][ll];
    }
  }
}

void AnalyzeThis::ConvertWMAP2BinaryTT(int *Progress) {
  ifstream Diag(ControlPanel::cmbeasyDir("/resources/wmap/tt_diag.dat").c_str());
  ifstream Offdiag(ControlPanel::cmbeasyDir("/resources/wmap/tt_offdiag.dat").c_str());

  if (! Diag || !Offdiag) throw Bad_Error("AnalyzeThis::ConvertWMAP2BinaryTT(): Data file not found");

  ofstream DiagBin(ControlPanel::cmbeasyDir("/resources/wmap/tt_diag_binary.dat").c_str());
  ofstream OffdiagBin(ControlPanel::cmbeasyDir("/resources/wmap/tt_offdiag_binary.dat").c_str());

  int dummy;
  for (int l = 2; l <= tarrsiz; l++) { 
    double cl_data,neff,fskyeff;
    Diag >> dummy >> cl_data >> neff >> fskyeff; 
    write<double>(DiagBin, cl_data);
    write<double>(DiagBin, neff);
    write<double>(DiagBin, fskyeff);
  }
  
  for (int l = 2; l <= tarrsiz; l++) {
    double off_diag,r_off_diag;
        if (Progress) *Progress = (int) ( (100* l) / (double) tarrsiz ); 
    for (int ll=l+1; ll <= tarrsiz; ll++) {
      Offdiag >> dummy >> dummy >> off_diag >> r_off_diag;
      write<double>(OffdiagBin, off_diag);
      write<double>(OffdiagBin, r_off_diag);
    }
  }
  ofstream Check(ControlPanel::cmbeasyDir("/resources/wmap/tt_conversion_completed.txt").c_str());
  Check << "conversion completed" << endl;
}


#define MAP_TE_LTOT 450
void AnalyzeThis::ConvertWMAP2BinaryTE(int *Progress) {
  ifstream Diag(ControlPanel::cmbeasyDir("/resources/wmap/te_diag.dat").c_str());  
  ifstream Offdiag(ControlPanel::cmbeasyDir("/resources/wmap/te_offdiag.dat").c_str());

  if (! Diag || !Offdiag) throw Bad_Error("AnalyzeThis::ConvertWMAP2BinaryTE(): Data file not found");

  
  ofstream DiagBin(ControlPanel::cmbeasyDir("/resources/wmap/te_diag_binary.dat").c_str());
  ofstream OffdiagBin(ControlPanel::cmbeasyDir("/resources/wmap/te_offdiag_binary.dat").c_str());
  

  int dummy;
  for (int l = 2; l <= TE_ARRSIZE; l++) {
    double te_data,cltt_data,ntt,nee;
    Diag >> dummy >> te_data >> cltt_data >> ntt >> nee;
    write<double>(DiagBin,te_data);
    write<double>(DiagBin,cltt_data);
    write<double>(DiagBin,ntt);
    write<double>(DiagBin,nee);
  }
  
  for (int l = 2; l <= TE_ARRSIZE; l++) {
    double te_off_diag;
    if (Progress) *Progress = (int) ( (100* l) / TE_ARRSIZE ); 
    for (int ll=l+1; ll <=  TE_ARRSIZE; ll++) {
      Offdiag >> dummy >> dummy >> te_off_diag;
      write<double>(OffdiagBin,te_off_diag);
    }
  }
  ofstream Check(ControlPanel::cmbeasyDir("/resources/wmap/te_conversion_completed.txt").c_str());
  Check << "conversion completed" << endl;
}



void AnalyzeThis::initWMAPCommonTE(int *Progress) {
  //  cout << "COMMONLIKELI TE" << endl;
  ifstream Diag(ControlPanel::cmbeasyDir("/resources/wmap/te_diag.dat").c_str());  
  ifstream Offdiag(ControlPanel::cmbeasyDir("/resources/wmap/te_offdiag.dat").c_str());

  if (! Diag || !Offdiag) throw Bad_Error("AnalyzeThis::initMap(): Data file not found");
  
  ifstream DiagBin(ControlPanel::cmbeasyDir("/resources/wmap/te_diag_binary.dat").c_str());
  ifstream OffdiagBin(ControlPanel::cmbeasyDir("/resources/wmap/te_offdiag_binary.dat").c_str());
  ifstream Check(ControlPanel::cmbeasyDir("/resources/wmap/te_conversion_completed.txt").c_str());

  if (!DiagBin || ! OffdiagBin || !Check) {
    //    cout << "TE: Binary does not exist" << endl;
    ConvertWMAP2BinaryTE(Progress);
    initWMAPCommonTE(Progress);
    return;
  } //else cout << "InitWMAPCOMMONTE::: Binary exists !!!!" << endl;

  //  int dummy;
  // for (int l = 2; l <= TE_ARRSIZE; l++) Diag >> dummy >> te_data[l] >> cltt_data[l] >> ntt[l] >> nee[l];
  for (int l = 2; l <= TE_ARRSIZE; l++) {
    te_data[l] = read<double>(DiagBin);
    cltt_data[l] = read<double>(DiagBin);
    ntt[l] = read<double>(DiagBin);
    nee[l] = read<double>(DiagBin);
    //    cout << "READIN IN : " << l << " ntt: " << ntt[l] << "  nee: " << nee[l] << endl;
  }
  // Dynamically resize the arrays.
  te_off_diag.resize(TE_ARRSIZE+1);
  for (int l = 2; l <= TE_ARRSIZE; l++) {
    te_off_diag[l].resize(TE_ARRSIZE+1);
  }
 
  //for (int l = 2; l <= MAP_TE_LTOT; l++) {
  //  for (int ll=l+1; ll <= MAP_TE_LTOT; ll++) {
 for (int l = 2; l <= TE_ARRSIZE; l++) {
   if (Progress) *Progress = (int) ( (100* l) / TE_ARRSIZE ); 
    for (int ll=l+1; ll <=  TE_ARRSIZE; ll++) {
      //      Offdiag >> dummy >> dummy >> te_off_diag[l][ll];
      te_off_diag[l][ll] = read<double>(OffdiagBin);
      te_off_diag[ll][l] = te_off_diag[l][ll];
    }
 }
}

/*!      diagon TT map - likelihood: This is to prevent 
 wrong answers from the real likelihood: If a model is really way off, then the
approximations break down and you will get a wrong (possibily positive) log
likelihood back. So here is only the diagonal part. This is fast, albeit not
really correct
*/
double AnalyzeThis::WMAPDiag_TT(Spline &ts) {
 
  int lmax = tarrsiz, lmin = 2;
  
// prepare to compute the offset longormal likelihood as in Bond Jaffe Knox.
// with the difference that the transformation we do on the curvature
// matrix is using clth not cltdata
// this is closer to the equal variance approx (see Bond Jaffe Knox again)  
// form more details see Verde et al .2003 
// here Fdiag is the diagonal term of the covariance matrix 
// Fisher denotes the curvature matrix.

  double Fdiag[tarrsiz+1];

  double z[tarrsiz+1], zbar[tarrsiz+1];

  vector< vector<double> > Fisher(lmax+1), off_log_curv(lmax +1);
  for (int l = 2; l <= lmax; l++) {
    off_log_curv[l].resize(lmax +1);
    Fisher[l].resize(lmax+1);
  }
  //  int larr[tarrsiz+1];
  for (int l = 2; l <= lmax; l++) {
    //larr[l] = l;  // What ? 
    // careful: l in the denominator is in the original larr(l)
    // however, as larr[l] = l, that doesn't matter
    Fdiag[l] = 2*(ts(l) + neff[l])*(ts(l) + neff[l]) / ( (2.0 * l + 1.0) * fskyeff[l]*fskyeff[l]);
    Fisher[l][l] = 1.0/Fdiag[l];
    z[l] = log(cl_data[l] + neff[l]);
    zbar[l] = log(ts(l) + neff[l]);
  }
  for (int l = 2; l <= lmax; l++)   off_log_curv[l][l]=(ts(l)+neff[l]) *Fisher[l][l]*(ts(l)+neff[l]);
    
  double loglike = 0.0,dloglike;
  for (int l = lmin; l <= lmax; l++) {
      // 
      //  to correct for residual 0.5% bias around the peak
      //  see Verde et.al. 2003 for more details
      // 
      //  this is an interpolation between Bond Knox Jaffe and Gaussian Likelihood.
      //  works extremely well on sims (again see paper for details)
      // 
      dloglike = 2./3.*(z[l]-zbar[l])*off_log_curv[l][l] *  (z[l]-zbar[l]);
      dloglike += 1./3.*(ts(l)-cl_data[l]) *Fisher[l][l]*(ts(l)-cl_data[l]);
      loglike +=dloglike;
  }
  // return -loglike/2.0;
  return loglike; // this is chi2
}

double AnalyzeThis::WMAP_TT(Spline &ts) {
 
  int lmax = tarrsiz, lmin = 2;
  
// prepare to compute the offset longormal likelihood as in Bond Jaffe Knox.
// with the difference that the transformation we do on the curvature
// matrix is using clth not cltdata
// this is closer to the equal variance approx (see Bond Jaffe Knox again)  
// form more details see Verde et al .2003 
// here Fdiag is the diagonal term of the covariance matrix 
// Fisher denotes the curvature matrix.

  double Fdiag[tarrsiz+1];

  double z[tarrsiz+1], zbar[tarrsiz+1];

  vector< vector<double> > Fisher(lmax+1), off_log_curv(lmax +1);
  for (int l = 2; l <= lmax; l++) {
    off_log_curv[l].resize(lmax +1);
    Fisher[l].resize(lmax+1);
  }
  //  int larr[tarrsiz+1];
  for (int l = 2; l <= lmax; l++) {
    //larr[l] = l;  // What ? 
    // careful: l in the denominator is in the original larr(l)
    // however, as larr[l] = l, that doesn't matter
    Fdiag[l] = 2*(ts(l) + neff[l])*(ts(l) + neff[l]) / ( (2.0 * l + 1.0) * fskyeff[l]*fskyeff[l]);
    Fisher[l][l] = 1.0/Fdiag[l];
    z[l] = log(cl_data[l] + neff[l]);
    zbar[l] = log(ts(l) + neff[l]);
  }
  for (int l = 2; l <= lmax; l++) {
    for (int ll = 2; ll <= lmax; ll ++) {
      if ( l != ll) {	
	//
	//here the two contributions to the off diagonal terms to the curvature matrix 
	//are treated separately
	//
	//see Verde et.al. 2003 for details.
	//
	Fisher[l][ll] =r_off_diag[l][ll]/sqrt(Fdiag[l]*Fdiag[ll]) + off_diag[l][ll]/(Fdiag[l]*Fdiag[ll]);

      } 
      //cout << "accessing: " << l << "  " << ll << endl;
      off_log_curv[l][ll]=(ts(l)+neff[l]) *Fisher[l][ll]*(ts(ll)+neff[ll]);
    }
  }
  

  double loglike = 0.0,dloglike;
  for (int l = lmin; l <= lmax; l++) {
    for (int ll = lmin; ll <= lmax; ll ++) { 
      // 
      //  to correct for residual 0.5% bias around the peak
      //  see Verde et.al. 2003 for more details
      // 
      //  this is an interpolation between Bond Knox Jaffe and Gaussian Likelihood.
      //  works extremely well on sims (again see paper for details)
      // 
      dloglike = 2./3.*(z[l]-zbar[l])*off_log_curv[l][ll] *  (z[ll]-zbar[ll]);
      dloglike += 1./3.*(ts(l)-cl_data[l]) *Fisher[l][ll]*(ts(ll)-cl_data[ll]);
      loglike +=dloglike;
    }    
  }
  //  return -loglike/2.0;
  return loglike; // this is chi2
}


// tt , te , ee , l0
double AnalyzeThis::WMAP_TE(Spline &ts, Spline &cs, Spline &es) {
  double fsky = 0.85;
  double loglike = 0;
  int lmax = MAP_TE_LTOT;
  
  double Fdiag[MAP_TE_LTOT+1];
  vector< vector<double> > Fisher(lmax+1);

  for (int l = 2; l <= lmax; l++) {
    Fisher[l].resize(lmax+1);
    Fdiag[l] = ((ts(l) + ntt[l])*(es(l) + nee[l]) + cs(l)*cs(l) )  / ( (2.0 * l + 1.0) *fsky*fsky/1.14);
    Fisher[l][l] = 1.0/Fdiag[l];
  }
  for (int l = 2; l<= lmax; l++) {
   for (int ll = 2; ll<= lmax; ll++) {
     if (l != ll) {
       Fisher[l][ll] = te_off_diag[l][ll] / sqrt(Fdiag[l] * Fdiag[ll]);
     }
   }
  }
   for (int l = 2; l<= lmax; l++) {
     for (int ll = 2; ll<= lmax; ll++)  loglike += (cs(l) - te_data[l]) * Fisher[l][ll] * (cs(ll)-te_data[ll]);
   }
   // return -loglike/2.0;
   return loglike; // this is chi2
}

double AnalyzeThis::WMAPDiag_TE(Spline &ts, Spline &cs, Spline &es) {
  double fsky = 0.85;
  double loglike = 0;
  int lmax = MAP_TE_LTOT;
  
  double Fdiag[MAP_TE_LTOT+1], Fisher[MAP_TE_LTOT+1];

  for (int l = 2; l <= lmax; l++) {
    Fdiag[l] = ((ts(l) + ntt[l])*(es(l) + nee[l]) + cs(l)*cs(l) )  / ( (2.0 * l + 1.0) *fsky*fsky/1.14);
    Fisher[l] = 1.0/Fdiag[l];
  }
   for (int l = 2; l<= lmax; l++) {
     loglike += (cs(l) - te_data[l]) * Fisher[l] * (cs(l)-te_data[l]);
   }
 
   //   return -loglike/2.0;
   return loglike; // this is chi2
}


/* ************************************************************************************************************************************
**
**
** SDSS
**
**
**
**
**
************************************************************************************************************************************** */
double AnalyzeThis::SDSS_chiSquared(Spline *powerSpline, double bias, double Kbreak) {
  float window[22][97], Pk_conv[22], k_orig[97];
  
  //initialize
  for (int i=0; i<22; i++) Pk_conv[i]=0;
  
  // read in the window matrix, a 22 x 97 matrix
  // and the k-Values corresponting to the columns and lines 
  ifstream winFile(ControlPanel::cmbeasyDir("/resources/sdss/sdss_windows.dat").c_str());
  ifstream kfile(ControlPanel::cmbeasyDir("/resources/sdss/sdss_kbands.dat").c_str());
  
  if (! winFile) throw Bad_Error("AnalyzeThis::SDSS_chiSquared() cannot find window matrix sdss_windows.dat" );
  if (!kfile) throw Bad_Error("AnalyzeThis::SDSS_chiSquared() cannot find kbands sdss_kbands.dat" );

  for (int j=0; j<97; j++){
    kfile >> k_orig[j];
  }  

  for(int i=0; i<22; i++) {
    for (int j=0;j<97 ; j++){
      winFile >> window[i][j];
    }
  }

  //Multiply the data with the window matrix
  for (int i=0;i<22; i++){
    for (int j=0; j<97; j++){
      Pk_conv[i]+=window[i][j]*SDSS_Pk(k_orig[j],powerSpline);
    } 
  }

  // so we have in Pk_conv[] now the convoluted powerspectrum, which we can compare with the data.
  // read the data
  ifstream dataFile(ControlPanel::cmbeasyDir("/resources/sdss/sdss_measurements.dat").c_str());
  if (!dataFile) throw Bad_Error("AnalyzeThis::SDSS_chiSquared() cannot find data sdss_measurements.dat" );

  float PkExp[22], deltaPkExp[22], k[22],dummy;
  // please note that in the file, as well as in the theoretical spectrum, we store k/h
  // so the k[i] read below is in fact k/h
  for (int i=0; i<22; i++){
    dataFile >> k[i] >> dummy >> dummy >> PkExp[i] >> deltaPkExp[i] >> dummy;
  }


  // compute chi^2 statistic
  double chi2=0;
  int counter=0;
  for(int i=0; i < 22 ; i++ ){
    if(Kbreak < k[i]) continue;
    chi2+=(PkExp[i]-bias*Pk_conv[i])*(PkExp[i]-bias*Pk_conv[i])/(deltaPkExp[i]*deltaPkExp[i]);
    counter++;
  }
  
  cout << "SDSS: " << counter << " data points: " << chi2 << " chi2" << endl;
  return chi2;
}

double AnalyzeThis::SDSS_bestBias(Spline *powerSpline, double Kbreak){
  float window[22][97], Pk_conv[22], k_orig[97];
  
  //initialize
  for (int i=0; i<22; i++) Pk_conv[i]=0;
  
  // read in the window matrix, a 22 x 97 matrix
  // and the k-Values corresponting to the columns and lines 
  ifstream winFile(ControlPanel::cmbeasyDir("/resources/sdss/sdss_windows.dat").c_str());
  ifstream kfile(ControlPanel::cmbeasyDir("/resources/sdss/sdss_kbands.dat").c_str());
  
  if (! winFile) throw Bad_Error("AnalyzeThis::SDSS_chiSquared() cannot find window matrix sdss_window.dat" );
  if (!kfile) throw Bad_Error("AnalyzeThis::SDSS_chiSquared() cannot find kbands sdss_kbands.dat" );

  for (int j=0; j<97; j++){
    kfile >> k_orig[j];
  }  

  for(int i=0; i<22; i++) {
    for (int j=0;j<97 ; j++){
      winFile >> window[i][j];
    }
  }

  //Multiply the data with the window matrix
  for (int i=0;i<22; i++){
    for (int j=0; j<97; j++){
      Pk_conv[i]+=window[i][j]*SDSS_Pk(k_orig[j],powerSpline);
    } 
  }

  // so we have in Pk_conv[] now the convoluted powerspectrum, which we can compare with the data.
  // read the data
  ifstream dataFile(ControlPanel::cmbeasyDir("/resources/sdss/sdss_measurements.dat").c_str());
  if (!dataFile) throw Bad_Error("AnalyzeThis::SDSS_chiSquared() cannot find data sdss_measurements.dat" );

  float PkExp[22], deltaPkExp[22], k[22],dummy;
  for (int i=0; i<22; i++){
    dataFile >> k[i] >> dummy >> dummy >> PkExp[i] >> deltaPkExp[i] >> dummy;
  }

  // find best bias
  double sum1=0;
  double sum2=0;
  double bestBias;
  for(int i=0;i<22; i++){
    if(Kbreak < k[i]) continue;
    sum1+=PkExp[i]*Pk_conv[i]/(deltaPkExp[i]*deltaPkExp[i]);
    sum2+=Pk_conv[i]*Pk_conv[i]/(deltaPkExp[i]*deltaPkExp[i]) ;
  }
  
  bestBias=sum1/sum2;
  //  cout << "SDSS BestBias " << bestBias << endl;
  return bestBias;
}


/*! Return the theoretical power spectrum or 0, if the value is
  not within the boundaries of the theoretical spline.
  \paramters k: the wavenumber (actually, this is k/h)
*/
double AnalyzeThis::SDSS_Pk(float k, Spline *powerSpline)  
{ 
  // excluded as the large k values contribute almost nothing in the window function
  /*
  if (k > powerSpline->stop()){
    double PatMaxK= powerSpline->fastY(powerSpline->stop());
    return PatMaxK*pow(powerSpline->stop(),3)/pow(k,3);
  }
  */
  
  if(k< powerSpline->start()) return 0.0;
  if (k > powerSpline->stop()) return 0.0;
  return (*powerSpline)(k);
}


/* ************************************************************************************************************************************
**
**
**  2dF stuff
**
**
**
**
**
************************************************************************************************************************************** */


/*!
  Use cdm powerspectrum spline p and return 
  bias parameter that leads to best fit to 2dF-data
*/
double AnalyzeThis::TwoDF_bestBias(Spline *p) {
  Spline pwr(*p,"analyzethis::2dfpower");
  double best=1.0;
  double best_chi2 = 1e5;
  
  double left = 0.6, right = 1.5001;
  double b_step = (right-left)/9;
  for (int i =0; i < 2; i++) {
    Spline like(100,"analyzethis::2dfpowerlike");
    for (double b = left; b <= right; b += b_step) {
      pwr *= b*b;
      pwr.arm();
      double chi2 = TwoDF_convolutePowerSpectrum(&pwr);
      like.set(b, chi2);
      pwr.disarm();
      pwr *= 1/(b*b);
    }
    like.arm();
    best = like.minimum();
    best_chi2 = like(best);
    left = best - b_step;
    right  = best+b_step;
    b_step *= 0.1;
  }
  
  cout << "best bias: " << best << endl;
  cout << "likelihood of best bias: " << best_chi2 << endl;

  return best;
}

double AnalyzeThis::TwoDF_convolutePowerSpectrum( Spline *powerSpline) {  
  float PkmcW, k_conv[33], k_orig[101], win[33][101], pk_orig[101];
  int i,j;
  vector<double> convTheoreticalData(32), k_value(32);
  double convData;

  /* read in window matrix and relevant k-values */
  ifstream winFile(ControlPanel::cmbeasyDir("/resources/2df/2dFGRS_pk_win.txt").c_str());

  if(!winFile)  throw Bad_Error("AnalyzeThis::convolutePowerSpectrum() cannot find window matrix");

  for(i=1;i<=32;i++) winFile >> k_conv[i];
  for(j=1;j<=100;j++) {
    winFile >> k_orig[j];
    for(i=1;i<=32;i++) winFile >> win[i][j];
  }
  winFile >> k_orig[j];
  for(i=1;i<=32;i++) winFile >> win[i][j];

  //  ofstream monitor("convolved2df.dat");

  /* output numerically calculated convolved P(k) */
  /* and convolved P(k) calculated using the window matrix */  
  for(j=1;j<=100;j++) {pk_orig[j] = TwoDF_Pk(k_orig[j], powerSpline); }
  for(i=1;i<=32;i++) {
    PkmcW = 0.0;
    for(j=1;j<=100;j++) PkmcW += win[i][j]*pk_orig[j];  // Matrix multiplication of convolution matrix with P(k)
    // printf("k=%g, P(k) conv W(k) = %g, %g\n",k_conv[i],PkconvW(k_conv[i], &powerSpline),PkmcW);
    convData = TwoDF_PkconvW(k_conv[i], powerSpline);
    //    monitor << k_conv[i] << "   " <<convData << endl;
    //cout << k_conv[i] << "   " << convData << endl;
    convTheoreticalData[i-1] = convData;
    //convTheoreticalData.push_back(convData);
    k_value[i-1] = k_conv[i];
    //k_value.push_back(k_conv[i]);
  }
  return TwoDF_chiSquared(convTheoreticalData, k_value);
}


/*! return normalised Power spectrum */
float AnalyzeThis::TwoDF_Pk(float k, Spline *spline)  
{  
  //float gam,pow_n;
  //float q, pn, tf, tfsq;
  //float tk_eh(float);
  
  if (k > spline->stop()) cout << "PK out of bound: " << k << endl;

  if(k< spline->start()) return 0.;
  if(k > spline->stop()) return 0;

  return (*spline)(k);
}

/*!  perform 3D convolution in k */
float AnalyzeThis::TwoDF_PkconvW(float rk, Spline *powerSpline)
{
  float s, sp, oa, a;
  int ik, imu, nint1=1, nint2=18;
  float dk, dmu, rkmax=0.3;
  float rmu, xk, xksq, rksq, yk;


  dmu  = 2.0/(float)nint2;
  rksq = rk*rk; 
  s = 0.0;
 
  /* first data point in trapezoidal rule */
  dk=0.5*rkmax;
  for(ik=1;ik<=2;ik++) {
    if(ik==1) { xk=xksq=0.0; } else { xk=rkmax; xksq=xk*xk; }
    for(sp=0.0,imu=1;imu<=nint2;imu++) {
      rmu = -1.0+((float)imu-0.5)*dmu;
      yk  = sqrt(rksq+xksq-2.*xk*rk*rmu);
      sp += TwoDF_Pk(yk, powerSpline);
    }
    s += sp*xksq*TwoDF_Wksq(xk);
  }
  a = s*dk;
  
  /* loop through next iterations until convergence */
  do {
    oa = a; s = 0.0;
    dk = rkmax/(float)nint1;
    for(ik=1;ik<=nint1;ik++) {
      xk   = ((float)ik-0.5)*dk;
      xksq = xk*xk;
      for(sp=0.0,imu=1;imu<=nint2;imu++) {
	rmu = -1.0+((float)imu-0.5)*dmu;
	yk  = sqrt(rksq+xksq-2.*xk*rk*rmu);
	sp  += TwoDF_Pk(yk, powerSpline);
      }
      s += sp*xksq*TwoDF_Wksq(xk);
    }
    a = 0.5*(oa+s*dk);
    nint1 *= 2;
  } while( fabs(oa-a)>1.0e-4*a );
  
  return a*dmu;
}

/*! return normalised window function (W_k)^2 */
float AnalyzeThis::TwoDF_Wksq(float k)
{
  float ksq, kp1, kp2, norm;
  ksq  = k*k;
  kp1  = 1.17011e-05; 
  kp2  = 9.33106e-09; 
  norm = 1150292.2;

  return norm/(1.+(ksq/kp1)+(ksq*ksq/kp2));
}

double AnalyzeThis::TwoDF_chiSquared(vector<double> convTheoreticalData, vector<double> k_theor) {
  //cout << "CHISQUARE" << endl;
  ifstream infile(ControlPanel::cmbeasyDir("/resources/2df/inverse_cov.dat").c_str());
  if (!infile) throw Bad_Error("AnalyzeThis::2dfchisquare(): Could not open inverse covariance matrix");
 
  int i,j;
  double sum, temp;
  double inverseCovariance[32][32];
  vector<double> dfData, dfk;
  
  // read the inverse covariance matrix
  for(i=0; i< 32; i++) {
    for(j=0; j< 32; j++) {
      if (infile.eof()) throw Bad_Error("scheiss");
      infile >> inverseCovariance[i][j];
      //cout << i << "  " << j << "   " <<  inverseCovariance[i][j] << endl;
    }
  }
  infile.close();

  ifstream infile2;

  // read the 2dfData
  infile2.open(ControlPanel::cmbeasyDir("/resources/2df/2dfdataShort.dat").c_str());
  if (!infile2) throw Bad_Error("AnalyzeThis::2dfchisquare(): Could not open 2dfDataFile");
     
  for(i=0; i<32; i++) {
    infile2 >> temp;
    //cout << i << "  " << temp;
    dfk.push_back(temp);
    infile2 >> temp;
    //cout << "   " << temp << endl;
    dfData.push_back(temp);
    // cout << dfk[i] << "  " << dfData[i] << "bla " <<  endl;
  }

  // consistency check: data at same k?
  for(i=0; i<32; i++ ) {
    if (k_theor[i]-dfk[i]/(k_theor[i]+dfk[i]) > 0.001)  
      throw exception();
  }
  
  sum=0; 
  // now, finally, we can calculate the chi squared thingie
  for(i=0; i < 32 ; i++ ) {
    for (j=0; j < 32 ; j++ ) {
      sum += (dfData[i]- convTheoreticalData[i])*inverseCovariance[i][j]*(dfData[j]- convTheoreticalData[j]);
    }
  }
  //cout << "Chi squared is: " << sum << endl;
  return sum;
}


/************************************************************************
 **                                                                    **
 **                 Baryon Acoustic Peak stuff                         **
 **                                                                    **
 ************************************************************************/
//Right now only the Eistenstein 05 A parameter is supported

/*!
  Baryon Acoustic Peak measurement from Eistenstein et al. 2005,
  ApJ, 633, 560.  This version is based on the A parameter,
  generalized for possible non-flatness and a different value for
  the scalar spectral index ns.

\param ns: scalar spectal index
\param Ode: early dark energy fraction (default = 0)

*/
double AnalyzeThis::SDSS_BAP_A_chiSquared(baseCosmos& cosmos, 
					  double ns, double Ode) const {

  const double zval = 0.35;  //Redshift from Eisenstein '05
  const double aval = 0.469; //Measured value of A from Eisenstein '05
  const double daval = 0.017; //Error in A
  const double ns_as = 0.98; //Value of n_s assumed by Eisenstein
  const double flattol = 1e-5; //If abs(omega_k) is less than this, assume flat
  double ciH0 = 2997.92458 / cosmos.h();  //c/H0 in Mpc units

  double effective_aval = aval * pow( ns / ns_as, -0.35 );  // A including n_s correction
  effective_aval *= sqrt(1 - Ode);  // Early dark energy correction

  //Avoid Universes that had no big bang
  if ( cosmos.isUniverseRebounding() )
    return numeric_limits<double>::infinity();

  //Calculate A value for current cosmological parameters
  double abs_ok = fabs( cosmos.omega_k() );   
  double Ez = ciH0 / cosmos.Z2iH( zval ); //E at zval (H/H0)
  double Apred; //Predicted value of A

  // Branch based on flatness.  Technically these equations are the
  // same, but the cancellation relies on limits as arguments go to
  // zero, which we can't necessarily rely on in the world of imperfect
  // numerical precision
  if ( abs_ok < flattol ) {
    //Assume flat
    Apred = sqrt( cosmos.omega_m() ) * pow( Ez, -1.0/3.0 ) *
      pow( 1.0/(zval*ciH0) * cosmos.propermotionDistance( zval ), 2.0/3.0 );
  } else {
    Apred = sqrt( cosmos.omega_m() ) * pow( abs_ok * Ez, -1.0/3.0 ) *
      pow( 1.0/(zval*ciH0) * cosmos.propermotionDistance( zval ), 2.0/3.0 );
  }

  return (effective_aval - Apred)*(effective_aval - Apred) / (daval * daval);
}

/*
 * Angular BAP distance for different dataset, the Baryon Acoustic Peak is 
 * detected at different angles and different redshifts. 
 * The actual value of the BAP depends on the fiducial model; along 
 * with every model we quote the corresponding value.
 * */

double AnalyzeThis::Angular_BAP_A_chiSquared(baseCosmos& cosmos, 
					  double ns, double Ode) const {

  const double zval = 0.35;  //Redshift from Eisenstein '05
  const double aval = 0.469; //Measured value of A from Eisenstein '05
  const double daval = 0.017; //Error in A
  const double ns_as = 0.98; //Value of n_s assumed by Eisenstein
  const double flattol = 1e-5; //If abs(omega_k) is less than this, assume flat
  double ciH0 = 2997.92458 / cosmos.h();  //c/H0 in Mpc units

  double effective_aval = aval * pow( ns / ns_as, -0.35 );  // A including n_s correction
  effective_aval *= sqrt(1 - Ode);  // Early dark energy correction

  //Avoid Universes that had no big bang
  if ( cosmos.isUniverseRebounding() )
    return numeric_limits<double>::infinity();

  //Calculate A value for current cosmological parameters
  double abs_ok = fabs( cosmos.omega_k() );   
  double Ez = ciH0 / cosmos.Z2iH( zval ); //E at zval (H/H0)
  double Apred; //Predicted value of A

  //Branch based on flatness.  Technically these equations are the
  // same, but the cancellation relies on limits as arguments go to
  // zero, which we can't necessarily rely on in the world of imperfect
  // numerical precision
  if ( abs_ok < flattol ) {
    //Assume flat
    Apred = sqrt( cosmos.omega_m() ) * pow( Ez, -1.0/3.0 ) *
      pow( 1.0/(zval*ciH0) * cosmos.propermotionDistance( zval ), 2.0/3.0 );
  } else {
    Apred = sqrt( cosmos.omega_m() ) * pow( abs_ok * Ez, -1.0/3.0 ) *
      pow( 1.0/(zval*ciH0) * cosmos.propermotionDistance( zval ), 2.0/3.0 );
  }

  return (effective_aval - Apred)*(effective_aval - Apred) / (daval * daval);
}


/* ************************************************************************************************************************************
**
**
**  VSA, ACBAR, BOOMERANG and CBI stuff follows
**
**
**
**
**
************************************************************************************************************************************** */


/*!
  BOOMERANG chi^2
  written by Christian M. Mueller,
  modifications by Georg Robbers

  This is the initilization routine, needs to
  be called before any of the chi^2 methods.
   \param cl The Cl's in units of muK^2.
   Note that  you need the Cl-spectrum up to 2000.
*/
void AnalyzeThis::BOOMERANGInit( CL &cl ){
  const int bands=47;
  Spline ClTT(*cl.ts[0],"BOOMERANGwithcalibration_ts");
  Spline ClTTtensor(*cl.tt[0],"BOOMERANGwithcalibration_tttensor");
  Spline ClTE(*cl.cs[0],"BOOMERANGwithcalibration_tte");
  Spline ClEE(*cl.es[0],"BOOMERANGwithcalibration_tee");
  Spline ClBB(*cl.bt[0],"BOOMERANGwithcalibration_tbb");
  ClTT.arm();
  ClTE.arm();
  ClEE.arm();
  ClBB.arm();

  map<int,bool> &deselect= BOOMERANGBandDeselect;

  std::fill_n( BOOMERANGthband, bands, 0. );

  const string DataFileName="B03_NA_21July05.newdat";
  NewdatReader boomerang((ControlPanel::cmbeasyDir("/resources/boomerang/")+DataFileName).c_str());
  BOOMERANGexpData = boomerang.allBands();


  for(int i=0; i<bands; ++i){
    if (deselect.find(i+1) == deselect.end() ){
      BOOMERANGthband[i] = BOOMERANGWindowConv(i+1,ClTT, ClTE, ClEE, ClBB );
      //cout << "Thband " << thband[i] << endl;
    }
  }

  ClTT.disarm();
  ClTE.disarm();
  ClEE.disarm();
  ClBB.disarm();


// The following block reads the inverse Fisher matrix directly from the
// data file. For speed reasons, and because we don't want to link to the GSL
// here, we read it from a binary file that was constructed using the code below
// and is shipped with CMBEASY.
//X   NewdatReader::Matrix invFisherMatrix = boomerang.inverseFisherMatrix();
//X   //The Fisher matrix of the bandpowers is transformed as...(see paper)
//X   for(int i=0; i<bands; i++)
//X       for(int j=0; j < bands; j++)
//X       {
//X         if ( i < boomerang.nrOfTTBands() )
//X           invFisherMatrix[i][j] /= ( BOOMERANGexpData[i].power+BOOMERANGexpData[i].offset );
//X 
//X         if ( j < boomerang.nrOfTTBands() )
//X           invFisherMatrix[i][j] /= ( BOOMERANGexpData[j].power+BOOMERANGexpData[j].offset );
//X       }
//X 
//X   //cout << "invFisherafterdividing: " << invFisherMatrix[ 0 ][ 0 ] << endl;
//X   gsl_permutation* perm = gsl_permutation_alloc( bands );
//X   gsl_matrix* inverse = gsl_matrix_alloc( bands, bands );
//X   gsl_matrix* fisher = gsl_matrix_alloc( bands, bands );
//X 
//X   for ( int i = 0; i < bands; ++i )
//X     for ( int j = 0; j < bands; ++j )
//X       gsl_matrix_set( inverse, i, j, invFisherMatrix[i][j] );
//X 
//X 
//X   int status;
//X   gsl_linalg_LU_decomp( inverse, perm, &status );
//X   gsl_linalg_LU_invert( inverse, perm, fisher );
//X 
//X   BOOMERANGFisherMatrix=NewdatReader::Matrix( bands );
//X   for ( int i = 0; i < bands; ++i)
//X     BOOMERANGFisherMatrix[ i ].resize( bands );
//X 
//X 
//X   for ( int i = 0; i < bands; ++i )
//X     for ( int j = 0; j < bands; ++j )
//X       BOOMERANGFisherMatrix[i][j] = gsl_matrix_get( fisher, i, j );
//X 
//X   string fisherFileName = ControlPanel::cmbeasyDir("/resources/boomerang/FisherMatrix.dat");
//X   ofstream out( fisherFileName.c_str() );
//X   for ( int i = 0; i < bands; ++i )
//X   {
//X     for ( int j = 0; j < bands; ++j )
//X       out << std::scientific << BOOMERANGFisherMatrix[i][j] << " ";
//X     out << endl;
//X   }
//X   out.close();

  BOOMERANGFisherMatrix=NewdatReader::Matrix( bands );
  for ( int i = 0; i < bands; ++i)
    BOOMERANGFisherMatrix[ i ].resize( bands );

  string fisherFileName = ControlPanel::cmbeasyDir("/resources/boomerang/FisherMatrix.dat");
  ifstream in( fisherFileName.c_str() );
  for ( int i = 0; i < bands; ++i )
    for ( int j = 0; j < bands; ++j )
    {
      in >> std::scientific >> BOOMERANGFisherMatrix[i][j];
    }
  in.close();
}

/*!
   BOOMERANGChi2 returns the chi^2 w.r.t. BOOMERANG03, optionally taking
   a beam error and/or calibration error into account. If you want the marginalized
   version ( over beam error and calibration uncertainty ) of chi^2, use 
   BOOMERANGChi2WithCalibration() instead.
   \warning call BOOMERANGInit() first!
   \param beamErrorAdjust add beamErrorAdjust * beamError(i) to
                          theoretical band power(i) before calculating chi^2
   \param calibrationFactor multiply theoretical bandpowers by this factor before comparing
*/
double AnalyzeThis::BOOMERANGChi2( double beamErrorAdjust, double calibrationFactor ){
  const int bands = 47;
  double sum=0;
  double ClTheoretical[bands];
  double ClExperimental[bands];

  map<int, bool>& deselect = BOOMERANGBandDeselect;

  double thband[bands];
  std::fill_n( ClTheoretical, bands, 0. );
  std::fill_n( ClExperimental, bands, 0. );
  std::copy( BOOMERANGthband, BOOMERANGthband+bands, thband );

  for(int i=0; i<bands; ++i){
    if (deselect.find(i+1) == deselect.end() ){
      if ( beamErrorAdjust != 0 )
        thband[i] += thband[i]*beamErrorAdjust*BOOMERANGexpData[i].beamError;
      thband[i] *= calibrationFactor;
    }
  }
  for(int i=0; i<bands; i++){
    if (deselect.find(i+1) == deselect.end()){
      if (BOOMERANGexpData[i].offsetFlag){
        ClTheoretical[i]=log(thband[i]+BOOMERANGexpData[i].offset);
        ClExperimental[i]=log(BOOMERANGexpData[i].power+BOOMERANGexpData[i].offset);
      }
      else {
        ClTheoretical[i]=thband[i];
        ClExperimental[i]=BOOMERANGexpData[i].power;
      }
    }
  }
  //copy( ClTheoretical, ClTheoretical + bands, ostream_iterator<double>(cout, "\n") );

  for(int i=0; i<bands; i++){
    if ( deselect.find(i+1) == deselect.end()) // if the bands have not been de-selected
      for(int j=0; j < bands; j++){
        if ( deselect.find(j+1) == deselect.end()) {
          double contrib =(-ClExperimental[i]+ClTheoretical[i])*BOOMERANGFisherMatrix[i][j]*(-ClExperimental[j]+ClTheoretical[j]);
          sum += contrib;
        }
      }
  }
  return sum;
}

/*!
  Return BOOMERANG chi^2.
  BOOMERANGChi2WithCalibration() returns the chi^2 with regard to the BOOMERANG03 data
  ( NA pipeline ) marginalized over the beam error and the calibration uncertainty of +- 4%.
  \warning Call BOOMERANGInit() first!
  \returns the marginalized chi^2 for BOOMERANG03
*/
double AnalyzeThis::BOOMERANGChi2WithCalibration() {
  //marginalize over beam and calibration errors, assuming gaussian error
  vector<double> likeCalib, likeBeam, weights;
  const double steps = 5;
  //allow for +/- 3 sigma
  for ( double beam = -steps; beam <= steps; ++beam )
    weights.push_back( exp( -pow( beam*3./steps, 2. ) /2. ) );
  double norm = std::accumulate( weights.begin(), weights.end(), 0. );

  //for (  double beam = -3.0; beam <= 3.0; beam += 0.5 ){
  for ( double beam = -steps; beam <= steps; ++beam ){
    for (double f = -steps; f <= steps; ++f) {
      //0.04 is BOOMERANG03 calibration uncertainty
      double boomerang = BOOMERANGChi2(beam*3/steps,1+0.04*f*3./steps);
      likeCalib.push_back( boomerang );
      //cout << f << " calibration factor and " <<(  beam*3./steps ) << " beam Error gives chi^2: " << boomerang << endl;
    }
    double calibMin = *std::min_element( likeCalib.begin(), likeCalib.end() ) ;
    for ( unsigned int i = 0; i < likeCalib.size(); ++i )
      likeCalib[i] = exp( -( likeCalib[i]-calibMin )/2. )*weights[i];

    likeBeam.push_back( -2. * log( accumulate( likeCalib.begin(), likeCalib.end(), 0. )/norm ) + calibMin );
    likeCalib.clear();
  }
  double beamMin = *std::min_element( likeBeam.begin(), likeBeam.end() ) ;
  for ( unsigned int i = 0; i < likeBeam.size(); ++i )
    likeBeam[i] = exp( -( likeBeam[i]-beamMin )/2. )*weights[i];

  //cout << "boomerang min chi^2 says: " << beamMin << endl;
  double res = -2. * log( accumulate( likeBeam.begin(), likeBeam.end(), 0. )/norm ) + beamMin;
  cout << "boomerang like says: " << res << endl;
  return res;
}



/*!Function for convoluting theoretical Cl-spectra with BOOMERANG03 window functions
  Written by Christian M. Mueller
 */
double AnalyzeThis::BOOMERANGWindowConv(int windowNumber,Spline &ClTT, Spline &ClTE, Spline &ClEE, Spline &ClBB){
  // still true though: cout << "WARNING: AnalyzeThis::BOOMERANGWindowConv(): implement Tensors and BB scalar!" << endl;

  //open WindowFile
  string windowName="B03_NA_21July05_";
  stringstream s;
  s << windowNumber;
  string windowFileName = ControlPanel::cmbeasyDir("/resources/boomerang/windows/")+windowName+s.str();
  ifstream windowFile(windowFileName.c_str());
  if (!windowFile) cout << "Error in BOOMERANGWindowConv: Could not open " << windowFileName << endl;

  //read all data from this file
  // the number of bins in the associated window functions
  int bins=2000;
  // note that the l value increases by one for every binning.; hence we only need the window function and
  //the beginning l value
  // i is just a dummy variable
  int i=0;
  int begin;
  vector<double> TTwindow(bins+1);
  vector<double> TEwindow(bins+1);
  vector<double> EEwindow(bins+1);
  vector<double> BBwindow(bins+1);

  windowFile >> begin >> TTwindow[2] >> TEwindow[2] >> EEwindow[2] >> BBwindow[2] ;
  TTwindow[2] *= 2; TEwindow[2] *= 2; EEwindow[2] *= 2; BBwindow[2] *=2; // Windows in file are W_l/l
  //now we can read the rest of the values
  for (int j=3; j < bins+1; j++){
    windowFile >> i >> TTwindow[j] >> TEwindow[j] >> EEwindow[j] >> BBwindow[j];
    TTwindow[j] *= j; TEwindow[j] *= j; EEwindow[j] *= j; BBwindow[j] *=j; // Windows in file are W_l/l
  }

  // to save time, we only multiply with the relevant entries
  int TT=24;
  //int EE=31;
  int BB=38;
  int TE=47;

  //now, multiply theoretical data with the window function 

  double sum=0;
  double norm=0;
  //  cout << windowNumber << endl;
  //cout << "Warning: only TT" << endl;
  //if ( windowNumber > TT ) return 0;
  for (int j=2; j<bins+1; j++) {
    if(windowNumber <= TT) {
      double bareCj = ClTT(j) /(j*(j+1))*2.*M_PI;
      sum+=TTwindow[j]*( j+0.5 )*bareCj/( 2.*M_PI );
      norm+=TTwindow[j]*(j+0.5)/(j*(j+1));
    }
    else if (TT < windowNumber && windowNumber <= BB ){
      double bareCj = ClEE(j) /(j*(j+1))*2.*M_PI;
      sum+=EEwindow[j]*( j+0.5 )*bareCj/( 2.*M_PI );
      bareCj = ClBB(j) /(j*(j+1))*2.*M_PI;
      sum+=BBwindow[j]*( j+0.5 )*ClBB(j)/( 2.*M_PI );
      norm+=(EEwindow[j]+BBwindow[j])*(j+0.5)/(j*(j+1));
    }
    else if ( BB < windowNumber && windowNumber <= TE){
      double bareCj = ClTE(j) /(j*(j+1))*2.*M_PI;
      sum+=TEwindow[j]*( j+0.5 )*bareCj/( 2.*M_PI );
      norm+=TEwindow[j]*(j+0.5)/(j*(j+1));
      // cout << "sum " << sum << " norm " << norm << endl;
    }
  }
  // cout << "sum: " << sum << " norm " << norm <<  endl;  // Windows are indeed (almost) normalized for Boomerang
  return sum/norm;
}

void AnalyzeThis::BOOMERANGDeselect(bool TT, bool TE, bool EE, bool BB){
  if (TT)  for (int i = 1; i <= 24; i ++) BOOMERANGBandDeselect[i] = true;
  if (EE)  for (int i=25; i<=31; i++)  BOOMERANGBandDeselect[i] = true;
  if (BB)  for (int i=32; i<=38; i++)  BOOMERANGBandDeselect[i] = true;
  if (TE) for (int i=39; i<=47; i++)  BOOMERANGBandDeselect[i] = true;
}

/*! BOOMERANGDeselect designates  particular bands to be excluded from the analysis.
    Useful for  combining the B03 dataset with e.g. WMAP data.
    McTavish et. al( astro-ph/0507503 ) say in section 3.1:

    \verbatim
    The cosmic variance of the WMAP and B03 data sets is correlated in the low multipole range( ... ).
    To account for this, we cut the lower multipoles of the B03 TT spectrum ( l < 375 ) when
    combining the B03 data with WMAP data.
    \endverbatim

    Cutting l < 375 corresponds to deselecting bands 1 to 6:

    \code
    AnalyzeThis ai;
    for ( int i = 1; i <= 6; ++i )
      ai.BOOMERANGDeselect( i );
    \endcode
*/
void AnalyzeThis::BOOMERANGDeselect(int band){
  BOOMERANGBandDeselect[band] = true;
}
// ------------------------ END OF BOOMERANG -------------------




/*!Function for convoluting theoretical Cl-spectra with VSA window functions
  Written by Christian M. Mueller
 */

double AnalyzeThis::VSAWindowConv(int windowNumber, Spline &Cl, double cut_l){
  //open WindowFile
  string windowName[16]={"VSAF1","VSAF2","VSAF3","VSAF4","VSAF5","VSAF6","VSAF7","VSAF8","VSAF9","VSAF10","VSAF11","VSAF12",
			   "VSAF13","VSAF14","VSAF15","VSAF16"};
  string windowFileName = ControlPanel::cmbeasyDir("/resources/vsa/")+windowName[windowNumber-1];
  ifstream windowFile(windowFileName.c_str());
  if (!windowFile) cout << "Error in VSAWindowConv: Could not open " << windowFileName << endl;
  //read all data from this file
  // the number of bins in the associated window functions
  int binnumbers[17]={0,194,199,377,374,396,362,396,380,374,371,404,404,397,434,434,531} ;
  int bins=binnumbers[windowNumber];
  // note that the l value increases by one for every binning.; hence we only need the window function and
  //the beginning l value
  // i is just a dummy variable
  int i=0;
  int begin;
  vector<double> window(bins);
  windowFile >> begin >> window[0];
  //now we can read the rest of the values
  for (int j=1; j < bins; j++){
    windowFile >> i >> window[j];
  }
 

  //now, multiply theoretical data with the window function 
  
  double sum=0;
  double normalize = 1.0; // pow(2.725e6,-2); // T_cmb ^ {-2}
  for (int j=0; j<bins; j++) {
    //  cout << "Binning now: " << begin+j << endl;
    if (begin+j < cut_l) sum+=window[j]*Cl(begin+j)*normalize;
  }
  // cout << "sum: " << sum << endl;
  return sum;
}

/*!
  VSA chi^2
  written by Christian M. Mueller
   \param Cl The Cl's in units of muK^2 
   Note that if using all the data points you need the Cl-spectrum up to 1725
*/
double AnalyzeThis::VSAChi2(Spline &Cl){
//  string band[16];
  double dummy;
  double FischerMatrix[16][16];
  double offset[16];
  double expData[16];
  double sum=0;
  string DataFileName="vsa.dat";
  string gaussianMatrixFileName="VSAF_mat.dat";
  string offsetFileName="VSAF_x.dat";
 
  // find windowed theoretical values
  double ClTheoretical[16];

  // Cl spectrum is needed up to 1725 if using all points;
  
  // we cut at 300 l's below the splines end. For the full and correct
  // vsa, use jlgen to generate jl's up to 1725 (internally, they will then
  // be up to 2025, and after subracting 300 here, we are fine)
  // on the other hand, you may not want to include all VSA and then
  // there is no point in going up to l = 1725, cause the window function
  // will attach no weight to high l for a low  l bin. For instance, if you omit
  // the last data point, you just need l=1500. 
  for (int i=0; i <16 ; i++) {
    ClTheoretical[i]= VSAWindowConv(i+1,Cl,Cl.stop()-300);
  }

  //read experimental values     
  ifstream dataFile((ControlPanel::cmbeasyDir("/resources/vsa/")+DataFileName).c_str());
  if (!dataFile) throw Bad_Error("AnalyzeThis::VSAChi2(): dataFile not found");
  for (int i=0; i<16; i++) {
    dataFile >> expData[i] >> dummy >> dummy;
    //    cout << "Experimental: " << expData[i] << endl;
  }
  dataFile.close();

  //read gaussianFischerMatrix
  ifstream matrixFile((ControlPanel::cmbeasyDir("/resources/vsa/")+gaussianMatrixFileName).c_str());
  if (!matrixFile) throw Bad_Error("AnalyzeThis::VSAChi2(): matrixFile not found");
  for (int i=0; i<16; i++) {
    for (int j=0; j<16; j++) {
      matrixFile >> FischerMatrix[i][j];
    }
  }

  dummy=0;
  
  // read offsets
  ifstream offsetFile((ControlPanel::cmbeasyDir("/resources/vsa/")+offsetFileName).c_str());
  if (!offsetFile) throw Bad_Error("AnalyzeThis::VSAChi2(): offsetFile not found");
  for (int i=0; i<16; i++) {
    offsetFile>> offset[i];
    //    cout << "offset: " << offset[i] << endl;
  }
  
  map<int,bool> &deselect = VSABandDeselect;//!< All integers appearing here will be considered de-selected from the VSA likelihood. Counting goes from 1 to 16.
// compute chi squared
  for(int i=0; i<16; i++){
    for(int j=0; j < 16; j++){
      // if the bands have not been de-selected
      if (deselect.find(i+1) == deselect.end() && deselect.find(j+1) == deselect.end()) { 
	double contrib =(log(expData[i]+offset[i])-log(ClTheoretical[i]+offset[i]))*FischerMatrix[i][j]*(log(expData[j]+offset[j])-log(ClTheoretical[j]+offset[j]));
	//cout << "adding: " << i << " : " << j << "  Fisher: " << FischerMatrix[i][j] << "  contibuing: " << contrib << endl;
	sum += contrib;
      }
    }
  }

return sum;
}

/*!
  Return VSA chi^2 using the calibration uncertainty of +- 3%
  \param Cl The Cl's in units of muK^2 
*/
double AnalyzeThis::VSAChi2WithCalibration(Spline& Cl) {  
  Spline like(100, "like");
  Spline ts(Cl,"vsawithcalibration_ts");

  for (double f = 0.97; f <= 1.035; f += 0.01) {
    ts *= f;
    ts.arm();
    double vsa = -0.5*VSAChi2(ts);
    like.set(f,vsa);
    ts.disarm();
    ts *= 1.0/f;
  }
  like.arm();
  // like.dump("vsalike");
  double x = like.maximum();
  cout << "vsa like says: " << x << endl;
  
  ts *= x;
  ts.arm();
  double res =  VSAChi2(ts);
  return res;
}

// CBI 2000+2001 follows

/*!Function for convoluting theoretical Cl-spectra with CBI00+01 window functions
  Written by Christian M. Mueller
 */

double AnalyzeThis::CBIWindowConv(int windowNumber, Spline &Cl, double cut_l){
  //open WindowFile
  string windowName[13]={"CBI2.0_2to141","CBI2.0_2to142","CBI2.0_2to143","CBI2.0_2to144","CBI2.0_2to145","CBI2.0_2to146",
	   "CBI2.0_2to147","CBI2.0_2to148","CBI2.0_2to149","CBI2.0_2to1410","CBI2.0_2to1411","CBI2.0_2to1412","CBI2.0_2to1413"};
  string windowFileName = ControlPanel::cmbeasyDir("/resources/cbi2/")+windowName[windowNumber-1];
  ifstream windowFile(windowFileName.c_str());
  if (!windowFile) cout << "Error in CBIWindowConv: Could not open " << windowFileName << endl;
  //read all data from this file
  // the number of bins in the associated window functions (for CBI the same in every window function file)
  int bins=3499;
  // note that the l value increases by one for every binning.; hence we only need the window function and
  //the beginning l value
  // i is just a dummy variable
  int i=0;
  int begin;
  vector<double> window(bins);
  windowFile >> begin >> window[0];
  //now we can read the rest of the values
  for (int j=1; j < bins; j++){
    windowFile >> i >> window[j];
  }
 

  //now, multiply theoretical data with the window function 
  // note that the window functions are not normalized, hence
  // we have to do that (see CBI webpage for details)
  
  double sum=0;
  double windowNorm=0;
  for (int j=0; j<bins; j++) {
    //  cout << "Binning now: " << begin+j << endl;
    if (begin+j < cut_l) {
      sum += window[j]*Cl(begin+j)/(begin+j);
      windowNorm += window[j]/(begin+j);
    }
  }
  return sum/windowNorm;
}

/*!
  CBI chi^2
  written by Christian M. Mueller
   \param Cl The Cl's in units of muK^2 
   Note that if using all the data points you need the Cl-spectrum up to 3500;
   consult astro-ph/0402359 for details and which points to choose.
*/
double AnalyzeThis::CBIChi2(Spline &Cl){ 
  //string band[13];
  double dummy;
  double FischerMatrix[13][13];
  double offset[13];
  double expData[13];
  double sum=0;
  string DataFileName="cbi.dat";
  string gaussianMatrixFileName="CBI2.0_2to14_matrix.dat";
  string offsetFileName="CBI2.0_2to14_x.dat";
 
  // find windowed theoretical values
  double ClTheoretical[13];

  // Cl spectrum is needed up to 3500 if using all points;
  
  // we cut at 300 l's below the splines end. For the full and correct
  // cbi, use jlgen to generate jl's up to 3500 (internally, they will then
  // be up to 3800, and after subracting 300 here, we are fine)
  // on the other hand, you may not want to include all CBI and then
  // there is no point in going up to l = 3500, cause the window function
  // will attach no weight to high l for a low  l bin. For instance, if you omit
  // the last two data point, you just need l=1800. 
  for (int i=0; i <13 ; i++) {  
    ClTheoretical[i]= CBIWindowConv(i+1,Cl,Cl.stop()-300);
    //  cout << ClTheoretical[i] << endl;
  }

  //read experimental values     
  ifstream dataFile((ControlPanel::cmbeasyDir("/resources/cbi2/")+DataFileName).c_str());
  if (!dataFile) throw Bad_Error("AnalyzeThis::CBIChi2(): dataFile not found");
  for (int i=0; i<13; i++) {
    dataFile >> expData[i] >> dummy >> dummy;
    //   cout << "Experimental: " << expData[i] << endl;
  }
  dataFile.close();

  //read gaussianFischerMatrix
  ifstream matrixFile((ControlPanel::cmbeasyDir("/resources/cbi2/")+gaussianMatrixFileName).c_str());
  if (!matrixFile) throw Bad_Error("AnalyzeThis::CBIChi2(): matrixFile not found");
  for (int i=0; i<13; i++) {
    for (int j=0; j<13; j++) {
      matrixFile >> FischerMatrix[i][j];
    }
  }

  dummy=0;
  
  // read offsets
  ifstream offsetFile((ControlPanel::cmbeasyDir("/resources/cbi2/")+offsetFileName).c_str());
  if (!offsetFile) throw Bad_Error("AnalyzeThis::CBIChi2(): offsetFile not found");
  for (int i=0; i<13; i++) {
    offsetFile>> offset[i] ;
  }
  
  map<int,bool> &deselect = cbiBandDeselect;//!< All integers appearing here will be considered de-selected from the CBI likelihood. Counting goes from 1 to 13.
// compute chi squared
  for(int i=0; i<13; i++){
    for(int j=0; j < 13; j++){
      // if the bands have not been de-selected
      if (deselect.find(i+1) == deselect.end() && deselect.find(j+1) == deselect.end()) { 
	double contrib =(log(expData[i]+offset[i])-log(ClTheoretical[i]+offset[i]))*FischerMatrix[i][j]*
	  (log(expData[j]+offset[j])-log(ClTheoretical[j]+offset[j]));
	sum += contrib;
	//	cout << "contributing exp " << expData[i]+offset[i]  << " for band " << i << " for theory: " <<ClTheoretical[i]+offset[i] << endl; 
      }
    }
  }
return sum;
}

/*!
  Return CBI chi^2 using the calibration uncertainty of +- 2.6%
  \param Cl The Cl's in units of muK^2 
*/
double AnalyzeThis::CBIChi2WithCalibration(Spline& Cl) {  
  Spline like(100, "cbilike");
  Spline ts(Cl,"cbiwithcalibration_ts");

  for (double f = 0.976; f <= 1.027; f += 0.01) {
    ts *= f;
    ts.arm();
    double cbi = -0.5*CBIChi2(ts);
    like.set(f,cbi);
    ts.disarm();
    ts *= 1.0/f;
  }
  like.arm();
  // like.dump("cbilike");
  double x = like.maximum();
  //  cout << "cbi like says: " << x << endl;
  ts *= x;
  ts.arm();
  double res =  CBIChi2(ts);
  return res;
}

/*! Function for convoluting theoretical Cl-spectra with the respective 
window functions of the experiment. For CBI DEEP and MOSAIC only.
 DEEP has binnumber 50 and binwidth 90, MOSAIC has binnumber 128 and binwidth 25
 Written by Christian M. Mueller
*/

double AnalyzeThis::WindowConv(string windowFileName, Spline &Cl, int binnumber, int binwidth, double cut_l) {
  // open WindowFile
  // cout << "BAND: " << windowFileName << " and cut at: " << cut_l << endl;
  windowFileName = ControlPanel::cmbeasyDir("/resources/cbi/") + windowFileName + ".dat";
  ifstream windowFile(windowFileName.c_str());
  if (!windowFile) cout << "Could not open " << windowFileName << endl;
  
  // read all data from this file
  // i is just a dummy for the uneeded binning value in the file
  int i=0;
  vector<double> window(binnumber);
  for(int j=0; j<binnumber;j++) {
    windowFile >> i >> window[j];
  }			
  
 
  //now, multiply theoretical data with the window function according to 
  // Sum W(l)/l l(l+1)Cl/2Pi to find the averaged bandpower
 
  double sum=0;
  
  // typedef map<int,double>::const_iterator IT;
  double normalize = pow(2.725e6,-2); // T_cmb ^ {-2}
  for (int j=0; j<binnumber; j++) {
    for(int k=1; k<=binwidth/2; k++) {
      double l = 2*k+binwidth*j;
      //cout << "window for l: " << l << "   is: " << window[j] << endl;
      if (l <= cut_l) 
	sum+=window[j]/(binwidth/2)*Cl(2*k+binwidth*j) * normalize;    // only even  l is present, i.e. 2,4,6..
    }
  }
  //cout << "sum: " << sum << endl;
  return sum;
}

/*!
 Same as WindowConv, only this time for ACBAR
  written by Christian M. Mueller 
*/

double AnalyzeThis::WindowConvACBAR( string windowFileName, Spline &Cl,double cut_l) {
   // open WindowFile
  ifstream windowFile(ControlPanel::cmbeasyDir("/resources/acbar/"+windowFileName).c_str());
  if (!windowFile) cout << "Could not open " << windowFileName << endl;
  //cout << "ACBAR band: " << windowFileName <<endl;

  // read all data from this file
  // i is just a dummy for the uneeded binning value in the file
  int i=0;
  vector<double> window(100);

  for(int j=0; j<100;j++) {
    windowFile >> i >> i >> window[j];
  }			
  
 
  //now, multiply theoretical data with the window function according to 
  // Sum W(l)/l l(l+1)Cl/2Pi to find the averaged bandpower
 
  double sum=0;
  
  // typedef map<int,double>::const_iterator IT;

  for (int j=0; j<100; j++) {
    for(int k=1; k<=15; k++) {
      double l = 2*k + 30*j + 15 ;
    
      if (l <= cut_l) 
	sum+=window[j]*Cl(2*k+30*j+15)*2;    // only even  l is present, i.e. 2,4,6..
      //cout << 2*k+90*j << " " << Cl[2*k+90*j] << endl;
      // cout << sum << endl;
      }
  }
  // cout << "Sum: " <<  sum*2.725*2.725*pow(10,12) << endl;
  
  return sum;
}

/*!
  Return CBI Mosaic chi^2 using the calibration uncertainty of +- 10%
  \param Cl The Cl's in units of muK^2 
*/
double AnalyzeThis::CBIMosaicChi2WithCalibration(Spline& Cl) {
  Spline like(100, "like");
  Spline ts(Cl,"cbimosaicwithcalibration_ts");

  for (double f = 0.9; f <= 1.101; f += 0.02) {
    ts *= f;
    ts.arm();
    double cbi = -0.5*CBIMosaicChi2(ts);
    like.set(f,cbi);
    ts.disarm();
    ts *= 1.0/f;
  }
  like.arm();
  like.dump("cbilike");
  double x = like.maximum();
  //  cout << "cbi like says: " << x << endl;
  
  ts *= x;
  ts.arm();
  double res =  CBIMosaicChi2(ts);
  return res;
}

double AnalyzeThis::ACBARChi2WithCalibration(Spline& Cl) {
  
  Spline like(100, "like");
  Spline ts(Cl,"acbarwithcalibration_ts");
  
  //ts.dump("acbardump");

  for (double f = 0.8; f <= 1.201; f += 0.02) {
    ts *= f;
    ts.arm();
    double acbar = -0.5*ACBARChi2(ts);
    //cout << "f:" << f<< "acbar like: " << acbar << endl;
    like.set(f,acbar);
    ts.disarm();
    ts *= 1.0/f;
  }
  like.arm();
  //like.dump("cbilike");
  double x = like.maximum();
  //  cout << "ACBAR like says: " << x << endl;
  
  ts *= x;
  ts.arm();
  double res =  ACBARChi2(ts);
  return res;
}
/*!
  CBI - MOSAIC chi^2
  written by Christian M. Mueller
   \param Cl The Cl's in units of muK^2 
*/
double AnalyzeThis::CBIMosaicChi2(Spline &Cl) {
  string band[14];
  double dummy;
  double FischerMatrix[14][14];
  double offset[14];
  double expData[14];
  double sum=0;
  string DataFileName="joint_final_iso_0.08_200_even_qb";
  string gaussianMatrixFileName="GaussianFischerMatrixMosaic";
  string offsetFileName="joint_final_iso_0.08_200_even_otherps";
 
  band[0]="joint_window_iso_0.08_200_even_deep_3level_feb_window_1";
  band[1]="joint_window_iso_0.08_200_even_deep_3level_feb_window_2";
  band[2]="joint_window_iso_0.08_200_even_deep_3level_feb_window_3";
  band[3]="joint_window_iso_0.08_200_even_deep_3level_feb_window_4";
  band[4]="joint_window_iso_0.08_200_even_deep_3level_feb_window_5";
  band[5]="joint_window_iso_0.08_200_even_deep_3level_feb_window_6";
  band[6]="joint_window_iso_0.08_200_even_deep_3level_feb_window_7";
  band[7]="joint_window_iso_0.08_200_even_deep_3level_feb_window_8";
  band[8]="joint_window_iso_0.08_200_even_deep_3level_feb_window_9";
  band[9]="joint_window_iso_0.08_200_even_deep_3level_feb_window_10";
  band[10]="joint_window_iso_0.08_200_even_deep_3level_feb_window_11";
  band[11]="joint_window_iso_0.08_200_even_deep_3level_feb_window_12";
  band[12]="joint_window_iso_0.08_200_even_deep_3level_feb_window_13";
  band[13]="joint_window_iso_0.08_200_even_deep_3level_feb_window_14";

  // find windowed theoretical values
  double ClTheoretical[14];

  // we cut at 300 l's below the splines end. For the full and correct
  // cbi mosaic, use jlgen to generate jl's up to 4500 (internally, they will then
  // be up to 4800, and after subracting 300 here, we are fine)
  // on the other hand, you may not want to include all cbi-points and then
  // there is no point in going up to l = 4500, cause the window function
  // will attach no weight to high l for a low  l bin.
  for (int i=0; i <14 ; i++) {
    ClTheoretical[i]=WindowConv(band[i],Cl,128,25, Cl.stop()-300);
  }

  //read experimental values                                                     joint_final_0.08_200_even_qb.dat
  ifstream dataFile(ControlPanel::cmbeasyDir("/resources/cbi/joint_final_iso_0.08_200_even_qb.dat").c_str());
  if (!dataFile) throw Bad_Error("AnalyzeThis::cbiMosaicLikelihood(): dataFile not found");
  for (int i=0; i<14; i++) {
    dataFile >> dummy >> expData[i] >> dummy;
    //cout << expData[i];
  }
  dataFile.close();

  //read gaussianFischerMatrix
  ifstream matrixFile(ControlPanel::cmbeasyDir("/resources/cbi/GaussianFischerMatrixMosaic.dat").c_str());
  if (!matrixFile) throw Bad_Error("AnalyzeThis::cbiMosaicLikelihood(): matrixFile not found");
  for (int i=0; i<14; i++) {
    for (int j=0; j<14; j++) {
      matrixFile >> FischerMatrix[i][j];
    }
  }

  dummy=0;
  
  // read offsets
  ifstream offsetFile(ControlPanel::cmbeasyDir("/resources/cbi/joint_final_iso_0.08_200_even_otherps.dat").c_str());
  if (!offsetFile) throw Bad_Error("AnalyzeThis::cbiMosaicLikelihood(): offsetFile not found");
  for (int i=0; i<14; i++) {
    offset[i]=0;
    for ( int j=0; j<4; j++ ){ 
      offsetFile >> dummy;
      offset[i]+=dummy;
    }
  }
  
// compute chi squared
  map<int,bool> &deselect = cbiMosaicBandDeselect;
  for(int i=0; i<14; i++){
    for(int j=0; j < 14; j++){
      // if the bands have not been de-selected.
      if (deselect.find(i+1) == deselect.end() && deselect.find(j+1) == deselect.end()) { 

	double contrib =(log(expData[i]+offset[i])-log(ClTheoretical[i]+offset[i]))*FischerMatrix[i][j]*(log(expData[j]+offset[j])-log(ClTheoretical[j]+offset[j]));
	//cout << "adding: " << i << " : " << j << "  Fisher: " << FischerMatrix[i][j] << "  contibuing: " << contrib << endl;
	sum += contrib;
      }
    }
  }
  
  return sum;
}

/*!
  This will give the CBI DEEP likelihood
  written by Christian M. Mueller
   \param Cl The Cl's in units of muK^2 
*/

double AnalyzeThis::CBIDeepChi2(Spline &Cl){ 
  // name of window functions
  string band[6];
  double dummy;
  double FischerMatrix[6][6];
  double offset[6];
  double expData[6];
  double sum=0;
 
  band[0]="joint_3deep_std2_best_window_1";
  band[1]="joint_3deep_std2_best_window_2";
  band[2]="joint_3deep_std2_best_window_3";
  band[3]="joint_3deep_std2_best_window_4";
  band[4]="joint_3deep_std2_best_window_5";
  band[5]="joint_3deep_std2_best_window_6";
  

  
  // find windowed theoretical values
  double ClTheoretical[6];

  for (int i=0; i <6 ; i++) {
    ClTheoretical[i]=WindowConv(band[i],Cl,50,90, Cl.stop()-300);
    //cout << ClTheoretical[i]<< endl;
  }

  //read experimental values
  ifstream dataFile(ControlPanel::cmbeasyDir("/resources/cbi/joint_3deep_std2_best_qb.dat").c_str());
  if (!dataFile) throw Bad_Error("AnalyzeThis::cbiDeepLikelihood(): dataFile not found");
  for (int i=0; i<6; i++) {
    dataFile >> dummy >> expData[i] >> dummy;
    //cout << expData[i];
  }
  dataFile.close();

  //read gaussianFischerMatrix
  ifstream matrixFile(ControlPanel::cmbeasyDir("/resources/cbi/GaussianFischerMatrixDeep.dat").c_str());
  if (!matrixFile) throw Bad_Error("AnalyzeThis::cbiDeepLikelihood(): matrixFile not found");
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++) {
      matrixFile >> FischerMatrix[i][j];
    }
  }

  dummy=0;
  
  // read offsets
  ifstream offsetFile(ControlPanel::cmbeasyDir("/resources/cbi/joint_3deep_std2_best_otherps.dat").c_str());
  if (!offsetFile) throw Bad_Error("AnalyzeThis::cbiDeepLikelihood(): offsetFile not found");
  for (int i=0; i<6; i++) {
    offset[i]=0;
    for ( int j=0; j<4; j++ ){ 
      offsetFile >> dummy;
      offset[i]+=dummy;
    }
  }
  
// compute chi squared
  for(int i=0; i<6; i++){
    for(int j=0; j < 6; j++){
      sum+=(log(expData[i]+offset[i])-log(ClTheoretical[i]+offset[i]))*FischerMatrix[i][j]*(log(expData[j]+offset[j])-log(ClTheoretical[j]+offset[j]));
    }
  }
  return sum;
}

/*!
  this will give the ACBAR chi^2
  written by Christian M. Mueller
*/
double AnalyzeThis::ACBARChi2(Spline &Cl){ 
  // name of window functions
  string band[14];
  double dummy;
  double offset[14];
  double expData[14];
  double variance[14];
  double sum=0;
  string DataFileName="acbar_data.dat";
 
  band[0]="acbar_window_1";  
  band[1]="acbar_window_2";
  band[2]="acbar_window_3";  
  band[3]="acbar_window_4";
  band[4]="acbar_window_5";  
  band[5]="acbar_window_6";
  band[6]="acbar_window_7";  
  band[7]="acbar_window_8";   
  band[8]="acbar_window_9";  
  band[9]="acbar_window_10";
  band[10]="acbar_window_11";  
  band[11]="acbar_window_12";
  band[12]="acbar_window_13";  
  band[13]="acbar_window_14";

  // find windowed theoretical values
  double ClTheoretical[14];

  for (int i=0; i <14 ; i++) {
    ClTheoretical[i]=WindowConvACBAR(band[i],Cl, Cl.stop()-300);
    //cout << ClTheoretical[i]<< endl;
  }

  //read experimental values
  ifstream dataFile(ControlPanel::cmbeasyDir("/resources/acbar/acbar_data.dat").c_str());
  for (int i=0; i<14; i++) {
    dataFile >> dummy >> expData[i] >> variance[i] >> offset[i];
  }
  dataFile.close();

// compute chi squared, a bit different from CBI since here the window functions are decorrelated, and as a result, the data is not correlated, so the Fischer Matrix is diagonal. Note that we transform from lognormal to normal here as well.
  for(int i=0; i<14; i++){
    if (acbarBandDeselect.find(i+1) == acbarBandDeselect.end()) {
      //cout << "acbar summing for band: " << i+1 << endl;
      sum+=(log(expData[i]+offset[i])-log(ClTheoretical[i]+offset[i]))*(log(expData[i]+offset[i])-log(ClTheoretical[i]+offset[i]))*(expData[i]+offset[i])/variance[i]*(expData[i]+offset[i])/variance[i];
    }
  }
  return sum;
}
/*!
  Compute chi2 of Lyman Alpha data using the routine of Pat McDonald.
  As the chi2 in Pat's routine is interpolated from a table, we use a gaussian
  approximation + penalty for values out of boundary. This should ensure
  convergence of monte carlo chains.  Please note that "gaussian approximation"
  is meant quite loosly. It is a very bad approximation for points closer to the
  boundary. It's sole purpose, as mentioned already is to guide the MCMC towards
  the region within the table boundaries. It will hence practically only be used during burn-in.
*/
double AnalyzeThis::lymanAlphaPatMcDonaldChi2(Cosmos& cosmos) {
  double z = 3;
  double tau3 = cosmos.z2tau(z), tau0=cosmos.tau_0();
  double k_pivoth = 0.009 * 100 * cosmos.tau2Hubble(tau3)/cosmos.tau2Hubble(tau0) / (1+z);  // k in Mpc{-1} / h
  double k_pivot = k_pivoth * cosmos.h(); // k in Mpc{-1}
  Anchor anchor; 
  double h = cosmos.h();
  Spline *Pz3 = cosmos.createPower(0,"hallo",cosmos.power_cdm(),&anchor,cosmos.z2tau(3));  
  Pz3->arm();

  Spline LnP(1000,"lnz", &anchor);  // ln P(k)
  Spline DlnP_Dlnk(&LnP,"dlnp");  // dln P / dlnk  = neff
  Spline D2lnP_Dlnk(&LnP,"dlnp2"); // d^2 lnP / dlnk = d neff / dnlk
  for (int i = 0; i < Pz3->size(); i++) {
    double k = Pz3->x(i) * h;  // power web stored as function of  k/h
    double lnk = log(k);
    double P = Pz3->y(i) / (h*h*h); // power stored with h^3 factored in, so we devide out
    double lnP = log(P);
    LnP.set(lnk,lnP);
  }
  LnP.arm();
  // LnP.dump("lnp",Spline::yes);
  LnP.derive(DlnP_Dlnk,Spline::yes);  // derive and arm()
  DlnP_Dlnk.derive(D2lnP_Dlnk,Spline::yes);  // derive and arm()
  // DlnP_Dlnk.dump("dlnp",Spline::yes);
 
  double Delta_L3 = cosmos.Pk2Delta2(Pz3, k_pivot);
  double neff = DlnP_Dlnk(log(k_pivot));
  double nrun = D2lnP_Dlnk(log(k_pivot));
  
  cout << "k_pivot [Mpc]: " << k_pivot << endl;
  cout << "Delta[k_pivot]: " << Delta_L3 << endl;
  cout << "neff[k_pivot]: " << neff << endl;
  cout << "running: " <<  nrun << endl;
  
  
  // first, we compute a gaussian approximation to the likelihood routine
  // we use this in cases where we operate outside the region covered by
  // the lyman alpha look up table of Pat McDonald. Please keep in mind
  // that this might constitute a severe theoretical bias, in case your model
  // deviates considerably from WMAP concordance cases.

  double x = -1.2*Delta_L3 + neff;   // rotate gaussian coordinates because of degeneracy
  double y = Delta_L3 + 0.8*neff;    // rotate gaussian coordinates because of degeneracy
   
  double chi2 = (x+2.84)*(x+2.84) / (2*0.05*0.05);
  chi2 += (y + 1.46)*(y + 1.46) / (2*0.11*0.11);
  //  chi2 += (nrun + 0.23)*(nrun + 0.23) /(2*0.05*0.05);  // penalty for models out of nrun range
  // add 504, i.e. minimum values of gaussian and full likelihood always off by penalty of about 410 
  chi2 += 504;                                   
  cout << "LYA chi2 gaussian approximation (including penalty): " << chi2 << endl;
  LyaFchi2Interpolator lyachi2;
  try {
    chi2 = lyachi2.chi2(Delta_L3,neff,nrun);
  } catch (  GlassInterpolatorOutOfBound x) {
    cout << "OUT OF BOUND" << endl;
  }
  cout << "final chi2: " << chi2  << endl;  
  return chi2;

}
