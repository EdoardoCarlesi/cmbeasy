#include "cbi2.h"

#include "spline.h"
#include "cl.h"
#include "controlpanel.h"

/*
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
*/

#include <numeric>
#include <algorithm>

CBI2::CBI2( CL &cl ): mCl(cl), initialized(false) {}

void CBI2::init()
{
  mSavedBandDeselect = bandDeselect;

  Spline ClTT(*mCl.ts[0],"CBI2Init_ts");
  // divide out the conventional l(l+1)/2*pi
  for (int i = ClTT.first(); i <= ClTT.last(); ++i ) {
      double l = ClTT.x(i);
      ClTT.setY(i, ClTT.y(i) * 2. * M_PI / (l * (l+1)));
  }

  ClTT.arm();

  const string dataFileName=ControlPanel::cmbeasyDir("/resources/CBI2.0_newdat/CBI2.0_2_to_14.newdat");
  vector<unsigned int> headerLinesToSkip;
  unsigned int lines[9] = {1}; // T_CMB ('old' format)

  std::copy(lines, lines+9, std::back_inserter(headerLinesToSkip));

  NewdatReader data(dataFileName.c_str(), headerLinesToSkip);
  mExpData = data.allBands();

  const int bands=data.totalNrOfBands();

  std::fill_n( mThband, bands, 0. );

  map<int,bool> &deselect = bandDeselect;
  for(int i=0; i<bands; ++i){
    if (deselect.find(i) == deselect.end()){
      mThband[i] = windowConv(i+1,ClTT, ClTT.stop()-300);
    }
  }
  ClTT.disarm();

  // denormalize cov matrix
  for(int i=0; i<bands; i++)
      mExpData[i].power = log(mExpData[i].power+mExpData[i].offset);


/*

// The following block reads the inverse Fisher matrix directly from the
// data file. For speed reasons, and because we don't want to link to the GSL
// here, we read it from a file that was constructed using the code below
// and is shipped with CMBEASY.

  NewdatReader::Matrix invCBI2Matrix = data.inverseFisherMatrix();

  for(int i=0; i<bands; i++)
      for(int j=0; j < bands; j++)
      {
        if ( i < data.nrOfTTBands() && j < data.nrOfTTBands() )
           // FISHER_T_CMB at top of file apparently means that this is in units of T_CMB, not \muK^2
          invCBI2Matrix[i][j] *= pow(2.725e6,4);
          invCBI2Matrix[i][j] /= (mExpData[i].power+mExpData[i].offset)*(mExpData[j].power+mExpData[j].offset);
      }

  gsl_permutation* perm = gsl_permutation_alloc( bands );
  gsl_matrix* inverse = gsl_matrix_alloc( bands, bands );
  gsl_matrix* fisher = gsl_matrix_alloc( bands, bands );

  for ( int i = 0; i < bands; ++i )
    for ( int j = 0; j < bands; ++j )
      gsl_matrix_set( inverse, i, j, invCBI2Matrix[i][j] );


  int status;
  gsl_linalg_LU_decomp( inverse, perm, &status );
  gsl_linalg_LU_invert( inverse, perm, fisher );

  mMatrix=NewdatReader::Matrix( bands );
  for ( int i = 0; i < bands; ++i)
    mMatrix[ i ].resize( bands );


  for ( int i = 0; i < bands; ++i )
    for ( int j = 0; j < bands; ++j )
      mMatrix[i][j] = gsl_matrix_get( fisher, i, j );


  string fisherFileName = ControlPanel::cmbeasyDir("/resources/CBI2.0_newdat/Matrix.dat");
  ofstream out( fisherFileName.c_str() );
  for ( int i = 0; i < bands; ++i )
  {
    for ( int j = 0; j < bands; ++j )
      //!!!out << std::scientific << CBI2Matrix[i][j] << " ";
      //out << std::scientific << ((mExpData[i].power+mExpData[i].offset)*mMatrix[i][j]*(mExpData[j].power+mExpData[j].offset))<< " ";
      out << std::scientific << mMatrix[i][j] << " ";
    out << endl;
  }
  out.close();
*/
//=================================================================

  // read the Matrix from file
  mMatrix = NewdatReader::Matrix( bands );
  for ( int i = 0; i < bands; ++i)
    mMatrix[ i ].resize( bands );

  string fisherFile = ControlPanel::cmbeasyDir("/resources/CBI2.0_newdat/Matrix.dat");
  ifstream in( fisherFile.c_str() );
  for ( int i = 0; i < bands; ++i )
    for ( int j = 0; j < bands; ++j )
    {
      in >> std::scientific >> mMatrix[i][j];
    }
  in.close();
}

double CBI2::windowConv(int windowNumber, Spline &Cl, double cut_l){
  //open WindowFile
  string windowName[15]={"cbi7_wins1","cbi7_wins2","cbi7_wins3","cbi7_wins4","cbi7_wins5","cbi7_wins6", "cbi7_wins7","cbi7_wins8","cbi7_wins9","cbi7_wins10","cbi7_wins11","cbi7_wins12",
                         "cbi7_wins13","cbi7_wins14","cbi7_wins15"};

  string windowFileName = ControlPanel::cmbeasyDir("/resources/CBI2.0_newdat/windows/")+windowName[windowNumber-1];
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


  // now, multiply theoretical data with the window function 

  double sum=0;
  //double windowNorm=0;
  for (int j=0; j<bins; j++) {
    if (begin+j < cut_l) {
      sum += window[j]*Cl(begin+j)*(begin+j)*(begin+j+0.5)/(2*M_PI);
    }
  }
  return sum;
}

/*!
 *  chi2 returns the chi^2 w.r.t. CBI2 03, optionally taking
 *  a beam error and/or calibration error into account. If you want the marginalized
 *  version ( over beam error and calibration uncertainty ) of chi^2, use 
 *  chi2WithCalibration() instead.
 *
 *  \param beamErrorAdjust add beamErrorAdjust * beamError(i) to
 *                         theoretical band power(i) before calculating chi^2
 *  \param calibrationFactor multiply theoretical bandpowers by this factor before comparing
 */
double CBI2::chi2( double beamErrorAdjust, double calibrationFactor )
{
  if ((!initialized) || (mSavedBandDeselect != bandDeselect)) init();
  initialized = true;
  const int bands = 15;
  double sum=0;
  double ClTheoretical[bands];
  double ClExperimental[bands];

  map<int, bool>& deselect = bandDeselect;

  double thband[bands];
  std::fill_n( ClTheoretical, bands, 0. );
  std::fill_n( ClExperimental, bands, 0. );
  std::copy( mThband, mThband+bands, thband );

  for(int i=0; i<bands; ++i){
    if (deselect.find(i) == deselect.end() ){
      if ( beamErrorAdjust != 0 )
        thband[i] += thband[i]*beamErrorAdjust*mExpData[i].beamError;
      thband[i] *= calibrationFactor;
    }
  }
  for(int i=0; i<bands; i++){
    if (deselect.find(i) == deselect.end()){
      if (mExpData[i].offsetFlag){
        ClTheoretical[i]=log(thband[i]+mExpData[i].offset);
        ClExperimental[i]=mExpData[i].power;
      }
      else {
        ClTheoretical[i]=thband[i];
        ClExperimental[i]=mExpData[i].power;
      }
    }
  }
  //copy( ClTheoretical, ClTheoretical + bands, ostream_iterator<double>(cout, "\n") );

  for(int i=0; i<bands; i++){
    if ( deselect.find(i) == deselect.end()) // if the bands have not been de-selected
      for(int j=0; j < bands; j++){
        if ( deselect.find(j) == deselect.end()) {
          double contrib =(-ClExperimental[i]+ClTheoretical[i])*mMatrix[i][j]*(-ClExperimental[j]+ClTheoretical[j]);
          sum += contrib;
        }
      }
  }
  return sum;
}


/*!
 *  returns the chi^2 w.r.t. CBI2 03, marginalized over the calibration error.
 */
double CBI2::chi2WithCalibration()
{
  //marginalize over calibration error, assuming gaussian error
  vector<double> likeCalib, weights;
  const double steps = 5;
  //allow for +/- 3 sigma
  for ( double beam = -steps; beam <= steps; ++beam )
    weights.push_back( exp( -pow( beam*3./steps, 2. ) /2. ) );
  double norm = std::accumulate( weights.begin(), weights.end(), 0. );

  for (double f = -steps; f <= steps; ++f) {
    //0.026 is CBI2 calibration uncertainty
    double cbi2 = chi2(0.,1+0.026*f*3./steps);
    likeCalib.push_back( cbi2 );
  }
    double calibMin = *std::min_element( likeCalib.begin(), likeCalib.end() ) ;
    for ( unsigned int i = 0; i < likeCalib.size(); ++i )
      likeCalib[i] = exp( -( likeCalib[i]-calibMin )/2. )*weights[i];
  return ( -2. * log( accumulate( likeCalib.begin(), likeCalib.end(), 0. )/norm ) + calibMin );
}
