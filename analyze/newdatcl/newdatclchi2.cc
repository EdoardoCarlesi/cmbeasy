#include "newdatclchi2.h"

#include "spline.h"
#include "cl.h"
#include "controlpanel.h"


#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"


#include <numeric>
#include <algorithm>
#include <sstream>

NewdatClChi2::NewdatClChi2(CL& cl, const std::string& dataDir, const std::string& dataFile)
                    : mCl(cl), initialized(false)
{
  mDataDir = dataDir ;
  if (mDataDir.size()>1 && mDataDir[mDataDir.size()-1] != '/')
    mDataDir.push_back('/');

  mDataFileName = mDataDir + dataFile;

}

void NewdatClChi2::init()
{

  Spline ClTT(*mCl.ts[0],"NewdatClChi2Init_ts");
  // divide out the conventional l(l+1)/2*pi
  for (int i = ClTT.first(); i <= ClTT.last(); ++i ) {
      double l = ClTT.x(i);
      ClTT.setY(i, ClTT.y(i) * 2. * M_PI / (l * (l+1)));
  }

  ClTT.arm();

  vector<unsigned int> headerLinesToSkip;

  //unsigned int lines[9] = {5,6,7,8,9,10};
  //std::copy(lines, lines+9, std::back_inserter(headerLinesToSkip));

  NewdatReader data(mDataFileName.c_str(), headerLinesToSkip);
  mExpData = data.allBands();

  if (data.nrOfTTBands() != data.totalNrOfBands())
    throw Bad_Error("NewdatClChi2::init(): fix deselection of non-TT bands");

  for (unsigned int i = 1; i <= data.nrOfTTBands(); ++i) {
    if (i<data.bandSelectionTT().min || i>data.bandSelectionTT().max) {
      bandDeselect[i] = true;
    }
  }

  mSavedBandDeselect = bandDeselect;

  mWindowFilesRootName = data.windowFileNameRoot();
  mCalibrationError = data.calibrationError();
  mHaveBeamError = data.beamErrorFlag();

  mBandsCount = data.totalNrOfBands();
  int bands = mBandsCount;

  mThband.resize(bands);
  std::fill( mThband.begin(), mThband.end(), 0. );

  map<int,bool> &deselect = bandDeselect;
  for(int i=0; i<bands; ++i){
    if (deselect.find(i) == deselect.end()){
      mThband[i] = windowConv(i+1,ClTT, ClTT.stop()-300);
    }
  }
  ClTT.disarm();




// The following block reads the inverse Fisher matrix directly from the
// data file. For speed reasons, and because we don't want to link to the GSL
// here, we read it from a file that was constructed using the code below
// and is shipped with CMBEASY.
/*
  cout << "WARNING: regenerating ACBAR 08 Fisher matrix" << endl;
  NewdatReader::Matrix invNewdatClChi2Matrix = data.inverseFisherMatrix();

  for(int i=0; i<bands; i++)
      for(int j=0; j < bands; j++)
      {
        if ( i < data.nrOfTTBands() && j < data.nrOfTTBands() )
//X            // only if FISHER_T_CMB at top of newdat file (means that this is in units of T_CMB, not \muK^2)
//X           invNewdatClChi2Matrix[i][j] *= pow(2.725e6,4);
          invNewdatClChi2Matrix[i][j] /= (mExpData[i].power+mExpData[i].offset)*(mExpData[j].power+mExpData[j].offset);
      }

  gsl_permutation* perm = gsl_permutation_alloc( bands );
  gsl_matrix* inverse = gsl_matrix_alloc( bands, bands );
  gsl_matrix* fisher = gsl_matrix_alloc( bands, bands );

  for ( int i = 0; i < bands; ++i )
    for ( int j = 0; j < bands; ++j )
      gsl_matrix_set( inverse, i, j, invNewdatClChi2Matrix[i][j] );


  int status;
  gsl_linalg_LU_decomp( inverse, perm, &status );
  gsl_linalg_LU_invert( inverse, perm, fisher );

  mMatrix=NewdatReader::Matrix( bands );
  for ( int i = 0; i < bands; ++i)
    mMatrix[ i ].resize( bands );


  for ( int i = 0; i < bands; ++i )
    for ( int j = 0; j < bands; ++j )
      mMatrix[i][j] = gsl_matrix_get( fisher, i, j );


  string fisherFileName = mDataDir+"Matrix.dat";
  ofstream out( fisherFileName.c_str() );
  for ( int i = 0; i < bands; ++i )
  {
    for ( int j = 0; j < bands; ++j )
      //!!!out << std::scientific << NewdatClChi2Matrix[i][j] << " ";
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

  string fisherFile = mDataDir + "Matrix.dat";
  ifstream in( fisherFile.c_str() );
  if (!in) {
    throw Bad_Error("NewdatClChi2::init() - Could not read file:\n\t"+fisherFile);
  }
  for ( int i = 0; i < bands; ++i )
    for ( int j = 0; j < bands; ++j )
    {
      in >> std::scientific >> mMatrix[i][j];
      //cout << i << " " << j << ": " << mMatrix[i][j] << endl;
    }
  in.close();
  // convert to lognormal
  for(int i=0; i<bands; i++)
    if (mExpData[i].offsetFlag)
      mExpData[i].power = log(mExpData[i].power+mExpData[i].offset);

}

double NewdatClChi2::windowConv(int windowNumber, Spline &Cl, double cut_l)
{
  stringstream wFileName;
  wFileName << mDataDir << "windows/" << mWindowFilesRootName << windowNumber;
  string windowFileName = wFileName.str();
  ifstream windowFile(windowFileName.c_str());
  if (!windowFile)
    cout << "Error in NewdatClChi2WindowConv: Could not open " << windowFileName << endl;

  map<unsigned int, double> window;
  typedef map<unsigned int, double>::const_iterator WindowIterator;
  unsigned int l;
  double val;
  while ( windowFile >> l >> val) {
    window[l] = val;
  }


  // now, multiply theoretical data with the window function 

  double sum=0;
  //double windowNorm=0;
  WindowIterator it, windowsEnd;
  windowsEnd = window.end();

  for ( it = window.begin(); it != windowsEnd; ++it) {
    int l = it->first;
    double windowVal = it->second;
    if (l > cut_l) {
      // throw Bad_Error("NewdatClChi2::windowConv() - l_cut too small to cover entire window function");
      // assume no contribution
      cout << "warning: NewdatClChi2::windowConv() - " << mDataFileName << endl;
      cout << "window " << windowNumber << ", l>cut_l at l=" << l << endl;
      break;
    }
    sum += windowVal*Cl(l)*(l)*(l+0.5)/(2.*M_PI);
  }
  return sum;
}

/*!
 *  chi2 returns the chi^2 w.r.t. to the cl dataset, optionally taking
 *  a beam error and/or calibration error into account. If you want the marginalized
 *  version ( over beam error and calibration uncertainty ) of chi^2, use 
 *  chi2WithCalibration() instead.
 *
 *  \param beamErrorAdjust add beamErrorAdjust * beamError(i) to
 *                         theoretical band power(i) before calculating chi^2
 *  \param calibrationFactor multiply theoretical bandpowers by this factor before comparing
 */
double NewdatClChi2::chi2( double beamErrorAdjust, double calibrationFactor )
{
  if ((!initialized) || (mSavedBandDeselect != bandDeselect))
    init();
  initialized = true;
  const int bands = mBandsCount;
  double sum=0;
  std::vector<double> ClTheoretical(bands, 0.);
  std::vector<double> ClExperimental(bands, 0.);

  map<int, bool>& deselect = bandDeselect;

  std::vector<double> thband(mThband);

  //copy( mThband.begin(), mThband.end(), ostream_iterator<double>(cout, "\n") );
  //copy( thband, thband+bands, ostream_iterator<double>(cout, "\n") );

  for(int i=0; i<bands; ++i){
    if (deselect.find(i) == deselect.end() ){
      //cout << thband[i] << "  " << calibrationFactor << endl;
      if ( beamErrorAdjust != 0 )
        thband[i] += thband[i]*beamErrorAdjust*mExpData[i].beamError;
      thband[i] *= calibrationFactor;
    }
  }
  for(int i=0; i<bands; i++){
    if (deselect.find(i) == deselect.end()){
      if (mExpData[i].offsetFlag){
        ClTheoretical[i]=log(thband[i]+mExpData[i].offset);
        //cout << i << " " << ClTheoretical[i] << " - " << thband[i] << " " << mExpData[i].offset << endl;
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
 *  returns the chi^2 w.r.t. data from a newdat file, marginalized over the beam+calibration error.
 */
double NewdatClChi2::chi2WithCalibration()
{
  if ((!initialized) || (mSavedBandDeselect != bandDeselect))
    init();
  initialized = true;
  using namespace std;
  //marginalize over beam+calibration error, assuming gaussian error
  vector<double> likeCalib, likeBeam, weights;
  const double steps = 5;
  //allow for +/- 3 sigma
  for ( double beam = -steps; beam <= steps; ++beam )
    weights.push_back( exp( -pow( beam*3./steps, 2. ) /2. ) );
  double norm = accumulate( weights.begin(), weights.end(), 0. );

  for ( double beam = -steps; beam <= steps; ++beam ) {
    likeCalib.clear();
    for (double f = -steps; f <= steps; ++f) {
      double val = chi2(beam*3./steps, 1.+mCalibrationError*f*3./steps);
      likeCalib.push_back( val );
    }

    double calibMin = *std::min_element( likeCalib.begin(), likeCalib.end() ) ;
    for ( unsigned int i = 0; i < likeCalib.size(); ++i )
      likeCalib[i] = exp( -( likeCalib[i]-calibMin )/2. )*weights[i];

    likeBeam.push_back( -2. * log( accumulate( likeCalib.begin(), likeCalib.end(), 0. )/norm ) + calibMin );
  }
  double beamMin = *std::min_element( likeBeam.begin(), likeBeam.end() ) ;
  for ( unsigned int i = 0; i < likeBeam.size(); ++i )
    likeBeam[i] = exp( -( likeBeam[i]-beamMin )/2. )*weights[i];

  return ( -2. * log( accumulate( likeBeam.begin(), likeBeam.end(), 0. )/norm ) + beamMin );
}
