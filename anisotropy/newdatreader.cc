#include "newdatreader.h"

#include "global.h"

#include <iterator>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>

NewdatReader::NewdatReader( const char* fileName, const std::vector<unsigned int>& skipLines )
    : mWindowFileNameRoot( std::string() ), mNrOfTTBands( 0 ), mNrOfEEBands( 0 ),
      mNrOfBBBands( 0 ), mNrOfEBBands( 0 ), mNrOfTEBands( 0 ),
      mNrOfTBBands( 0 ), mCalibrationFlag( false ), mBeamErrorFlag( false ),
      mRecalibrationFactor( 0 ), mCalibrationError( 0 ), mBeamSize( 0 ),
      mBeamError( 0 ), mTTBands( std::vector<Band>() ), mEEBands( std::vector<Band>() ),
      mBBBands( std::vector<Band>() ), mEBBands( std::vector<Band>() ),
      mTEBands( std::vector<Band>() ), mTBBands( std::vector<Band>() ),
      mInverseFisherMatrix( Matrix() )
{
  setSkipHeaderLines(skipLines);
  read( fileName );
}

std::vector<NewdatReader::Band> NewdatReader::allBands() const
{
  std::vector<Band> retVal;
  std::copy( mTTBands.begin(), mTTBands.end(), std::back_inserter( retVal ) );
  std::copy( mEEBands.begin(), mEEBands.end(), std::back_inserter( retVal ) );
  std::copy( mBBBands.begin(), mBBBands.end(), std::back_inserter( retVal ) );
  std::copy( mEBBands.begin(), mEBBands.end(), std::back_inserter( retVal ) );
  std::copy( mTEBands.begin(), mTEBands.end(), std::back_inserter( retVal ) );
  std::copy( mTBBands.begin(), mTBBands.end(), std::back_inserter( retVal ) );

 return retVal;
}

/*
static NewdatReader::Band dummyBand()
{
  NewdatReader::Band dummy;
  dummy.power = 0;
  dummy.offset = 0;
  dummy.offsetFlag = false;
  return dummy;
}
*/

void NewdatReader::readBands( std::ifstream& file, std::vector<Band>& bands, const int nr )
{
  //cout << "-------enter" << endl;
  int bandNr = 1;
  while ( bandNr < nr )
  {
    NewdatReader::Band current;
    file >> bandNr >> current.power >> current.plusError >> current.minusError >> current.offset >> current.lMin
         >> current.lMax;
    double l_effm1 = current.lMin;
    double l_effp1 = current.lMax;
    if (likelihoodType()==2) {
      file >> current.offsetFlag;
    } else {
      current.offsetFlag = ((likelihoodType()==1)?true:false);
    }

    const double l_eff = l_effm1 + ( l_effp1 - l_effm1 ) / 2.;
    current.lEff = l_eff;

    if ( mBeamErrorFlag )
    {
      //beam error according to eq. (4) of astro-ph/0507494; only first order (delta^2/delta is ~.02)
      //fractional error delta is mBeamError/mBeamSize
      // mBeamError and mBeamSize are in arcminutes, therefore conversion factor is 2*pi/(360*60)
      // FWHM factor is 1/sqrt(8*ln(2)), that squared times conversion factor squared gives 1.526e-8
      current.beamError = fabs( exp( - l_eff * ( l_eff + 1 ) * 1.526e-8 * 2. * mBeamError * mBeamSize ) - 1. );
    }
    else
    {
      current.beamError = 0;
    }
    //cout << "read band number " << bandNr << " of " << nr << endl;
    //cout << bandNr << " " << current.power << current.offset << endl;
    bands.push_back( current );
  }
  //cout << "------leave" << endl;
}

//! helper function for reading newdat files
static void seekToNextLine(unsigned int& currentLine,
                           const vector<unsigned int>& linesToSkip,
                           std::ifstream& file)
{
    ++currentLine;
    while (std::find(linesToSkip.begin(), linesToSkip.end(), currentLine)
                                            != linesToSkip.end())
    {
        file.ignore(4096, '\n');
        ++currentLine;
    }
}

void NewdatReader::read( const char* fileName)
{
  std::ifstream dataFile( fileName );
  if (!dataFile)
    throw Bad_Error("NewdatReader: dataFile " + std::string(fileName) + " not found");

  unsigned int currentHeaderLine = 0;
  seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
  dataFile >> mWindowFileNameRoot;

  seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
  dataFile >> mNrOfTTBands >> mNrOfEEBands >> mNrOfBBBands >> mNrOfEBBands >> mNrOfTEBands >> mNrOfTBBands;


  dataFile.ignore(4096, '\n');
  seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
  std::stringstream line;
  std::string lineStr;
  std::getline(dataFile, lineStr);
  line << lineStr;

  if (line.str() == "BAND_SELECTION") {
     seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
     dataFile >> mMinMaxTT.min >> mMinMaxTT.max;                        // TT
     seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
     dataFile >> mMinMaxEE.min >> mMinMaxEE.max;                        // EE
     seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
     dataFile >> mMinMaxBB.min >> mMinMaxBB.max;                        // BB
     seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
     dataFile >> mMinMaxEB.min >> mMinMaxEB.max;                        // EB
     seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
     dataFile >> mMinMaxTE.min >> mMinMaxTE.max;                        // TE
     seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
     dataFile >> mMinMaxTB.min >> mMinMaxTB.max;                        // TB

     seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
     dataFile >> mCalibrationFlag >> mRecalibrationFactor >> mCalibrationError;
   } else {
     line     >> mCalibrationFlag >> mRecalibrationFactor >> mCalibrationError;
   }

  //cout << " read calib: " << mCalibrationFlag << " " << mRecalibrationFactor << " " << mCalibrationError << endl;

  seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);
  dataFile >> mBeamErrorFlag >> mBeamSize >> mBeamError;

  seekToNextLine( currentHeaderLine, mHeaderLinesToSkip, dataFile);

  dataFile >> mLikeType;
  //cout << mLikeType << " is the likelihoodtype" << endl;
  dataFile.ignore( 4096, '\n' );


  if ( mNrOfTTBands > 0 )
  {
    std::string polarizationType;
    dataFile >> polarizationType;
    //cout << polarizationType << " is the polarizationType" << endl;
    if ( polarizationType != "TT" )
      throw Bad_Error( "NewdatReader: Error parsing data File: TT powers expected. " );
    readBands( dataFile, mTTBands, mNrOfTTBands );
    // skip band to band correlation
    int skippedLines = 0; while ( skippedLines++ <= mNrOfTTBands ) dataFile.ignore( 8192, '\n' );
  }

  if ( mNrOfEEBands > 0 )
  {
    std::string polarizationType;
    dataFile >> polarizationType;
    //cout << polarizationType << " is the polarizationType" << endl;
    if ( polarizationType != "EE" )
      throw Bad_Error( "NewdatReader: Error parsing data File: EE powers expected. " );

    readBands( dataFile, mEEBands, mNrOfEEBands );
    // skip band to band correlation
    int skippedLines = 0; while ( skippedLines++ <= mNrOfEEBands ) dataFile.ignore( 8192, '\n' );
  }

  if ( mNrOfBBBands > 0 )
  {
    std::string polarizationType;
    dataFile >> polarizationType;
    if ( polarizationType != "BB" )
      throw Bad_Error( "NewdatReader: Error parsing data File: BB powers expected. " );

    readBands( dataFile, mBBBands, mNrOfBBBands );
    // skip band to band correlation
    int skippedLines = 0; while ( skippedLines++ <= mNrOfBBBands ) dataFile.ignore( 8192, '\n' );
  }

  if ( mNrOfEBBands > 0 )
  {
    std::string polarizationType;
    dataFile >> polarizationType;
    if ( polarizationType != "EB" )
      throw Bad_Error( "NewdatReader: Error parsing data File: EB powers expected. " );

    readBands( dataFile, mEBBands, mNrOfEBBands );
    // skip band to band correlation
    int skippedLines = 0; while ( skippedLines++ <= mNrOfEBBands ) dataFile.ignore( 8192, '\n' );
  }

  if ( mNrOfTEBands > 0 )
  {
    std::string polarizationType;
    dataFile >> polarizationType;
    if ( polarizationType != "TE" )
      throw Bad_Error( "NewdatReader: Error parsing data File: TE powers expected. " );

    readBands( dataFile, mTEBands, mNrOfTEBands );
    // skip band to band correlation
    int skippedLines = 0; while ( skippedLines++ <= mNrOfTEBands ) dataFile.ignore( 8192, '\n' );
  }

  if ( mNrOfTBBands > 0 )
  {
    std::string polarizationType;
    dataFile >> polarizationType;
    if ( polarizationType != "TB" )
      throw Bad_Error( "NewdatReader: Error parsing data File: TB powers expected. " );

    readBands( dataFile, mTBBands, mNrOfTBBands );
    // skip band to band correlation
    int skippedLines = 0; while ( skippedLines++ <= mNrOfTBBands ) dataFile.ignore( 8192, '\n' );
  }

  // read inverse Fisher matrix
  const unsigned int totalBandNr = totalNrOfBands();

  for ( unsigned int i=0; i < totalBandNr; ++i )
  {
    vector<double> row;

    for ( unsigned int j=0; j < totalBandNr; ++j )
    {
      double entry;
      dataFile >> entry;
      row.push_back( entry );
    }

    mInverseFisherMatrix.push_back( row );
  }

  dataFile.close();
}

