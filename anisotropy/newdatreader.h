#ifndef NEWDATREADER_H
#define NEWDATREADER_H

#include <vector>
#include <string>

/*! \class NewdatReader
 *  \brief class reads files in the .newdat format.
 *  The NewdatReader class handles data files in the .newdat format.
 *  The data files of e.g. ithe BOOMERANG03 NA pipeline
 *  (http://oberon.roma1.infn.it/boomerang/b2k/data/na.tgz)
 *  are shipped in this format.
 *  Extensive documentation of this format can be found in the README.newdat
 *  file that comes with that dataset.
 *  Usage of this class is straightforward:
 *
 *  \code
 *  NewdatReader reader( "filename.newdat" );
 *  // allBands are TT, EE, BB, EB, TE, TB, in this order
 *  std::vector<NewdatReader::Band> allBands = reader.allBands();
 *  double calibrationError = reader.calibrationError();
 *  // etc.
 *  \endcode
 *
 * For an example of how to use this class, see the BOOMERANG03
 * chi^2 handling in AnalyzeThis ( \ref AnalyzeThis::BOOMERANGInit() ).
 *
 * \warning The NewdatReader class does not handle the reading of the
 * corresponding window files.
 */

class NewdatReader
{
  public:
    //! NewdatReader constructor
    /*! Construct and initialize the NewdatReader instance.
     *  The data file is automatically read on construction.
     *
     *  \param fileName: the filename of the .newdat file,
     *  "e.g. /home/staff/someone/cmbeasy/resources/boomerang/B03_NA_21July05.newdat"
     *  \param skipLines: a vector of line numbers in the file header that will be ignored
     *                    this is useful for CBI2, which has some additional fields
     */
    NewdatReader( const char* fileName, const std::vector<unsigned int>& skipLines = std::vector<unsigned int>() );

    //! The root filename of the window files
    /*!  \return the root part of the filenames of the window functions,
     *   as read from the data file, so usually without any path information.
     */
    std::string windowFileNameRoot() const { return mWindowFileNameRoot; }

    //! likelihood type: 0==Gaussian, 1==all bands offset lognormal, 2==select bands offset lognormal
    int likelihoodType() const { return mLikeType; }


    int nrOfTTBands() const { return mNrOfTTBands; } //!< returns number of TT Bands
    int nrOfEEBands() const { return mNrOfEEBands; } //!< returns number of EE Bands
    int nrOfBBBands() const { return mNrOfBBBands; } //!< returns number of BB Bands
    int nrOfEBBands() const { return mNrOfEBBands; } //!< returns number of EB Bands
    int nrOfTEBands() const { return mNrOfTEBands; } //!< returns number of TE Bands
    int nrOfTBBands() const { return mNrOfTBBands; } //!< returns number of BB Bands
    /*! returns total number of Bands*/
    int totalNrOfBands() const { return mNrOfTTBands + mNrOfEEBands + mNrOfBBBands
                                     + mNrOfEBBands + mNrOfTEBands + mNrOfTBBands; }
    struct MinMaxBands
    {
      unsigned int min, max;
    };

    MinMaxBands bandSelectionTT() const { return mMinMaxTT; }
    MinMaxBands bandSelectionEE() const { return mMinMaxEE; }
    MinMaxBands bandSelectionBB() const { return mMinMaxBB; }
    MinMaxBands bandSelectionEB() const { return mMinMaxEB; }
    MinMaxBands bandSelectionTE() const { return mMinMaxTE; }
    MinMaxBands bandSelectionTB() const { return mMinMaxTB; }


    bool calibrationFlag() const { return mCalibrationFlag; } //!< return the calibration flag
    double recalibrationFactor() const { return mRecalibrationFactor;} //!< return the recalibration factor (in temperature)
    double calibrationError() const { return mCalibrationError; } //!< return the calibration error (in power)

    bool beamErrorFlag() const { return mBeamErrorFlag; } //!< return the beamErrorFlag
    double beamSize() const { return mBeamSize; } //!< return the beam size(in arcmins)
    double beamError() const { return mBeamError; } //!< return the beam error, converted to power (in readBand())

    /*! \struct Band
     * \brief Struct representing one Band
     * The Band struct has fields representing
     * the power in the band, errors, offsets and the like.
     */
    struct Band
    {
      double lMin, lMax, lEff; //!< minimal, maximal and effective ell contributing to this band
      double power;  //!< The power in this band
      double beamError; //<! the beam error(in power); for conversion see readBands()
      double plusError; //!< +err of the power
      double  minusError; //!< -err of the power
      bool offsetFlag; //!< whether this band has a lognormal offset */
      double offset; //!< lognormal offset for the lognormal approximation for band power likelihood */
    };

    std::vector<Band>& TTBands() { return mTTBands; } //!< returns a vector of the TT bands
    std::vector<Band>& EEBands() { return mEEBands; } //!< returns a vector of the EE bands
    std::vector<Band>& BBBands() { return mBBBands; } //!< returns a vector of the BB bands
    std::vector<Band>& EBBands() { return mEBBands; } //!< returns a vector of the EB bands
    std::vector<Band>& TEBands() { return mTEBands; } //!< returns a vector of the TE bands
    std::vector<Band>& TBBands() { return mTBBands; } //!< returns a vector of the TB bands

    std::vector<Band> allBands() const; //!< returns a vector of all bands

    typedef std::vector<std::vector<double> > Matrix; //!< Matrix is just a typedef for a vector< vector<double> >.
    Matrix& inverseFisherMatrix() { return mInverseFisherMatrix; } //!< returns the inverse Fisher matrix read from the file.

  protected:
    void read( const char* fileName ); //!< this method coordinates the acutal reading of the file; is called from the constructor
    void readBands( std::ifstream& file, std::vector<Band>& bands, const int nr ); //!< readBands does the actual reading of a read the power Band 

    //! set a vector of lines in the file header that will be ignored
    void setSkipHeaderLines( const std::vector<unsigned int>& l ) { mHeaderLinesToSkip = l; }

  private:
    std::string mWindowFileNameRoot;
    int mNrOfTTBands, mNrOfEEBands, mNrOfBBBands;
    int mNrOfEBBands, mNrOfTEBands, mNrOfTBBands;

    MinMaxBands mMinMaxTT, mMinMaxEE, mMinMaxBB,
                mMinMaxEB, mMinMaxTE, mMinMaxTB;

    int mLikeType;

    bool mCalibrationFlag, mBeamErrorFlag;
    double mRecalibrationFactor, mCalibrationError;
    double mBeamSize, mBeamError;

    std::vector<Band> mTTBands, mEEBands, mBBBands, mEBBands;
    std::vector<Band> mTEBands, mTBBands;

    Matrix mInverseFisherMatrix;

    std::vector<unsigned int> mHeaderLinesToSkip;
};
#endif // NEWDATREADER_H

