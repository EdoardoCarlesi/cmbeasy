#ifndef NEWDATCLCHI2_H
#define NEWDATCLCHI2_H

#include "newdatreader.h"

#include <map>
#include <vector>

class CL;
class Spline;

/*! \class NewdatClChi2
 *
 *   Class to compute the chi2 of a given Cl-spectrum for data in the .newdat file format.
 *
 *   Example usage:
 *   \code
 *   NewdatClChi2 data(cl, "datadir/", "datafile.newdat");
 *   // deselect some bands
 *   data.bandDeselect[0] = true;
 *   data.bandDeselect[14] = true;
 *   double chisq = data.chi2WithCalibration();
 *   \endcode
 */
class NewdatClChi2
{
  public:
    NewdatClChi2(CL&, const std::string& dataDir, const std::string& dataFile);
    double chi2( double beamErrorAdjust=0, double calibrationFactor=1 );
    double chi2WithCalibration(); //!< Cl's in muK^2 

    //! All integers appearing here will be considered de-selected from the likelihood.
    std::map<int,bool> bandDeselect;

  protected:
    double windowConv(int windowNumber,Spline &Cl, double cut_l);
    void init();

  private:
    const CL& mCl;
    std::vector<double> mThband;
    bool initialized;
    std::map<int,bool> mSavedBandDeselect;
    std::vector<NewdatReader::Band> mExpData;
    NewdatReader::Matrix mMatrix;
    std::string mDataFileName, mDataDir, mWindowFilesRootName;
    int mBandsCount;
    double mCalibrationError;
    double mHaveBeamError;
};

#endif // NEWDATCLCHI2_H
