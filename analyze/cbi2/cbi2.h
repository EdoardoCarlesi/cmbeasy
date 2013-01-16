#ifndef CBI2_H
#define CBI2_H

#include "newdatreader.h"

#include <map>
#include <vector>

class CL;
class Spline;

/*! \class CBI2
 *
 *   Class to compute the CBI2 likelihood of a given Cl-spectrum.
 *
 *   Example usage:
 *   \code
 *   CBI2 cbi2(cl);
 *   // deselect the very first and last bands
 *   cbi2.bandDeselect[0] = true;
 *   cbi2.bandDeselect[14] = true;
 *   double chisq = cbi2.chi2WithCalibration();
 *   \endcode
 */
class CBI2
{
  public:
    CBI2(CL&);
    double chi2( double beamErrorAdjust=0, double calibrationFactor=1 );
    double chi2WithCalibration(); //!< Cl's in muK^2 

    //! All integers appearing here will be considered de-selected from the CBI2000+01 likelihood. Counting goes from 0 to 14.
    std::map<int,bool> bandDeselect;

  protected:
    double windowConv(int windowNumber,Spline &Cl, double cut_l);
    void init();

  private:
    const CL& mCl;
    double mThband[15];
    bool initialized;
    std::map<int,bool> mSavedBandDeselect;
    std::vector<NewdatReader::Band> mExpData;
    NewdatReader::Matrix mMatrix;
};

#endif // CBI2_H
