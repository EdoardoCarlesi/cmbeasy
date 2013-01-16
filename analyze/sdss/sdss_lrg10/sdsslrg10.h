#ifndef SDSSLRG10_H
#define SDSSLRG10_H

#include "controlpanel.h"
#include <string>
#include <vector>

class baseCosmos;
class Spline;

/*! \class SdssLrg10
 *  \brief likelihood from SDSS Luminos Red Galaxies
 *
 *  This class implements the likelihood routine for
 *  the SDSS lrg data from 
 *
 *  http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2966.2010.16276.x/full
 *
 *  Example usage:
 *     \code
 *      Spline* cdm=cosmos.createPower(0, "cdmPower", cosmos.power_cdm(),0,cosmos.z2tau(0));
 *      SdssLrg s(cosmos);
 *      double chisq = s.chi2margedOverb(*cdm);
 *     \endcode
 */

class SdssLrg10
{
  public:
    /*! \brief construct an SdssLrg object
     *
     * @param cosmos is the Cosmos object that produced the power spectrum
     * @param dataDir path to the sdss lrg data
     */
    SdssLrg10(const baseCosmos& c, std::string dataDir = ControlPanel::cmbeasyDir("/resources/sdss_lrg10/"));
    double chi2(Spline& pk, double bias = 1.); //!< returns the chi2 for sdss lrg
    //!< scans the bias range (firstBias, lastBias), returning the best chi2
    double chi2bestBias(Spline& pk, double firstBias = 1., double lastBias = 3., unsigned int steps = 20 );
    double chi2margedOverb(Spline& pk); //!< returns chi2 marginalized analytically over b^2 (flat prior)
    /*! returns chi2 marginalized analytically over b^2*Q (flat prior)
     * this code is directly ported from cosmomc, the comment there says:
     * > Marginalize analytically with flat prior on b^2 and b^2*Q
     * > as recommended by Max Tegmark for SDSS
     */
    double chi2margedOverb2Q(Spline& pk, double Ag = 1.4);
    void setMaxkBand( unsigned int n ) { mMaxBand = n; } //!< set max k band; default is 12; for the last 6 (of 20) "nonlinear modeling is definitely required" (see paper)

  protected:
    double aFactor();  //!< d_V(z)/d_V^{fiducial} from App. 4 of the paper
    std::vector<double> shiftedPower(Spline& pk, double bias);

  private:
    const baseCosmos& mCosmos;
    struct Band { double keff, klow, khigh, obs, sdev; };
    typedef std::vector<Band>::iterator BandIterator;
    std::vector<Band> mBands;
    typedef std::vector<std::vector<double> > Matrix;
    Matrix mMatrix;
    const std::string mDataDir;
    std::vector<double> mkvals;
    unsigned int mkbands;
    unsigned int mMaxBand;
};

#endif // SDSSLRG10_H
