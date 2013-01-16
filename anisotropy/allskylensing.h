#ifndef ALLSKYLENSING
#define ALLSKYLENSING

#include "spline.h"
#include "cmbcalc.h"

#include <string>
#include <vector>

class Cosmos;
class ControlPanel;
class CL;

/*! \class AllSkyLensing
 *  \brief All-Sky-Lensing from astro-ph/0502425
 *
 *  This class implements the lensing of the CMB power spectra using the all-sky
 *  correlation function method of astro-ph/0502425.
 *  The old flat-sky method from cmbfast is still available in the \ref FlatSkyLensing class.
 *
 *  The default accuracy settings should be sufficient for most purposes; if you need
 *  higher accuracy, change the #defines in allskylensing.cpp.
 *
 *  You can either use the lensedCls() method to get the lensed Cls or call
 *  dumpLensedCls() with a file name to have the lensed spectra written to disc.
 *
 */

class AllSkyLensing
{
  public:
    /*! \brief construct an AllSkyLensing object for lensing cl-spectra
     *
     *  @param cosmos is a Cosmos object that provides the lensing potential power spectrum
     *  @param controlPanel a ControlPanel holding the ( lensing ) settings
     *  @param cl the unlensed spectra
     *  @param cmbCalc used to write out the lensed spectra
     */
    AllSkyLensing( Cosmos& cosmos, const ControlPanel& controlPanel, const CL& cl, const CmbCalc& cmbCalc );

    /*! \brief calculate the lensed spectra
     *
     *  @return the lensed version of the CL passed to the constructor
     */
    CL* lensedCls();

    /*! \brief calculate and write the lensed spectra to disc
     *
     * This is a convenience function. It calls lensedCls()
     * and writes the result to the file fileName.
     *
     * @param fileName the file to write to
     */
    void dumpLensedCls( const string fileName );


  private:
    struct LensingDifference;
    struct DlSplines;
    struct Ximn;

    /*! \brief calculate the lensing contribution from one theta
     *
     *  Calculates the lensing contribution for the different spectra and adds them
     *  to the corresponding members of LensingDifference diffs. They are added to
     *  the original ( unlensed ) spectra in lensedCls( ).
     *  lensedCls() loops over the angle range in order to avoid redundant computations
     *  of DlSplines.
     *
     *  @param n is the spectral index of the spectra the lensing correction is being
     *  computed for
     *  @param theta the angle to compute the lensing correction for
     *  @param diffs the LensingDifference struct to hold the calculated corrections
     *
     *  @see lensedCls()
     *  @see LensingDifference
     */
    void oneAngle( int n, double theta, LensingDifference * diffs );

    /*! \brief return sigma^2(beta)
     *
     *  @param n the spectral index
     *  @param beta the angle to compute sigma^2 for
     *  @param s precalculated DlSplines for beta
     *
     *  @return sigma^2 = C_gl( beta )-C_gl( 0 )
     */
    double sigma2( int n, double beta, const DlSplines& s );

    /*! \brief return C_gl,2(beta)
     *
     *  @param n the spectral index
     *  @param beta the angle to compute C_gl,2 for
     *  @param s precalculated DlSplines for beta
     *
     *  @return C_gl,2( beta )
     */
    double C_gl2( int n, double beta, const DlSplines& s );

    /*! \brief compute lensing contribution for ts-spectrum
     *
     *  computes xiLensed( beta ) - xiUnlensed( beta ) for spectral index n
     *
     *  @param n the spectral index
     *  @param beta angle to compute the difference for
     *  @param cgl2 C_gl,2( beta )
     *  @param s precomputed DlSplines for angle beta
     *  @param x precomputed X_imns for angle beta
     *
     *  @return lensing contribution for angle beta
     */
    double xiLensed( int n, double beta, double cgl2, const DlSplines& s, const Ximn& x );

    /*! \brief compute lensing contribution for e + b
     *
     *  essentially the same as above, but for e+b
     *  @see xiLensed
     *  @return lensing contribution for e+b
     */
    double xiLensed_plus( int n, double beta, double cgl2, const DlSplines& s, const Ximn& x );

    /*! \brief compute lensing contribution for e - b
     *
     *  essentially the same as above, but for e+b
     *  @see xiLensed
     *  @return lensing contribution for e+b
     */
    double xiLensed_minus( int n, double beta, double cgl2, const DlSplines& s, const Ximn& x );

    /*! \brief compute lensing contribution for cs
     *
     *  essentially the same as above, but for cross-correlation T-E
     *  @see xiLensed
     *  @return lensing contribution for cross-correlation
     */
    double xiLensed_X( int n, double beta, double cgl2, const DlSplines& s, const Ximn& x );

    //! conventional l-factor
    inline double norml2( double l ) { return 16 * M_PI * M_PI * l * (l + 1); }

  private:
    Cosmos& mCosmos;
    const CmbCalc& mCmbCalc;
    const ControlPanel& mControlPanel;
    const CL& mCl;
    vector<Spline*> mCpsi_l;
    vector<Spline*> mSigma2;
    vector<Spline*> mC_gl2;

    typedef struct {unsigned int n; double value; } FactorStruct;
    //! Table of the few factorials needed for lensing
    static const FactorStruct factorial[12];

    /*! \struct LensingDifference
     *  \brief internal struct to keep track of lensing corrections
     *
     *  holds a spline with lensing corrections for each CMB spectrum
     */
    struct LensingDifference
    {
      LensingDifference( const int size );
      ~LensingDifference();
      Spline *ts, *cs;
      Spline *e, *b;
    };

    /*! \struct DlSplines
     *  \brief struct to compute and pass around the beta-term of the Wigner-D function
     *
     *  This struct is used by AllSkyLensing to compute the dlmns for each
     *  angle and pass them around.
     */
    struct DlSplines
    {
        //! compute the needed dlmns for angle theta up to l=lmax
        //! \warning call initializeFactorTables first!
        DlSplines( int lmax, double theta, Anchor* a = 0 );
        const Spline *dl00, *dl1m1, *dl2m2, *dl22, *dl31;
        const Spline *dl11, *dl3m3, *dl40, *dl4m4, *dl02;
        const Spline *dl3m1, *dl4m2, *dl20, *dlm24;

        //! fill the mFactor(1,2,3) maps that are used by dl();
        //! call this before constructing the first DlSplines
        static void initializeFactorTables( unsigned int lmax);

      private:
        /*! \brief compute dlmn of alpha
         *
         * calculate d^l_mn( beta ) up to l=lmax using the recurrence relations
         * from  Blanco, Florez, Bermejo.
         * This is the most expensive part of AllSkyLensing.
         *
         * @return pointer to Spline with Anchor a filled 
         * with x=2..lmax, y=d^2_mn(alpha)...d^lmax_mn(alpha).
         */
        inline const Spline * dl( int m, int n, double alpha, int lmax, Anchor* a = 0 );

        typedef pair<int, int> poi;
        typedef vector<poi> voip;

        //! internal helper function to fill mFactor1, mFactor2, mFactor3 maps
        static void addToVector( voip& v, int i, int j );

        //! Map that holds a vector for each mn-pair.
        //! The vector is filled with l(2l-1)/(sqrt(l^2-m^2)*sqrt(l^2-n^2)).
        static map< poi, vector<double> > mFactor1;

        //! Map that holds a vector for each mn-pair.
        //! The vector is filled with m*n*(2l-1)/(l-1)/(sqrt(l^2-m^2)*sqrt(l^2-n^2)).
        static map< poi, vector<double> > mFactor2;

        //! Map that holds a vector for each mn-pair.
        //! The vector is filled with (sqrt(l-1)^2-m^2)*sqrt((l-1)^2-n^2))/((l-1)*sqrt(l^2-m^2)*sqrt(l^2-n^2))
        static map< poi, vector<double> > mFactor3;
    };



   /* \struct Ximn
    * \brief small struct to precompute X_imn( beta ) and pass it around
    *
    * one Ximn struct for each angle is constructed in oneAngle() and passed
    * to the methods that need them. Redundant computation is avoided.
    */
   struct Ximn
   {
     /*! calculate X_imns for theta
      *
      *  the constructor fills the member vectors with the needed X_imns for theta from
      *  lstart to lend.
      */
     Ximn( const double theta, const double sigma2, const double c_gl2, const int lstart, const int lend );

     vector<double> X_000, X_220, X_022, X_121, X_132;
     vector<double> X_242, Xprime_000, Xprime_022;
   };
};
#endif //ALLSKYLENSING
