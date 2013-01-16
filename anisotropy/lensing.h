#ifndef LENSING
#define LENSING

class Cosmos;
class CL;
class ControlPanel;
class CmbCalc;

#include <string>

/*! \class Lensing
 *  \brief computes the lensing effect
 *
 *  The Lensing class computes the lensing effect on the CMB spectra using using one of the
 *  available lensing methods @see LensingMethod
 *
 */


class Lensing
{
  public:
    /*! \brief Constructs a Lensing object
     *
     *  @param cosmos is a Cosmos object that provides the lensing potential power spectrum
     *  @param cp a ControlPanel holding the ( lensing ) settings
     *  @param cl the unlensed spectra
     *  @param cmbCalc the CmbCalc object that was used to compute the cl
     */
    Lensing( Cosmos& cosmos, const ControlPanel& cp, CL& cl, const CmbCalc& cmbCalc );

    /*! \brief sets mBesselFileName
     *
     * mBesselFileName is the name of the file that contains the precalculated bessel functions.
     * You do not need to call this method if you do not use the "flatsky" 
     * lensing method.
     *
     * @param name the name (and path, if not in the current working directory) of the file that contains the precalculated bessel functions
     */
     void setBesselFileName( const std::string& name ) { mBesselFileName = name; }

    /*! The available lensing methods
     *
     *  the only available method now is  the full-sky-method of astro-ph/0502425
     *  ('FullSky'); the original CMBFAST lensing ('FlatSky') is no longer part
     *  of the official cmbeasy distribution
     *
     *  The full-sky-correction-function method is implemented by the class AllSkyLensing.
     *
     *  @see AllSkyLensing
     *  @see dumpLensedCls
     *
     */
    enum LensingMethod { Automatic, FlatSky, AllSky };

    /*! \brief lense the CL spectra handed to the constructor
     *
     *  The lense method computes the lensed spectra.
     *  This method directly modifies the CL reference passed
     *  to the constructor by overwriting the unlensed spectra
     *  with the lensed versions.
     *  
     *  @see LensingMethod
     */
    void doLensing( const LensingMethod method = Lensing::Automatic );

    /*! \brief compute the lensed spectra
     *
     *  The lensedCls method computes the lensed spectra.
     *  Unlike the above method, it leaves the CLs passed to the
     *  constructor untouched.
     *
     *  @see LensingMethod
     *  @returns a pointer to a CL object holding the lensed spectra
     */
    CL* lensedCls( const LensingMethod method = Lensing::Automatic );


    /*! \brief compute and write the lensed spectra to disc
     *
     * This functions computes the lensed spectra and writes them to the file fileName.
     * It leaves the CL reference passed to the constructor untouched.
     * For available lensing methods: @see LensingMethod
     * 
     * @param fileName the file to write the lensed spectra to
     * @param method the method used to compute the lensing effect;
     * if set to Lensing::Automatic, the method will be determined
     * from the ControlPanel passed to the constructor.
     */
    void dumpLensedCls( const std::string& fileName,  const LensingMethod method = Lensing::Automatic );

  protected:

    /* \brief return the lensing method used
     *
     * methodUsed determines the lensing method to use
     * if the LensingMethod method is passed to one of
     * the functions lensed(), lensedCls() or dumpLensedCls().
     * The returned value is only different from the the parameter
     * m if m is Lensing::Automatic.
     *
     * @param m the lensing method, can be Lensing::Automatic
     * @returns the lensing method, is guaranteed to not equal Lensing::Automatic
     * @see LensingMethod
     */
     LensingMethod methodUsed( LensingMethod m );



  private:
     Cosmos& mCosmos;                     //!< the Cosmos object passed to the constructor
     const ControlPanel& mControlPanel;   //!< the ControlPanel object passed to the constructor
     CL& mCl;                       //!< the unlensed Cls
     const CmbCalc& mCmbCalc;             //!< the CmbCalc object passed to the constructor
     std::string mBesselFileName;   //!< the name of the file containing the precalculated bessel functions
};
#endif //LENSING
