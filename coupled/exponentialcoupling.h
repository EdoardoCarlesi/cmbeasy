#ifndef EXPONENTIALCOUPLING_H
#define EXPONENTIALCOUPLING_H

#include "coupling.h"
#include "quintessence.h"
#include "spline.h"
#include "basecosmos.h"

// debugging only
#include <iostream>

class Anchor;

/*!This class implements an exponential coupling of the form
 * \f[
 *    m(\chi) = m_0 \exp{-\frac{b \chi}{M}}
 * \f]
 */
class ExponentialCoupling : public Coupling
{
 public:
  ExponentialCoupling();  //!< constructor
  virtual ~ExponentialCoupling(); //!< destructor

  void setCouplingParameters(const vector<double>&);
  void setCouplingParameters(double beta);
  /*! Returns the effective value of the coupling
   *  for the given value of the field phi.
   *  In the case of leaping kinetic term quintessence,
   *  this is
   *  \f$ \beta(\phi) \equiv \left< \frac{b}{k(\chi(\phi))} \right> \f$.
   */

  // Now it returns a CONSTANT Beta - commented the phi-dependent part to avoid comp. errors!
  double phi2beta(const double phi) const { 
//cout << " mb: " << mb << endl;
	  return mb; }

  /*! Returns the derivate of beta with regard to the field.
   *  @see phi2beta()
   *
   *  @warning: must only be called after postPrepare() (which is
   *            usually called at the end of CoupledQuintCosmos::history() )
   */
  double phi2betaPrime(const double phi) const {
                return (mbetaPrimeSpline?mbetaPrimeSpline->fastY(phi):0.); }

  /*! Returns the second derivate of beta with regard to the field.
   *  @see phi2betaPrime()
   *
   *  @warning: must only be called after postPrepare() (which is
   *            usually called at the end of CoupledQuintCosmos::history() )
   */
  double phi2betaDoublePrime(const double phi) const {
                return (mbetaDoublePrimeSpline?mbetaDoublePrimeSpline->fastY(phi):0.); }

  /*! Returns the source term S,
   *  \f$ S(\phi, \rho_c) \equiv \frac{\beta(\phi)}{M} \rho_c \f$.
   *  for the given value of the field phi and the density of cold
   *  dark matter.
   */
  virtual double phi2S(const double phi, const double rho_c) const {
                return -phi2beta(phi)*rho_c/ baseCosmos::M_p(); 
  }

  /*! Returns the tau derivative of S.
   */
  virtual double phi2Sprime(const double phi, const double phiprime,
                            const double rho_c, const double adotoa) const {
                              const double beta = phi2beta(phi);
                              double Sprime = -phi2betaPrime(phi) * phiprime * rho_c;
                              Sprime += 3.*beta*adotoa*rho_c
                                       +beta *phi2Q_c_0(phi, phiprime, rho_c);
                             Sprime /= baseCosmos::M_p();
                              return Sprime;
                            }

  /*! Returns the value of the gauge invariant perturbation
   *  for the perturbed source
   * \f$ DS = -\frac{\beta (\phi)}{M} \rho_c \left[ \Delta_c
   *          - 3(1+\frac{\beta (\phi) \phi \prime )}{ 3 M  \mathcal{H}} \Phi  \right]
   *          - \frac{\beta_{,\phi}}{M} \rho_c X \f$
   */
  virtual double phi2DS(const double phi, const double phiprime, const double rho_c,
                        const double Delta_c, const double adotoa, const double Phi,
                        const double X) const {
                          const double M_p = baseCosmos::M_p();
                          const double beta = phi2beta(phi);
                          const double betaprime = phi2betaPrime(phi);
                          double DS = Delta_c - 3.*(1.+beta*phiprime/(3.*M_p*adotoa))*Phi;
                          DS *= -beta * rho_c / M_p;
                          DS -= betaprime*rho_c*X/M_p;
                          return DS;
                        }

  virtual double phi2Q_phi_0(const double phi, const double phidot, const double rho_c) const {
                                                           return phi2S(phi, rho_c)*phidot; }

  /*! Returns the original constant b
   *  in the exponential of the coupling. The method phi2beta(), on
   *  the other hand, returns the effective coupling that is
   *  actually used. This might differ from the value returned
   *  by this  function e.g. for quintessence models with
   *  a non-standard kinetic term after rewriting to a canonical term
   *  with a different potential.
   *
   *  @see betaSpline()
   *  @see postPrepare()
   */
  double b_original() const {return mb;}

  /*! Returns a pointer to a spline holding the effective value of the coupling
   *  as a function of phi.
   */
  Spline* betaSpline() const {return mbetaSpline;}

  /*! Returns the derivative of betaSpline() with respect to the field
   */
  Spline* betaPrimeSpline() const {return mbetaPrimeSpline;}

  /*! Returns the derivative of betaPrimeSpline() with respect to the field
   */
  Spline* betaDoublePrimeSpline() const {return mbetaDoublePrimeSpline;}


  /*! Updates the beta-related splines w/ the knowldedge of the background evolution.
   *  Automatically called from CoupledQuintCosmos::history().
   */
  virtual void postPrepare();

  /*! Returns a vector of the parameters for this coupling.
   *  Since there is only one parameter, the first and only entry
   *  of the returned vector will be equal to b_original().
   */
  virtual vector<double> getParameters() const;

  /*! Returns the type of this coupling, i.e. exponential */
  virtual Type type() {return exponential;}

  /*!  out status information */
  virtual void printStatus(std::ostream& out=std::cout) {
    out << "ExponentialCoupling, coupling set to: " << mb << endl;
  }


 private:
  double mb;  //!< the coupling constant b. Can be read by calling b_original().
  Spline* mbetaSpline;   //!< holds the effective coupling parameter beta as a function of the field.
  Spline* mbetaPrimeSpline;   //!< dbeta/dphi
  Spline* mbetaDoublePrimeSpline;   //!< d^2 beta/d phi^2
};

#endif // EXPONENTIALCOUPLING_H
