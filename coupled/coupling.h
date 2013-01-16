#ifndef COUPLING_H
#define COUPLING_H

#include <vector>
#include <iostream>

#include "quintessence.h"

/*!\class Coupling
 * \brief Base class for couplings of dark energy and dark matter.
 *
 * Pure virtual base class for specifying the coupling between
 * quintessence and dark matter.
*/
class Coupling {
 public:
  Coupling(): mQuint(0) {}
  virtual ~Coupling() {}

  /*! the different types of couplings */
  enum Type {none,constant,exponential, varexponential};

  /*! This may be used to set the parameters of the coupling.*/
  virtual void setCouplingParameters(const std::vector<double>&) =0;

  /*! Returns the value of the source term
   * \f$ S \equiv < m_{, \phi} \bar{\psi}\psi > \f$
   * at the field value phi for this coupling.
   */
  virtual double phi2S(const double phi, const double rho_c) const = 0;

  /*! Returns the dS/dtau.
   *  @see phi2S()
   */
  virtual double phi2Sprime(const double phi, const double phiprime,
                            const double rho_c, const double adotoa) const = 0;

  /*! Returns the value of the gauge invariant perturbation
   *  for the perturbed source
   * \f$ \tilde{S} \f$
   */
  virtual double phi2DS(const double phi, const double phiprime, const double rho_c,
                        const double Delta_c, const double adotoa, const double Phi,
                        const double X) const = 0;

  /*! Returns the source Q_phi_0(phi, phiprime, rho_c) of the conservation equation
   *  for phi, \f$ \nabla_\nu T_{(\phi) \mu}^\nu = Q_{(\phi)\mu}  \f$.
   *  phiprime is dphi/dtau.
   */
  virtual double phi2Q_phi_0 (const double phi, const double phiprime, const double rho_c) const = 0;

  /*! Returns tau derivative of Q_phi_0.
   */
  virtual double phi2Q_phi_0prime(const double phi, const double phiprime, const double rho_c,
                                  const double a, const double adotoa) const {
                                    const double S = phi2S(phi, rho_c);
                                    const double Sprime = phi2Sprime(phi, phiprime, rho_c, adotoa);
                                    const double Q_phi_0_prime = 3.*adotoa/a*S*phiprime
                                                                 + a*S*S
                                                                 - Sprime*phiprime/a
                                                                 + a*S*mQuint->Vprime(phi, a, adotoa);
                                    return Q_phi_0_prime;
                                  }

  /*! Returns tau derivative of Q_c_0.
   */
  virtual double phi2Q_c_0prime(const double phi, const double phiprime, const double rho_c,
                                  const double a, const double adotoa) const {
                                    return -phi2Q_phi_0prime(phi, phiprime, rho_c, a, adotoa);
                                  }

  /*! Returns the source Q_c_0 = - Q_phi_0 of the conservation equation
   *  for cdm, \f$ \nabla_\nu T_{(c) \mu}^\nu = Q_{(c)\mu}  \f$.
   */
  virtual double phi2Q_c_0(const double phi, const double phiprime, const double rho_c) const {
                         return -phi2Q_phi_0(phi, phiprime, rho_c); }


  /*! returns a vector of the parameters of the coupling */
  virtual std::vector<double> getParameters() const = 0;

  /*! returns the type of the coupling */
  virtual Type type() = 0;

  /*!  out status information */
  virtual void printStatus(std::ostream& out=std::cout) = 0;

  /*! Set the instance of a coupled quintessence class
   * this coupling is operating on. */
  virtual void setQuintessence(Quintessence *q) {mQuint=q;}

  /*! This method is called at the end of CoupledQuintCosmos::history.
   *  Overwrite it if your coupling needs to be updated after the background
   *  evolution is known.
   */
  virtual void postPrepare() {}

 protected:
  Quintessence *mQuint;  //!< holds a pointer to the quintessence class that is using this coupling.
};

#endif
