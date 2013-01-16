#ifndef QUINTESSENCE_H
#define QUINTESSENCE_H 

#define Q_MV -9.341256734881e35

class QuintCosmos;

using namespace std;
#include <vector>
#include <string>
#include "gauge.h"
#include "global.h"

/*! 
  Small struct for giving explicit 
  information about the parameters. Needed
  for nice frontend.
*/
struct QPName {
  string name;
  string tooltip;
  string whatsthis;
  bool determined;
};

/*!
  Quintessence parent class.
  
  Provides basic functionality common to all quintessence realizations.
  
  Potential energy and initial conditions are pure virtual and have
  to be reimplemented in subclasses, such as Exponential, Ratra or voidQuintessence.
  
  In principle, construction a Quintessence field is easy:
  \code
  Exponential myQuint(myCosmos);
  \endcode
  If you specify some scalefactor a and would like the field (if it is one) to
  initialize according to your specs, use
  \code
  myQuint.setInitial(a);
  \endcode
  Now myQuint.q() and myQuint.qDot() will give you the field and 
  first derivative of the field with respect to conformal time

  The scale factor is usually not needed as most of the models will not
  depend on it anyway (the only use is in kinetic energy calculations).

  A call to V(phi, a) wil return the potential, a call to rho_kin(a) the
  kinetic energy.
  
  For convenience a stub V() is provided that returns the
  potential at scalefactor 1 (most models do not depend on it anyway)
  and field value q(), which you can set via
  \code
  myQuint.setQ(q);
  \endcode
  By accessing V(), q(), qdot() and delQ() etc, you will get values in powers
  of Mpc. Internally, however, the classes may calculate q in units
  of planc mass etc. 

*/
  
class Quintessence { 
 private:
  double Q, Qdot;
  double SpeedOfSound2; //!< Rest frame speed of sound^2, only for general fluid
 protected:
  QuintCosmos& cosmos;
  double Gpi8;
  double M_p;  
  vector<double> param;
 public:
  bool callFromHistory; // crossoverfield needs it 
  enum InitialConditions {adiabatic};
  enum Type {exponential, leaping, ipl ,none, arbitrary,celestine ,crossover,crossoverfield,corasaniti,tristar,omegastep,binomega, arbitraryomega,
             constantomega,cuscutan,coupledexp,coupledleaping,decelerationde,arbitraryfield};
 
  Quintessence(QuintCosmos&);
  virtual ~Quintessence() {};
  virtual void reset() {};  //!< try get into intial state, except for user defined settings (i.e.) make splines new etc.
  virtual void prepare() {}; //!< For Quintessence that do need preparation. 
  virtual void postPrepare() {}; //!< history() calls this at the end of the history calculation. Used by Arbitrary. 
  virtual Type type() { return none; }

  virtual double V(const double,const double=1) const = 0;                 //!< Potential
  virtual double Vprime(const double, const double=1,const double=1) const =0;        //!< dV/dphi
  virtual double Vprime2(const double,const double=1,const double=1,const double =1 ) const =0;      //!< d^2 V/ dphi^2

  // double V() const { return V(q(),1); }  //!< convenience using q()
  //double Vprime() const { return Vprime(q()); }  //!< convenience using q()
  //double Vprime2() const { return Vprime2(q()); }  //!< convenience using q()
  virtual double initialQ(const double) const =0;       //!< phi(a_initial )
  virtual double initialQDot(const double) const =0;  //!< d phi/ dtau (a_initial)

  virtual bool needsPropagation() { return true;} //!< If the field needs to be evolved for the background, this returns true. Arbitrary Quintessence, where one specifies w(a) does not need this. So it will return false. 
  virtual void setInitial(const double); //!< set intitial conditions

  virtual void touch(const double); //!< asks cosmos for the values of q and qdot at time tau

  virtual void setQ(const double q) { Q = q;}
  virtual void setQDot(const double qd) {Qdot = qd;}

  void setQs(const double q, const double qdot) { setQ(q); setQDot(qdot); }

  //void setPhi(const double phi) { setQ(phi);}          //!< convenience wrapper for setQ
  //void setPhiDot(const double phidot) { setQDot(phidot);}   //!< convenience wrapper for setQDot

  virtual double q(double a) const {return Q;}    //!< return field value Q which is stored in this Quintessence class
  virtual double qDot(double a) const {return Qdot;}  //!< return dQ/dtau 

  //double phi() const { return q();}   //!< convenience wrapper for q()
  //double phiDot() const {return qDot();}   //!< convenience wrapper for qDot()

  virtual double rho(const double a) const ;  //!< return energy density
  virtual double rho_kin(const double) const ;  //!< return energy density
  // double rho_pot() const { V(q()); }  //!< return energy density
  virtual double p(const double a) const ;  //!< return pressure
  virtual double w(const double a) const ;  //!< return equation of state  p = w * rho

  virtual double lambda() const {return param[0];} //!< return most prominent constant, param[0], i.e  in wetterich case: exp (-lambda * phi/ M_p) 

  virtual double getKinetic(double q) const { throw Bad_Error("Quintessence::getKinet() called, but not reimplemented in this subclass."); }

  virtual const vector<double> parameters() const { return param;}
  virtual vector< QPName > parameterNames() const { return vector<QPName>();}


  double speedOfSound2() { return SpeedOfSound2; }  //!< return c_s^2 of general dark energy fluid
  //! In the case of a general dark energy fluid, you can set the rest-frame speed of sound squared, i.e c_s^2 here
  void setSpeedOfSound2(double cs2) { SpeedOfSound2 = cs2; }

  // virtual void setLambda(const double l, const double l2=0) {}; //!< set parameter lambda. This is the most prominent one and for each model has a different meaning, l2 is the second most prominent one and so on :-)

  virtual void setParameters(const double p1, const double p2=Q_MV,const double p3=Q_MV,const double p4=Q_MV,const double p5=Q_MV); //!< built a vector and call setParameters(vector<double>). Avoid parameter having MagicValue Q_MV

  virtual void setParameters(const vector<double> & l); //!< copy to parameter list of this Quintessence. Only the first 'n' values are copied (if Quintessence needs n parameters).

  virtual void parameterResize();
  //!< if you specify *less* parameters than are needed, quintessence will call this member function. Now if the specific realization, like e.g. ratra, does not care about it (i.e. has some defaults or things like this), it should overwrite this behaviour to just return.



  virtual void tuneQuintessence(double omebar=0.0) {return;}

  virtual string name() const { return "Base";}

  //
  // Status Information
  //
  virtual void printStatus();

  virtual void getAlpha() { };
  virtual bool needsToUpdatePhi() { return false; } //!< true, if history should update Tau2Phi after postPrepare()



};

#endif 
