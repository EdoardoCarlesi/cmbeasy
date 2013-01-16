#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

#include "quintessence.h"

/*!
  Quintessence subclass for classical exponential
  potential.

  As a special feature, you may use

  \code
  setImitateConstant(true);
  \endcode

  This forces the initial values of q (and qdot, which is
  set to zero) to be chosen such that the Quintessence
  field contributes just as a ordinary cosmological
  constant.

  By default, this is feature is off.
*/

class Exponential : public Quintessence
{
 private:
  bool imitateConstant;
  double mV0;
  bool mHasInitialQ, mHasInitialQDot;
  double mInitialQ, mInitialQDot;
  bool mHasLambda;
  double mLambda;

 public:
  Exponential(QuintCosmos &c);

  virtual double V(const double,const double=1) const;                                        //!< Potential
  virtual double Vprime(const double,const double =1,const double=1) const;                  //!< dV/dphi
  virtual double Vprime2(const double,const double =1,const double=1,const double=1) const; //!< d^2 V/ dphi^2

  virtual double V0() const  { return mV0;}  //!< potential constant, Mp^4 by default

  // set the initial parameters individually
  void setInitialQ(const double);   //!< set the initial field value
  void setInitialQDot(const double); //!< set the initial value of dphi/dtau

  double initialQFromRadiationAttractor(const double a) const; //!< returns q on the rad-attractor at scale factor a

  void setLambda (const double);
  void setV0 (const double);

  virtual double initialQ(const double) const;     //!< phi(a_initial )
  virtual double initialQDot(const double) const;  //!< d phi/ dtau (a_initial) 
  virtual double lambda() const;                   //!< the potential is V0() * exp(-lambda phi/ M_p() ) where M_p is reduced planckmass, and V0 = Mp^4 by default

  //!< read the parameters and set them from the given vector
  //!< index 0 is lamda, 1 is initialQ, 2 is initialQDot, 3 is V0
  //!< "inf" as a value means use default values (i.e. on attractor/ v0= Mp^4)
  virtual void setParameters(const vector<double> &);

  void setImitateConstant(bool b) { imitateConstant = b;}

  virtual void printStatus();

  Type type() { return exponential;}
};


#endif 
