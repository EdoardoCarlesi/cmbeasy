#ifndef COUPLEDEXP_H
#define COUPLEDEXP_H

#include "exponential.h"
#include "coupledquintcosmos.h"

/*!
  Quintessence subclass for a coupled exponential
  potential.
*/

class CoupledExp : public Exponential {

 private:
  bool imitateConstant;
  double mQinitial, mQinitialFactor;
  bool mHasQInitial; //!< true if QInitial has been explicitely set to some value
  bool mHasQInitialFactor; //!< true if QInitialFactor has been explicitely set to some value


 public:
  CoupledExp(CoupledQuintCosmos&);
   
  virtual double lambda() const;  //!< the potential is M_p()^4 * exp(-lambda phi/ M_p() ) where M_p is reduced plankm.

  virtual void printStatus(); 

  virtual void   setInitialQ(const double qInitial) { mHasQInitial = true;  mQinitial = qInitial; }
  virtual void   setInitialQFactor(const double f) { mQinitialFactor = f; mHasQInitialFactor = (f!=1.); }
  virtual double initialQ(const double) const;       //!< phi(a_initial )
  virtual double initialQDot(const double) const;  //!< d phi/ dtau (a_initial)

  Type type() { return coupledexp;}
  double B;
};


#endif 
