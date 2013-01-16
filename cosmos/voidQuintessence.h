#ifndef VOIDQUINTESSENCE_H
#define VOIDQUINTESSENCE_H 

#include "quintessence.h"
#include <iostream>
/*!
  Quintessence class that has one very nice property: Whenever 
  you ask for something, you get 0 back !
  
  This allows elegant case independent programming in the
  cosmos class:

  e.g. totalRho = rho_m + rho_r + ... + quintessence->rho() 

  And if this quintessence is a VoidQuintessence, well then you will get back 0 all
  the way 
*/
  
class VoidQuintessence : public Quintessence {
 public:
  VoidQuintessence(QuintCosmos &c) : Quintessence(c)  {};

  double V(const double q,const double a=1) const {return 0;}
  double Vprime(const double, const double a=1,const double=1)  const {return 0;}  //!< dV/dphi
  double Vprime2(const double, const double a=1,const double adotoa=1,const double tau=1) const {return 0;}

  double initialQ(const double)   const   {return 0;}   //!< phi(a_initial )
  double initialQDot(const double)  const {return 0;}  //!< d phi/ dtau (a_initial)

   double p(const double a) const { return 0;}  //!< return pressure
   double w(const double a) const { return 0;}
   double rho_kin(const double) const { return 0;}

   void setQ(const double q) {};
   void setQDot(const double qd) {};

   double delQ() const {return 0.0;}
   double delQDot() const {return 0.0;}

   void printStatus()  { cout << "This is a VoidQuintessence\n"; Quintessence::printStatus(); }

};
#endif 
