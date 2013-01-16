#ifndef ARBITRARY_H
#define ARBITRARY_H

#include "anchor.h"
#include "quintessence.h"
#include "mathobject.h"
#include <cmath>
#include "spline.h"

/*!
  Base class for quintessence models with arbitrary given equation of state w(a).
  Just re-implement (or change in here) the functions w(), dwda() and d2wda2().
*/

class Arbitrary : public Quintessence,  public Mathobject {
  double Mpinv2;
 protected:
  bool ValidPhi;
 public:
  Anchor anchor;
  Spline *Rho, *PhiDot, *Phi;
  Arbitrary(QuintCosmos&);

  /*
  virtual double w(double a) const  { return -2.0/3.0 -0.3*a*a;}  //!< equation of state as function of scale factor
  virtual double dwda(double a) const { return -0.6*a;} //!< first  derivative wrt a:  d w(a)/ da
  virtual double d2wda2(double a) const { return -0.6;}  //!< second derivative wrt a:  d^2 w(a)/ da^2
  */
  
  virtual double w(double a) const  { return -0.6;}  //!< equation of state as function of scale factor
  virtual double dwda(double a) const { return 0;} //!< first  derivative wrt a:  d w(a)/ da
  virtual double d2wda2(double a) const { return 0;}  //!< second derivative wrt a:  d^2 w(a)/ da^2


  double wdot(double a, double adotoa) const { return dwda(a)*adotoa*a;}
  double wdot2(double a, double adotoa,double addotoa) const { return a*(d2wda2(a)*adotoa*adotoa*a + dwda(a)*addotoa);}

 
  virtual void prepare();
  void  integrateWa(const double x, const double *y, double *dy) { dy[1] = 3*(1+w(x))/x; }
  //  void  integrate3Wa(const double x, const double *y, double *dy) { dy[1] = (1+3*w(x))/x; }
  



   double initialQ(double a) const { return q(a);}
   virtual double initialQDot(const double a) const  { return PhiDot->fastY(a);}


   virtual double rho(const double a) const { return Rho->fastY(a); }
   virtual double p(const double a) const { return w(a)*rho(a);}
   virtual double V(const double q,const double a =1) const { return 0.5*(1-w(a))*rho(a);}
   virtual double rho_kin(const double a) const { return 0.5*(1+w(a))*rho(a); }	
   virtual double q(double a) const;  //!< Careful, call reconstructPhi() first, (but after cosmos.history()) to get sensible values
   virtual double qDot(double a) const { return PhiDot->fastY(a);}
   void reconstructPhi(); //!< integrate phidot to obtain phi
   void postPrepare() { reconstructPhi(); }
	      
   virtual double Vprime2(const double q,const double a =1,const double adotoa=1,const double tau=1)  const; 
   double Vprime(const double q,const  double a=1,const double adotoa=1) const;

 
   virtual bool needsPropagation() { return false;}
   virtual bool needsToUpdatePhi() { return true; } //!< true, if history should update Tau2Phi after postPrepare()
  virtual Type type() { return arbitrary;}

};


#endif 
