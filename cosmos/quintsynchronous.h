#ifndef QUINTSYNCHRONOUS_H
#define QUINTSYNCHRONOUS_H 

#include "synchronous.h"
#include "quintessence.h"

/*!
  Synchronous gauge Quintessence fluctuations. 
  Reimplements some functions and adds some
  for quintessence fields. 
  If you like to know about the time evolution of the perturbations,
  you may use the ofstream ofs defined here. Just put 
  a 
  \code
  #define OUTPUTFILE  
  \endcode
  somewhere in the quintsynchronous.h file
*/

class Quintessence;
class QuintCosmos;
class QuintSynchronous : public Synchronous {


 protected:
  QuintCosmos *quintcosmos;
  Quintessence *quint;
  int qidx; // the index position of the delta phi and delta phidot variable in the y[] array
 public:
  QuintSynchronous(QuintCosmos* c);
  virtual ~QuintSynchronous() {};

  virtual void fderivs(const double tau, const double *y, double *yprime);
  virtual void initialScalarPerturbations(const ControlPanel &control, const double tau);
 
  virtual void calcPerturbationNr(const ControlPanel &); //! calculate the number of scalar and tensor equations 
  
  virtual void getReady(const ControlPanel&);

  //double q(double a) { return quint->q(a); }   // field 
  //double qDot(double a) { return quint->qDot(a);}  // tau derivative of field
  //double Vprime() { return quint->Vprime();}  //  first derive of  potential w.r.t field 
  //double Vprime2() {return quint->Vprime2();}  // second deriv. of potential
};

#endif 
