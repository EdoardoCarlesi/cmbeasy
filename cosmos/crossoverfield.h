#ifndef CROSSOVERFIELD_H
#define CROSSOVERFIELD_H

#include "arthur.h"
#include "mathobject.h"
// #include "spline.h"

class Spline;

/*! 
  Crossover quintessence class. The kinetic
  term changes as a function of the field. This
  running is determined  by a renomralization
  group running of the kinetic term. 
*/
class CrossoverField : public Arthur {

 
 protected:
  Spline *delta;
  bool Tuning;

 public:  
  
  CrossoverField(QuintCosmos&);
 
  virtual Type type() { return crossoverfield; }
  virtual double estimatedRhoQ(const double a) const;
  virtual double initialQ(const double) const;       //!< phi(a_initial )
  virtual double initialQDot(const double) const;  //!< d phi/ dtau (a_initial)
  
  virtual void prepare();
  virtual void reset();

  virtual double kineticTerm(const double q) const;
  
  virtual void flow(const double, const double *, double*);
  void flowWrapper(const double x, const double *y, double*dy) { flow(x,y,dy);} 
  // virtual void printStatus();

  virtual void tuneQuintessence(double omebar=0.0);
  //virtual vector<QPName> parameterNames() const;
  virtual string name() const { return "CrossoverField";}

  void getAlpha();
  void betaFunction(const double, const double *, double*);

};

struct PhiCritTooLow {};

#endif 
