#ifndef RECOMBINATION_H
#define RECOMBINATION_H

#include "miscmath.h"
#include "spline.h"
#include "anchor.h"

class Cosmos;
class ControlPanel;

/*!
  Base class for all recombination versions. It defines the 
  interface for all recombination routines. The functions
  include the free electron fraction, the baryon temperature,
  the sound speed of baryons and the value of the
  fine structure constant. Please note that the absolute
  value of the fine structure constant does not play a role,
  only the relative change with redshift will be used by 
*/
class Recombination : public AnchorEnabled, public Mathobject {
 public:
  enum recombination { Recfast, RecfastAlpha,RecfastQuint};
  Spline *noteX, *noteTb , *noteTbDeriv,*noteCs2 ;
  Anchor convenient;
  Cosmos &cosmos;

  Recombination(Cosmos& c, ControlPanel& control, Anchor* a=0)
      : AnchorEnabled(a), noteX(0), noteTb(0), cosmos(c), mControl(control) {
    noteX = new Spline(10000,"Recombination x",&convenient);
    noteTb = new Spline(noteX, "Recombination Tb",&convenient);
    noteTbDeriv = new Spline(noteX, "Recombination TbDeriv",&convenient);
    noteCs2 = new Spline(noteX, "Recombination cs2",&convenient);
  };
  virtual void compute()=0;
  

  // please note that all functions are encoded as functions of -z, i.e. 
  // baryon temperature at z =3 is Tb(-3);
  double xe(double z) { return noteX->fastY(z); } //!< return free electron fraction as a function of -z
  double Tb(double z) { return noteTb->fastY(z);} //!< return baryon temperature as a function of -z
  double TbDeriv(double z) { return noteTbDeriv->fastY(z);}  //  d (ln T_baryon) / d (ln a) as a function of -z
  double cs2(double z) { return noteCs2->fastY(z);} //!< return baryon sound speed as a function of -z
  virtual double BoltzmannOverWeightAndC2()=0; //!< see the implementation
  double maxRedshift() {return - noteX->start(); } //!< maximum redshift covered by recombination
  double highZFraction() { return noteX->front(); } //!< x_e at highest redshift
  virtual double fineStructureConstant(double z) { return 1/137.0;}  //!< the fine structure constant, not needed except for varying alpha
  virtual double delAlpha(double z) { return fineStructureConstant(z)/fineStructureConstant(0) - 1.0; } // relative change in alpha

  Spline* delAphaSpline(Anchor* a,double zmax=100,int steps=1000) {
    Spline *s = new Spline(steps,"delAlpha",a);
    double logmax = log(zmax+1);
    double logstep=logmax/steps;
    for (double logz = 0; logz <= logmax+1e-7; logz += logstep) {
      double z = exp(logz) - 1.0;
      s->set(z,delAlpha(z));
    }
    return s;
  }

 protected:
  ControlPanel& mControl;

};

#endif 
