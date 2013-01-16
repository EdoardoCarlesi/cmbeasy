#include "crossover.h"  


Crossover::Crossover(QuintCosmos &c) : Arbitrary(c) , a_catch(0.001)  {
    param.resize(3);
    param[0] = -0.91;
    param[1] = -0.188;
    param[2] = -0.0091;
    

   
}

void Crossover::prepare() {
  CrossAnchor.kill();
  CW = new Spline(1000,"Crossover",&CrossAnchor);
  double step=0.1;
  cout << "crossover: prepare" << endl;
  for (double loga = -16; loga <= 0.2; loga += step) {
    double a = pow(10,loga);
    cout << a << "    " << CrossW(a,false) << endl;
    CW->set(a,CrossW(a,false));
  }
  CW->arm();


  a_catch = CW->maximum();
  if (a_catch < 1e-5 || a_catch > 1e-2) throw Bad_Error("Crossover::prepare() strange a_catch");
  
  cout << "THE MAXIMUM IS AT: " << a_catch << endl;


  CrossAnchor.kill();
  CW = new Spline(1000,"Crossover",&CrossAnchor); 
  CW1 = new Spline(CW, "Cw1",&CrossAnchor);
  CW2 = new Spline(CW,"Cw2",&CrossAnchor);
  
  step = 0.1;
  for (double loga = -18; loga <= 0.2; loga += step) {
    double a = pow(10,loga);
    cout << a << "    " << CrossW(a,true) << endl;
    CW->set(a,CrossW(a,true));
  }
  CW->arm();


  CW->derive(*CW1);
  for (int n=0; n < CW1->size(); n++) {
    if (CW1->x(n) < a_catch) CW1->setY(n,0);
  }
  CW1->arm();

  CW1->derive(*CW2,Spline::yes);

  CW->dump("crossw",false);
  CW1->dump("crossw1",false);
  CW2->dump("crossw2",false);

  Arbitrary::prepare();
}


  
 
