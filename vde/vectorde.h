#ifndef VECTORDE_H
#define VECTORDE_H 

#define Q_MV -9.341256734881e35

class VdeCosmos;

using namespace std;
#include <vector>
#include <string>
#include "gauge.h"
#include "grid.h"
#include "global.h"
#include "vdecosmos.h" 

// DEBUG
#include <iostream>

class Vde { 
private:
 /* Initial values for the vector field */
  double A0, rho_A0;
 /* Today's value for the energy fraction of VDE and matter */
 double Omega_A0;
 double Omega_M0;
 double Omega_k0;
 /* General value for the energy density */
  double rho_A;
  double SpeedOfSound2;
 double ww;

/* Store the omega, z and w in a grid to allow for interpolation */
 GridXYZ *grid;

protected:
  VdeCosmos& cosmos;
  double Gpi8;
  double M_p;  
  double rho_c; 
  double H_0;
public:
  enum InitialConditions {adiabatic};
 
  Vde(VdeCosmos&);
  ~Vde(){
delete grid;
};
  virtual void reset() {
grid->free_mem();
}  

 double initialRho(double);
 void setInitialA0();
//static double w_integ(double, void*);

void initialize_grid(string file, int rep, double w0);

  void set_rho(double v){rho_A=v;}

  double rho();  //!< return energy density
  double p(double);  //!< return pressure
  double w(double);
  void setGrid(GridXYZ *g){grid=g;}
  void set_rhoA0(double a){rho_A0=a;}
//!< return equation of state  p = w * rho
  void set_OmegaA0(double oA){Omega_A0=oA;}
  double get_OmegaA0(){return Omega_A0;}
  double getInitialA0();
  void set_OmegaM0(double oM){
//cout << "Vde:set_OmegaM0()" << oM << endl;	
Omega_M0=oM;}
  double get_OmegaM0(){return Omega_M0;}
  double get_Omegak0(){return Omega_k0;}
  void set_H_0(double hi){H_0=hi;}
  void setW(double wa){ww=wa;}
  double speedOfSound2(){return SpeedOfSound2;}  //!< return c_s^2 of general dark energy fluid
  //! In the case of a general dark energy fluid, you can set the rest-frame speed of sound squared, i.e c_s^2 here
  void setSpeedOfSound2(double cs2){SpeedOfSound2 = cs2;}
  // Status Information
  virtual void printStatus();
};

#endif 
