#include "vectorde.h"
#include "grid.h"
#include "global.h"
#include "vdecosmos.h"
#include <math.h>
#include <gsl/gsl_integration.h>
#include <iostream>

Vde::Vde(VdeCosmos &c) : SpeedOfSound2(1) , cosmos(c)
{
  Gpi8 = cosmos.Gpi8();
  M_p = cosmos.M_p();
  rho_c = cosmos.rho_0();
}

/* Interpolate from a table */
double Vde::w(double a) { 
double zp1 = 1./a;
double om = get_OmegaM0(); 
double ok = get_Omegak0();
//cout << "Omega_M0 " << get_OmegaM0();
//Omega_M0=0.3;
ww = grid->interpolate_z(om, zp1);
//if(a==1.) cout << "Vde->w(a): w, om, zp1 " << ww << " " << om << ", " << zp1 << endl;
//ww = -0.8;

return ww;
}

double Vde::rho(){
	return rho_A;
}

double Vde::p(double a){
	return w(a)*rho();
} 

void Vde::initialize_grid(string file, int rep, double w0){
grid = new GridXYZ; 
grid->set_file_name(file);
grid->set_repeat(rep);
grid->set_w0(w0);

/* Actually read the values and load them */
grid->initialize_grid(); 
}

void Vde::setInitialA0(){
A0=(1./H_0*sqrt(2*rho_A0/3));
// cout << "Vde::initialA0(), rho_A0, H_0: " << rho_A0 << " " << cosmos.H_0_cpm() << endl; 
}

double Vde::getInitialA0(){
return A0/M_p;
}

void Vde::printStatus() {
  cout << "======== Vde status =======\n";
  cout << "rho_A0() : " << rho_A0 << endl;
  cout << "A0       : " << getInitialA0()*pow(10,4) << " 10^-4 M_p" << endl;
  cout << "w(a=1)   : " << w(1.) << endl;
  cout << "\n=========================\n";
}
