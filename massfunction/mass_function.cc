/* mass_function.cc */
#include <iostream>
#include <math.h>
#include "spline.h"
#include "anchor.h"
#include "mass_function.h"

#define PI 3.14159265
#define inv_log_10 1./log(10)
// Size of the most used mass function splines
#define SIZE 1000

MassFunction::MassFunction(Spline* pk, double mmin, double mmax, double rho){
set_pk(pk);
set_m_min(mmin);
set_m_max(mmax);
set_rho0(rho);
}

double MassFunction::W_k(double k, double R){
return 3*(sin(k*R) - k*R*cos(k*R))*pow(k*R,-3);
}

double MassFunction::R_m(double m){
return pow((3*m)/(rho_0*4*PI), 1./3.);
}

double MassFunction::M_r(double r){
return pow(r,3)*rho_0*4*PI*(1./3.);
}

double MassFunction::norm(double M){
double lnslnm = sqrt(pow(dLnSigma_dLnM->fastY(log(M)),2)); 
double norm = rho_0/(M*M);
return norm*lnslnm*sqrt(2./PI);
}

double MassFunction::sigma8(){
double m = M_r(8);
//fprintf(stderr, "sigma8, M: %e, inv_ln_10: %lf\n", m, inv_log_10);
return Sigma->fastY(m);
}

void MassFunction::init(double norm){
//cout << "Initializing and generating mass function..." << endl;
initSplines();
fill();
normalize_sigma8(norm);
init_standard_param();
//cout << "End init() " << endl;
set_differential();
set_integral();
if(z==0.) {
cout << "Sigma_8: " << sigma8() << endl;
}}

void MassFunction::normalize_sigma8(double norm){
if(!normalization()) norm=1.;
for (int i=0; i<Sigma->size(); i++){
Sigma->mul(i,norm);
dLnSigma_dLnM->mul(i,norm);
}}

void MassFunction::initSplines(){
//cout << "MassFunction::init()" << endl; 
Sigma = new Spline(SIZE, "sigma", &SplineAnchor);
LnSigma_LnM = new Spline(SIZE, "ln_sigma", &SplineAnchor);
dLnSigma_dLnM = new Spline(SIZE, "dLnSigma_dLnM", &SplineAnchor);
dLnSigma_M = new Spline(SIZE, "dLnSigma_M", &SplineAnchor);
}

void MassFunction::fill(){
/* First fill all sigma(M) splines */
//cout << "MassFunction::fill()" << endl; 
double k, k_max, k_min;
double m;
int size_k = SIZE; 
int size_m = 200;
k_min = Pk->start(); k_max = Pk->stop();
//k_min = 0.1; //Pk->start(); 
//k_max = 1; //Pk->stop();
Pk->arm();

double step_k = log(k_max/k_min)/(size_k-1);
double step_m = log(M_max()/M_min())/(size_m-1);
// Fill for the masses
for(int j=0; j<size_m; j++){
Spline *integrand_sigma = new Spline(SIZE, "integrand_sigma", &SplineAnchor);
m = m_min*exp(step_m*j);
double r=R_m(m);
//cout << "R: " << r << " M: " << m << endl;
// Integrate on the k_s
for(int i=0; i<size_k; i++){
k = k_min*exp(step_k*i);
double i_s = k*k*Pk->fastY(k)*pow(W_k(k,r),2);
//cout << "MassFunction::fill() - k: " << k << " Pk: " << Pk->fastY(k) << " integrand: " << i_s << endl;
integrand_sigma->set(k,i_s);
}
//cout << "MassFunction::fill()" << endl; 
integrand_sigma->arm();
double i_i_s = integrand_sigma->integrate(k_min,k_max);
i_i_s /= 2*PI*PI;
//cout << "MassFunction::fill(). M: " << m << " R: " << R_m(m) <<  " integrand_sigma: " << sqrt(i_i_s) << endl;
double sig = sqrt(i_i_s);
Sigma->set(m,sig);
LnSigma_LnM->set(log(m),log(sig));
delete integrand_sigma;
}
Sigma->arm();
LnSigma_LnM->arm();
LnSigma_LnM->derive(*dLnSigma_dLnM);
dLnSigma_dLnM->arm();
}

double MassFunction::get_integral(double M){
return Integral->fastY(M);
}

double MassFunction::integral_number_density(double M){
return Differential->integrate(M,M_max());
}

void MassFunction::set_integral(){
//cout << "MassFunction::set_integral() " << endl;
double mmin = M_min(); double mmax = M_max();
int size=SIZE; double m, ind;
Integral = new Spline(size, "integral", anchor());
double step = log(mmax/mmin)/(size-1);
for(int i=0; i<size-1; i++){
m=mmin*exp(step*i);
ind=integral_number_density(m);
Integral->set(m, ind);
}
//cout << "MassFunction::set_integral() " << endl;
Integral->arm();
}

void MassFunction::set_differential(){
double mmin = M_min(); double mmax = M_max();
int size=SIZE; double m, dnd;
Differential = new Spline(size, "differential", anchor());
double step = log(mmax/mmin)/(size-1);
for(int i=0; i<size; i++){
m=mmin*exp(step*i);
dnd=differential_number_density(m);
Differential->set(m, dnd);
}
Differential->arm();
}

Tinker::Tinker(Spline *power, double mmin, double mmax, double rho) : MassFunction(power, mmin, mmax, rho){
// Constructor is the same as the one in MassFunction
;
}

void Tinker::init_standard_param(){
	p[0]=0.186; p[1]=1.47;
	p[2]=2.570; p[3]=1.19;
}

double Tinker::differential_number_density(double M){
double N = norm(M); 
double s = sigma(M);
//cout << "M: " << M << " Norm: " << N << " sigma: " << s << endl;
return N*p[0]*(pow(s/p[2], -p[1])+1)*exp(-p[3]/(s*s));
}

TinkerZ::TinkerZ(Spline *power, double mmin, double mmax, double rho) : Tinker(power, mmin, mmax, rho){
// Constructor is the same as the one in Tinker
;
}

double TinkerZ::differential_number_density(double z, double M){
double N = norm(M); 
double s = sigma(M);
//cout << "Norm: " << N << " sigma: " << s << endl;
return N*param(0,z)*(pow(s/param(2,z), -param(1,z))+1)*exp(-param(3,z)/(s*s));
}

double TinkerZ::param(int i, double z){
double pp[4]; double delta = 200.;
double alpha = exp(-pow(log(0.75/(delta/75)), 1.2));
pp[0] = p[0]*pow(1+z, -0.14);
pp[1] = p[1]*pow(1+z, -0.06);
pp[2] = p[2]*pow(1+z, alpha);
pp[3] = p[3];
return pp[i];
}

void TinkerZ::set_differential(double z){
double mmin = M_min(); double mmax = M_max();
int size=SIZE; double m, dnd;
Differential = new Spline(size, "differential", anchor());
double step = log(mmax/mmin)/(size-1);
for(int i=0; i<size; i++){
m=mmin*exp(step*i);
dnd=differential_number_density(z,m);
Differential->set(m, dnd);
}
Differential->arm();
}

ShethTormen::ShethTormen(Spline *power, double mmin, double mmax, double rho) : MassFunction(power, mmin, mmax, rho){
;
}

void ShethTormen::init_standard_param(){
// TODO
	p[0]=0.186; p[1]=1.47;
	p[2]=2.570;
}

double ShethTormen::differential_number_density(double M){
//TODO
double N = norm(M); 
double s = sigma(M);
return N*p[0];
}
