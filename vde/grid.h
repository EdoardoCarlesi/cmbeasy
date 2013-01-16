/* grid.h */
/* File containing different routines for i/o output and other stuff */
#ifndef gridxyz_H
#define gridxyz_H
#include <iostream>
#include <fstream>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

class GridXYZ {

/* File containing the table */
string fName;

/* Internal variables used for the interpolation */
double z_0,z_1,z00,z01,z10,z11;
double x0,x1,y0,y1;

/* Eps=x spacing, repeat=number of position after which the y is repeated */
double eps;
double w0;
int repeat;
int size;

/* Arrays storing the grid points */

double *x;
double *y;
double *z;

//gsl_vector *vx, *vy, *vz;

public:

GridXYZ();
~GridXYZ();

void set_file_name(string f){
fName=f;
//cout << "GridXYZ::set_file_name(). Setting file name to: " << f << endl;
}

void set_eps(double e){
eps=e;
//cout << "GridXYZ::set_eps(). Setting x step to: " << 3 << endl;
}

void set_repeat(int r){
repeat=r;
//cout << "GridXYZ::set_repeat(). Setting number of z repeating to: " << r << endl;
}

void set_w0(double w){
w0=w;
//cout << "GridXYZ::set_w0(). Setting w0 to: " << w << endl;
}

string get_file_name(){
return fName;
}

void initialize_grid();

void free_mem(){

delete [] x;
delete [] y;
delete [] z;
/*
gsl_vector_free(vx);
gsl_vector_free(vy);
gsl_vector_free(vz);
*/
}

/* Given x and y interpolate and get the corresponding z value */
double interpolate_z(double, double);

};
#endif
