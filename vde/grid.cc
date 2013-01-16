/* 
 * grid.cc
 */
using namespace std;

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <string>
#include <cstring>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include "grid.h"
#include <sys/stat.h>

/* Constructor - destructor */
GridXYZ::GridXYZ(){}
GridXYZ::~GridXYZ(){
delete [] x;
delete [] y; 
delete [] z;
}

void GridXYZ::initialize_grid(){
//cout << "GridXYZ::initialize(). Size: " << size << endl;
ifstream f (fName.c_str(), ios::in);
string line;

if(f.is_open()){
//cout << "File is opened"  << endl;
} else {
cout << "Error. File is NOT opened."  << endl;
}

int i=-1;

while(f.good()){
getline (f,line);
//cout << line << endl;
i++;
}
f.clear();
f.seekg(0,fstream::beg);
size=i;

x = new double[size];
y = new double[size];
z = new double[size];


double a,b,c;
for(i=0;i<size;i++){
f >> x[i] >> y[i] >> z[i];
}
f.close();
}

/* 
 * Linear interpolation routine
 */

double GridXYZ::interpolate_z(double xx, double yy){
//cout << "Grid, x, y: " << xx << " " << yy << endl;
int i=0;
int j=0;

double zz;
double m_x0, m_x1, m_y;
double b_x0, b_x1, b_y;
double dx, dy;
double Dx, Dy;
bool condition;

// At high redshift we interpolate to w_0
//if( yy > gsl_vector_get(vy,size-1) ) {

if( yy < 1. ){ 
//cout << "grid(), y out of range for parameters x,y:  " <<  xx  <<  "  "  << yy << endl; 
yy=1.; }
if( yy > y[size-1] ) {
zz = w0;

// Else make a linear interpolation 
} else {

do{
condition = xx>=x[i] && xx<x[i+repeat]; 
if(!condition)i++;
} while(!condition);

if(i+repeat<size){
x0 = x[i];
x1 = x[i+repeat];

dx = x1 - x0;
Dx = x1 + x0;
} else {
cout << "grid(), x out of range" << endl; 
}
//cout << "dx: " << dx << endl;

j=i -1;
do{
condition = yy>=y[j] && yy<y[j+1]; 
if(!condition)j++;
}while(!condition);

y0 = y[j];
y1 = y[j+1];

dy = y1 - y0;
Dy = y1 + y0;
//cout << "dy: " << dy << endl;
//cout << "y0, y1: " << y0 << " " << y1 << endl;

z00 = z[j];
z01 = z[j+1];
z10 = z[j+repeat];
z11 = z[j+1+repeat];

// Linear interpolation scheme
m_x0 = (1./dx)*(z10 - z00);
b_x0 = (0.5) * ((z10 + z00) - m_x0*Dx);
z_0  = m_x0 * xx + b_x0; 

m_x1 = (1./dx) * (z11 - z01);
b_x1 = (0.5)  * ((z11 + z01) - m_x1*Dx );
z_1  = m_x1  * xx + b_x1;

m_y  = (1./dy) * (z_1 - z_0);
b_y  = (0.5)  * ((z_1 + z_0) - m_y*Dy );
zz   =   m_y * yy + b_y;
}

return zz;
}

#ifdef TEST
// Debugging only
int main(){

GridXYZ grid;

double z;
char* fname="/home/edoardo/wtable.dat";

grid.set_file_name(fname);
//grid.set_eps(0.1);
grid.set_repeat(720);
//grid.set_repeat(3);
grid.set_w0(1./3.);
grid.initialize_grid();

double xx, yy;

cout << "Interpolate value at x, y: ";
cin >> xx >> yy ;
cout << endl;  

z=grid.interpolate_z(xx,yy);

cout << "Interpolated z: " << z << endl;

return 0;
}
#endif
