#include "integrator.h"
#include "cosmos.h"

Integrator::Integrator(Cosmos *c, double mi, double ma)  : cosmos(c), akmin(mi), akmax(ma) {


  if (akmin == 0 && akmax ==0) return;  // if both are zero, k-stepping is different and we return

  //cout << "akmin: " << akmin << "  akmax: " << akmax << endl;

  double tau0 = cosmos->tau_0();
  double dk=2.5/tau0;
  double dk0=1.5/tau0;

  int no=500;
  double  dlnk1=0.1;

  int  no1=int(log(10.0*dk0/akmin)/dlnk1)+1;

  if (akmax > (no*dk0)) {
    nko=int((akmax-no*dk0)/dk)+no;
  }  else {
    no=int((akmax-10.0*dk0)/dk0)+no1;
    nko=no;
  }
  K.resize(nko);
  dK.resize(nko);
  
  for (int k = 1; k <= nko; k++) {
    if (k <= no) {
      if (k <= no1) {
	K[k-1] = 10.0*dk0*exp(-(no1-k)*dlnk1);
	dK[k-1]=K[k-1]*dlnk1;
      } else {
	K[k-1]=K[no1-1]+(k-no1)*dk0;
	dK[k-1]=dk0;
      }
    } else {
      K[k-1]=K[no-1]+(k-no)*dk;
      dK[k-1]=dk;
    }  
  }
  /*
  ofstream kspace("kspace.dat");
  for (int i = 1; i < nko; i++) {
    kspace << K[i] << " " << 1 << endl;
    cout << "i: " << i << "  k: " << K[i] << endl;
  }
  */

  dK[0]=0.5*dK[0];
  dK[no1-1]=0.5*(dK[no1-1]+dk0);
  dK[nko-1]=0.5*dK[nko-1];
}

Integrator::~Integrator() {}
 
