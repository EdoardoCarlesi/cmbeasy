
#include "myutils.h"

namespace MU {

void myexit(char error_text[])
/* based on Numerical Recipes standard error handler */
{
        fprintf(stderr,"My run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}

void reversevectororder(double *v,int N){
  double *t;
  t=NR::dvector(1,N);
  for(int i=1;i<=N;i++) t[i]=v[i];
  for(int i=1;i<=N;i++) v[i]=t[N-i+1];
  NR::free_dvector(t,1,N);
  return;
}

}
