
#include "GlassInterpolator.h"
 
GlassInterpolator::GlassInterpolator(vector<ObjectInSpace> voP,
              const vector<double> scales,double SVDtol):
                      voP(voP),scales(scales),SVDtol(SVDtol),
           		NDim(voP.begin()->NDim()),Nd((2+NDim)*(1+NDim)/2+1){
  return;
}    

void GlassInterpolator::Interpolate(ObjectInSpace* P){
  int ma=(2+NDim)*(1+NDim)/2;
  int nsvduse=0;
  do { 
  vector<bool> upbound(NDim,false);
  vector<bool> downbound(NDim,false);
  partial_sort(voP.begin(),voP.begin()+Nd,voP.end(),
		   ObjectInSpace::CompareDistance(*P,scales));
  vector<ObjectInSpace> voPuse;
  for(int i=0;i<Nd;i++) {
    for(int j=0;j<NDim;j++) {
    //  printf("%d %g %g ",j,voP[i].Coord(j),P->Coord(j));
      if(voP[i].Coord(j)>=P->Coord(j)) upbound[j]=true; 
      if(voP[i].Coord(j)<=P->Coord(j)) downbound[j]=true; 
    }
     // printf("%d\n",i);
    voPuse.push_back(voP[i]);
  }
  for(int j=0;j<NDim;j++) {
    //should return this information instead of complaining
    if(!upbound[j]) { 
      printf("coord %d breaks contain above\n",j); 
      throw GlassInterpolatorOutOfBound();
    }
    if(!downbound[j]) {
      throw GlassInterpolatorOutOfBound();
      printf("coord %d breaks contain below\n",j); 
    }
  }
  //stuff below should be an object with first part setup and second part 
  //callable then this function would just have to decide if the training 
  //list had changed and a new setup was needed
  double* y = NR::dvector(1,Nd);
  double* sig = NR::dvector(1,Nd);
  double* a = NR::dvector(1,ma);
  double** T = NR::dmatrix(1,ma,1,Nd);
  double maxsep=P->Distance(*(voPuse.end()-1),scales)*(1.0+mytinyfloat);
  for(int i=0;i<Nd;i++){
    y[i+1]=voPuse[i].Value();
    double sep=P->Distance(voPuse[i],scales);
    if(DEQUAL(sep,0)) { 
      P->setValue(voPuse[i].Value()); 
      NR::free_dvector(a,1,ma);
      NR::free_dvector(y,1,Nd);
      NR::free_dvector(sig,1,Nd);
      NR::free_dmatrix(T,1,ma,1,Nd);
      return;
    }
    sig[i+1]=(maxsep*sep)/(maxsep-sep);
    int j=1; 
    T[j][i+1]=1;
    for(int k=0;k<NDim;k++) {
      double coord1=voPuse[i].Coord(k)/scales[k];
      T[(++j)][i+1]=coord1;
      for(int l=k;l<NDim;l++) 
        T[(++j)][i+1]= voPuse[i].Coord(l)*coord1/scales[l];
    }
  }

  double chisq=0, maxoff=0;
  NR::myweightedsvdlsq(T,y,sig,Nd,a,ma,&chisq,&maxoff,&nsvduse,SVDtol,false); 
 
  double result=0.0;
  int j=1;
  result+=a[j];
  for(int k=0;k<NDim;k++) {
    double coord1=P->Coord(k)/scales[k];
    result+=a[(++j)]*coord1;
    for(int l=k;l<NDim;l++) 
      result+=a[(++j)]*P->Coord(l)*coord1/scales[l];
  }
  P->setValue(result);
  NR::free_dvector(a,1,ma);
  NR::free_dvector(y,1,Nd);
  NR::free_dvector(sig,1,Nd);
  NR::free_dmatrix(T,1,ma,1,Nd);
  if(nsvduse<ma) Nd++;
  } while(nsvduse<ma);

  return; 
}

GlassInterpolator::~GlassInterpolator(){}

