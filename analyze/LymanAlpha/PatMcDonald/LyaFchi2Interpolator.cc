

#include "LyaFchi2Interpolator.h"
#include "mynrutils.h"
#include "controlpanel.h"

/*
extern "C"
{
  void coordssave_ (double *coords);
}
*/

const double LyaFchi2Interpolator::lkcent=log(0.009);
const double LyaFchi2Interpolator::lkstep=log(1.1);
 
LyaFchi2Interpolator::LyaFchi2Interpolator(string tf,double s0,double s1, 
                             double s2): tablefile(tf) {
  
  double chisq;
  vector<double> coords(3);
  vector<ObjectInSpace> voP;
  
  if (tablefile == "default") tablefile = ControlPanel::cmbeasyDir("/resources/LymanAlpha/PatMcDonald/lyafchisq_pcagrid.txt");

  //ifstream fp(tablefile.c_str());
  ifstream fp(tablefile.c_str());
  if(fp==NULL) {
    cout << tablefile << endl;
    throw Bad_Error("No chi^2 table file.");
  }
  string buf;
  while(getline(fp,buf)){
    istringstream(buf) >> coords[0] >> coords[1] >> coords[2] >> chisq;
    if(chisq>0) voP.push_back(ObjectInSpace(chisq,coords));
  }
  fp.close();
  coords[0]=s0;
  coords[1]=s1;
  coords[2]=s2;
/*
  coords[0]=(0.11995-0.096197)*100;
  coords[1]=(2.14193-2.11951)*100;
  coords[2]=0.05*100;
*/
  gi = new  GlassInterpolator(voP,coords,1.0e-6);
  return;
}    

/*!
  New chi2 routine for Lyman Alpha Data. Differs in the arguments it
  takes from the original one chi2_old(). 
  You have to supply directly the Delta^2, neff = dlnP/dlnk and running of the 
  effective spectral index d^2lnP/dlnk^2 at the pivot scale k =0.009 s/km.
  See AnalyzeThis::lymanAlphaPatMcDonaldChi2() for an example
*/
double LyaFchi2Interpolator::chi2(double Delta2, double neff, double dndlnk) {
  vector<double> coords(3);
  coords[0] = Delta2; 
  coords[1] = neff;
  coords[2] = dndlnk;

  ObjectInSpace P(0,coords);
  gi->Interpolate(&P);
  return P.Value();
}


//expects arrays to start at 0
double LyaFchi2Interpolator::chi2_old(int N, double *k, double *D2){
  vector<double> coords(3);
  //really need to objectify this splining procedure
  double *lk,*lD,*lDsp;
  lk=NR::dvector(1,N);
  lD=NR::dvector(1,N);
  lDsp=NR::dvector(1,N);
  for(int i=0;i<N;i++){
    lk[i+1]=log(k[i]);
    lD[i+1]=log(D2[i]);
  }
  NR::oldspline(lk,lD,N,1.0e31,1.0e31,lDsp);
 
  double Dcent;
  NR::oldsplint(lk,lD,lDsp,N,lkcent,&Dcent);
  double Dup;
  NR::oldsplint(lk,lD,lDsp,N,lkcent+lkstep,&Dup);
  double Ddown;
  NR::oldsplint(lk,lD,lDsp,N,lkcent-lkstep,&Ddown);
  coords[0]=(exp(Dcent)+exp(Ddown)+exp(Dup))/3.0;
  coords[1]=(Dup-Ddown)/(2.0*lkstep)-3;
  coords[2]=(Dup-2.0*Dcent+Ddown)/(lkstep*lkstep);
  
  cout << "Pat's coords look up: coords[0]: " << coords[0] << "  coord[1]:  "<< coords[1] << "  coords[2]: " << coords[2] << endl;

  // cout << "in coords:" << coords[0] << coords[1] << coords[2] << endl;
  //coordssave_(&coords[0]);

  ObjectInSpace P(0,coords);
  //cout << P ;
  gi->Interpolate(&P);
  double val=P.Value();
  NR::free_dvector(lk,1,N);
  NR::free_dvector(lD,1,N);
  NR::free_dvector(lDsp,1,N);

  return val;
}

LyaFchi2Interpolator::~LyaFchi2Interpolator(){
  delete gi;
}

