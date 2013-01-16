#include "model.h"
#include "cosmos.h"
#include <iostream>
#include "spline.h"
#include "cl.h"
#include "anchor.h"
#include "controlpanel.h"

Model::Model(Cosmos& cosmos, CL& cl,const ControlPanel& control,int pn,double stopcl) : o_qls(0)  {

  ts = new ClModelData();
  tt  = new ClModelData();

  luminosity = new LinearModelData();
  Spline lum(100, "Model::lum");
  for (double z = 0; z <= 2; z+= 0.05) 
    lum.set(z, cosmos.luminosityDistance(z)); 
  luminosity->compress(&lum, 30);

  ts->compress(cl.ts[pn],60,stopcl); 
  tt->compress(cl.tt[pn],60,stopcl);
  es.compress(cl.es[pn],60,stopcl);
  cs.compress(cl.cs[pn],60,stopcl);
  
  if (control.power_cdm) {
    Spline *pc = cosmos.createPower(pn,"cdm",cosmos.power_cdm());
    //Spline *pc = cosmos.power_cdm(pn)->createAlongX(cosmos.tau_0(),"cdm");
    powerCdm.compress(pc,120); delete pc;
  }

  sigma8 = cosmos.sigma8[pn]; 
  o_b = cosmos.omega_b();
 
  if (cosmos.validHistory()) 
  o_qls = cosmos.tau2AvOmega_q(cosmos.tau_ls());

  o_cdm = cosmos.omega_cdm();
  o_lambda = cosmos.omega_v();
  h = cosmos.h();
  n = cosmos.InitialPower[pn];
  optdlss = cosmos.optDistanceLss();
}



Model::Model(istream& i) :  ts(0), tt(0), luminosity(0), o_qls(0) {
  short version = read<short>(i);

  if (version < 2 || version > 3 || true)  throw InvalidModelVersion();
  

  if (version ==2) {
    o_b = read<float>(i);
    o_cdm = read<float>(i);
    o_lambda = read<float>(i);
    h = read<float>(i);
    sigma8 = read<float>(i);
    n =  read<float>(i);
    optdlss =  read<float>(i);   
    ts = new LinearModelData();
    ts->readData(i);
    powerCdm.readData(i); 

    //printStatus();

  }

  if (version == 3) {
    o_b = read<float>(i);
    o_cdm = read<float>(i);
    o_lambda = read<float>(i);
    h = read<float>(i);
    sigma8 = read<float>(i);
    n =  read<float>(i);
    optdlss =  read<float>(i);   

    n_t = read<float>(i);   
    o_quintessence = read<float>(i);   
    for (int k =0; k < 5; k++) mp[k] = read<float>(i);  
    
    ts = new ClModelData();
    tt  = new ClModelData();
    luminosity = new LinearModelData();

    ts->readData(i);
    tt->readData(i);
    es.readData(i);
    cs.readData(i);
    powerCdm.readData(i);   
    luminosity->readData(i);
  }
}

/*!
  Construct model from cosmos and splines. If the pointers to 
  the splines are 0, compress() will not process it and thus
  you may use these models to store ob etc with just
  a little overhead 
*/
Model::Model(Spline* cl, Cosmos& cosmos, Spline* cdm, WhichData which) : luminosity(0) , o_qls(0) {
  if (which == clts) ts = new ClModelData();
  if (which == linear) ts = new LinearModelData();
  ts->compress(cl,60);  
  powerCdm.compress(cdm,40);

  sigma8 = 0;
  
  if (cosmos.validHistory()) 
  o_qls = cosmos.tau2AvOmega_q(cosmos.tau_ls()); 
  
  o_b = cosmos.omega_b();
  o_cdm = cosmos.omega_cdm();
  o_lambda = cosmos.omega_v();
  h = cosmos.h();
  n = cosmos.InitialPower[0];
  optdlss = cosmos.optDistanceLss();
}

Model::~Model() {
  //cout << "ending model" << endl;
  if (ts) delete ts;
  //cout << "em1" << endl;
  if (tt) delete tt;
  //cout << "em2" << endl;
  if (luminosity) delete luminosity;
  //cout << "em3" << endl;
}


/*
Model::Model(Spline* clts, float ob, float ocdm, float olambda, float H,float s8 = 0,Spline* cdm=0) :   sigma8(s8) {
  
  if (clts) ts = new Spline(*clts);
  if (cdm) powerCdm = new Spline(*cdm);
  o_b = ob;
  o_cdm =ocdm;
  o_lambda = olambda;
  h = H;
}

Model::Model(Spline* clts,Cosmos& cosmos,float s8 = 0,Spline* cdm=0) : ts(0),powerCdm(0), sigma8(s8) {
  
  ts = new Spline(*clts);
  if (cdm)  powerCdm = new Spline(*cdm);
  o_b = cosmos.omega_b();
  o_cdm = cosmos.omega_cdm();
  o_lambda = cosmos.omega_v();
  h = cosmos.h();
}
  */


void Model::printStatus() {
  cout << "Model::printStatus()" << endl;
  cout << "o_b: " << o_b << endl;
  cout << "o_cdm: " << o_cdm << endl;
  cout << "h: " << h << endl;
  cout << "n: " << n << endl;
  cout << "optdlss: " << optdlss << endl;
}

/*!
  Output scalar cl's in CMBFAST - format
*/
void Model::outputScalarCl(string filename) {
  ofstream scalarFile(filename.c_str());
  Anchor anchor;
  Spline *t = cl_ts(&anchor);
  Spline *e = cl_es(&anchor);
  Spline *c = cl_cs(&anchor);
  t->arm();e->arm();c->arm();

  for (int l = 2; l <= t->stop(); ++l) {
      scalarFile << l << " " << (*t)(l);  
      scalarFile << " " << (*e)(l) << " " << (*c)(l) << endl;
  }


}

void Model::write(ostream& o) {
  //printStatus();
  
  write<short>(o,MODEL_VERSION);
  write<float>(o,o_b);
  write<float>(o,o_cdm);
  write<float>(o,o_lambda);
  write<float>(o,h);
  write<float>(o,sigma8);
  write<float>(o,n);
  write<float>(o,optdlss);
  
  write<float>(o,n_t); 
  write<float>(o,o_quintessence);
  
  for (int k =0; k < 5; k++)  write<float>(o,mp[k]);
  
  ts->writeData(o);
  tt->writeData(o);
  es.writeData(o);
  cs.writeData(o);
  powerCdm.writeData(o);
  luminosity->writeData(o);
}

/*
void Model::writeSpline(ostream& o, Spline* s,bool compress) {
  if (compress)  {
      } else { 
    write<short>(o,s->size());
    for (int i = 1; i <= s->size(); i++) { 
      write<float>(o,s->x(i));
      write<float>(o,s->y(i));
    }
  }
}
void Model::readSpline(istream& i, Spline** s) {
  if (*s !=0) delete s;
  *s = new Spline(1000,"ModelReadSpline");
  int n = read<short>(i);
  for (int j = 1; j <=n; j++) { 
    (*s)->setForce(read<float>(i),read<float>(i));
  }
}
*/

/*!
  Create Spline that contains the cl.ts data.
  The Spline necessarily should have an anchor
  in order to prevent memory holes, so please
  do not set a = 0
*/
Spline* Model::cl_ts(Anchor* a) {
  return spline(a, ts, "model_ts");
}
Spline* Model::cl_es(Anchor* a) {
  return spline(a, &es, "model_es");
}
Spline* Model::cl_cs(Anchor* a) {
  return spline(a, &cs, "model_cs");
}
 /*!
  Create Spline that contains the cl.tt data.
  The Spline necessarily should have an anchor
  in order to prevent memory holes, so please
  do not set a = 0
*/
Spline* Model::cl_tt(Anchor* a) {
  return spline(a, tt, "model_tt");
}

Spline* Model::d_l(Anchor* a) {
  return spline(a, luminosity, "model_luminosity");
}


Spline* Model::spline(Anchor *a, ModelData* daten,string name) {
  map<float,float> m;
  Spline *s = new Spline(daten->total, name,a);
  daten->data(m); 
  s->set(m);
  return s;
}
 

Spline* Model::power_cdm(Anchor *a) {
  map<float,float> m;
  //cout << "power_cdM() " << endl;
  Spline *s = new Spline(powerCdm.total, "model_cdm",a);
  //cout << "trying to get data" << endl;
  powerCdm.data(m);
  //cout << "have data" << endl;
  s->set(m);
  //s->printStatus();
  return s;
}

void ModelData::readData(istream &i) {
  min_x = Model::read<float>(i);
  max_x = Model::read<float>(i);
  min_y = Model::read<float>(i);
  max_y = Model::read<float>(i);
  //cout << "readData: " << min_x << "  " << max_x << "  " << min_y << "  " << max_y << endl;
  total = Model::read<unsigned short>(i);  
  // cout << "readData: total: " << total << endl;
  if (s) delete[] s;
  s = new unsigned short[total];
  for (short q = 0; q < total; q++) s[q] = Model::read<unsigned short>(i);
}

void ModelData::writeData(ostream &o) {
  Model::write<float>(o,min_x);
  Model::write<float>(o,max_x);
  Model::write<float>(o,min_y);
  Model::write<float>(o,max_y);
  Model::write<unsigned short>(o,total);
  for (short q = 0; q < total; q++)  Model::write<unsigned short>(o,s[q]);
}

void LinearModelData::data( map<float,float>&m ) {
  m.clear();
  if (total > 1 && s) {
    float xstep = (max_x - min_x)/(total-1);
    float ystep = (max_y - min_y)/(pow(256.0,(double)sizeof(unsigned short))-1);
    for (unsigned short q = 0; q < total; q++) {
      float x = min_x + q*xstep;
      m[x] = min_y + s[q] * ystep;
    }
  }
}

void LogModelData::data(map<float,float>&m ) {
  m.clear();
  // cout << "logmodel " <<  total << endl;
  if (total > 1 && s) {
    //cout << "logmodelData: in " << endl;
    float logXStep = (log(max_x) - log(min_x)) / (total-1);
    float logYStep  = (log(max_y) - log(min_y)) /(pow(256.0,(double)sizeof(unsigned short))-1);
    for (unsigned short q = 0; q < total; q++) {
      float x = exp( log(min_x) + logXStep*q );
      m[x] = exp( log(min_y) + logYStep*s[q]);
    }
  }
}

void ClModelData::data(map<float,float>&m ) {
  m.clear();
  if (total > 1 && s) { 
    float ystep = (max_y - min_y)/(pow(256.0,(double)sizeof(unsigned short))-1);
    for (unsigned short q = 0; q < total; q++) {
      m[X[q]] = min_y + s[q] * ystep;  
    }
  }
}

bool ClModelData::initializedX(false);
double ClModelData::X[100];

ClModelData::ClModelData() : ModelData() {
  if (initializedX) return;
  initializedX = true;
  double x=2;
  double step =1;
  for (int i = 0; i < 100; i++) {
    X[i] = x;
    if (x >= 10) step = 2;
    if (x >= 20) step = 10;
    if (x >= 90) step = 20;
    if (x  >= 150) step = 50;
    x += step;
  }
}





void LinearModelData::compress(Spline* spline, unsigned short t, double mXx) {
  if (spline) {
    if (! spline->isArmed()) spline->arm();
    total = t;
    if (s) delete[] s;
    s = new unsigned short[total];
    min_x = spline->start();
    max_x= min(spline->stop(),mXx);
    min_y = spline->front();
    max_y = spline->front();
    float xstep = (max_x - min_x)/(total-1);


    for (unsigned short q = 0; q < total; q++) {
      float x = min_x + q*xstep;
      float y = spline->fastY(x);
      if (y < min_y) min_y = y;
      if (y > max_y) max_y = y;
    }
 
    float ystep = (max_y - min_y)/(pow(256.0,(double) sizeof(unsigned short))-1);
  
    for (unsigned short q = 0; q < total; q++) {
      float x = min_x + q*xstep;
      float y = spline->fastY(x);
      s[q] =  (int) rint( (y - min_y) / ystep);
    }
  }
}

void LogModelData::compress(Spline* spline, unsigned short t,double mXx) {
  if (spline) {
    if (! spline->isArmed()) spline->arm();
    total = t;
    if (s) delete[] s;
    s = new unsigned short[total];
    min_x = spline->start();
    max_x= min(mXx,spline->stop());
    min_y = spline->front();
    max_y = spline->front();
    float xstep = (log(max_x) - log(min_x))/(total-1);
 
    for (unsigned short q = 0; q < total; q++) {
      float x = exp(log(min_x) + q*xstep);
      float y = spline->fastY(x);
      if (y < min_y) min_y = y;
      if (y > max_y) max_y = y;
    }
 
    float ystep = (log(max_y) - log(min_y))/(pow(256.0,(double)sizeof(unsigned short))-1);
  
    for (unsigned short q = 0; q < total; q++) {
      float x =exp( log(min_x) + q*xstep);
      float y = spline->fastY(x);
      s[q] = (int) rint( (log(y) - log(min_y)) / ystep );
    }
  }
}

void ClModelData::compress(Spline* spline, unsigned short t,double mXx) {
  total = 0;

  if (spline) {
    if (! spline->isArmed()) spline->arm();
    min_y = spline->front();
    max_y = spline->front();

    for (int i =0; i < 100; i++) {
      if (X[i] > mXx) break;
      total = i+1;
      float y = spline->fastY(X[i]);
      if (y < min_y) min_y = y;
      if (y > max_y) max_y = y;
    }
    s = new unsigned short[total];
    float ystep = (max_y - min_y)/(pow(256.0,(double) sizeof(unsigned short))-1);
  
    for (short i  = 0; i < total; i++) {
      float y = spline->fastY(X[i]);
      s[i] =  (int) rint( (y - min_y) / ystep);
    }
  }
}
