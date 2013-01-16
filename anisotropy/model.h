#ifndef MODEL_H
#define MODEL_H

#define MODEL_VERSION 3

class Spline;
class CL;
class Cosmos;
class Anchor;
class ControlPanel;

using namespace std;
#include <map>
#include <iostream>
#include <string>

//! Exception thrown if version of model stored in file is not known to Model
struct InvalidModelVersion {};

struct ModelData {
  ModelData() : total(0), s(0) {};
  virtual ~ModelData() { if (s) delete[] s; }
  float min_x, max_x;
  float min_y, max_y;
  unsigned short total;
  unsigned short *s;
  void readData(istream& i);
  void writeData(ostream& o);
  virtual void data( map<float,float> & ) {};
  virtual void compress(Spline*,unsigned short,double=1e100)=0;
};

struct LinearModelData : public ModelData {
  void data( map<float,float> & );
  void compress(Spline*,unsigned short,double=1e100 );
};

struct LogModelData : public ModelData {
  void data( map<float,float> & );
  void compress(Spline*,unsigned short,double=1e100);
};

struct ClModelData : public ModelData {
  ClModelData();
  static double X[100];
  static bool initializedX; 
  void data( map<float,float> & );
  void compress(Spline*,unsigned short,double=1e100 );
};



class Model {

  //Spline* ts;
  //Spline* powerCdm;
  
  
  ModelData  *ts,*tt;
  ModelData *luminosity;
 
  LogModelData  powerCdm;
  
  ClModelData es,cs; // E and cross-corrolation


 public:
  enum WhichData {linear, clts};
  float sigma8;

  float o_b;
  float o_cdm;
  float o_lambda;
  float o_qls;
 
  float h;
  float n;
  float optdlss;
 
  float n_t;  // tensor
  float o_quintessence;
  float mp[5];

 public:
  Model(Cosmos&, CL&,const ControlPanel&,int pn,double stopcl=1e100);
  Model(istream&);
  Model(Spline*,Cosmos&,Spline*, WhichData=linear);
  ~Model();  

  void write(ostream&);
  //void writeSpline(ostream&,Spline*);
  //void readSpline(istream&,Spline**);
  void outputScalarCl(string filename);

  void printStatus();

  //template<class T> static T  read(istream&);
  //template<class T> static void write(ostream&,const T&);
  template<class T> static T read(istream& i) {
    T s;
    i.read((char*) &s,sizeof(T));
    return s;
  }

  template<class T> static void write(ostream& o,const T& s) {
    o.write((const char*) &s,sizeof(T));
  }

  Spline *spline(Anchor*, ModelData*, string);
  Spline* cl_ts(Anchor*);
  Spline* cl_es(Anchor*);
  Spline* cl_cs(Anchor*);
  Spline* cl_tt(Anchor*);
  
  Spline* power_cdm(Anchor*); 
  Spline * d_l(Anchor*); //!< luminosity distance (little d_l)
};



#endif 
 
