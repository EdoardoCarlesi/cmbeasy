
#ifndef OBJECTINSPACE_FLAG
#define OBJECTINSPACE_FLAG

#include "CPP.h"
#include "myutils.h"
#include "mynrutils.h"

class ObjectInSpace {
 
  private:

  double value;
  vector<double> coords;

  public:
   
    ObjectInSpace(double value,vector<double>coords);
    ObjectInSpace(const ObjectInSpace& old);
    double Value() const;
    double Coord(int k) const;
    void setValue(double value);
    void setCoord(int k, double coord);
    int NDim() const;
    double Distance(ObjectInSpace& compOb,vector<double> scales);
    void Print();
    class CompareDistance {
        ObjectInSpace& compOb;
	vector<double> scales;
      public:
        CompareDistance(ObjectInSpace& compOb,
		 const vector<double> scales):compOb(compOb),scales(scales){
	}
	bool operator()(ObjectInSpace O1,ObjectInSpace O2){
          double dis1=O1.Distance(compOb,scales);   
          double dis2=O2.Distance(compOb,scales);  
	  if(dis1<dis2) return true;
	  return false; 
        }
    };

    friend ostream& operator<<(ostream& os, const ObjectInSpace& p);

    ~ObjectInSpace();
};

#endif


