
#include "ObjectInSpace.h"
#include "global.h"
 
ObjectInSpace::ObjectInSpace(double value,vector<double> coords):
         value(value),coords(coords){
  return;
}
ObjectInSpace::ObjectInSpace(const ObjectInSpace& old):
         value(old.value),coords(old.coords){
  return;
}

double ObjectInSpace::Value() const{
  return value;
}
double ObjectInSpace::Coord(int k) const{
  return coords[k];
}
void ObjectInSpace::setValue(double value){
  this->value=value;
  return;
}
void ObjectInSpace::setCoord(int k,double coord){
  this->coords[k]=coord;
  return;
}
int ObjectInSpace::NDim() const {
  return coords.size();
}
double ObjectInSpace::Distance(ObjectInSpace& compOb,vector<double> scales){
  if(compOb.NDim() != this->NDim()) throw Bad_Error("objects not comparable");
  double sum=0;
  for(int i=0;i<this->NDim();i++) 
	  sum+=pow((compOb.Coord(i)-this->Coord(i))/scales[i],2); 
  return sqrt(sum);
}
void ObjectInSpace::Print(){
  for(int i=0;i<this->NDim();i++) printf("%g ",this->Coord(i));
  printf("%g\n",this->Value());
  return;
}

ostream& operator<<(ostream& os,const ObjectInSpace& p){
  for(int i=0;i<p.NDim();i++) os << p.Coord(i) << " ";
  os << p.Value();
 
  return os;
}


ObjectInSpace::~ObjectInSpace(){}

