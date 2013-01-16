#include "data.h"
#include "newdatreader.h"
#include "global.h"

#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <map>
#include <iostream>
#include <algorithm>

Data::Data(const char* fileName, const char* name, double normal) : Name(name) , FileName(fileName) , Normal(normal) {
  read(XYerror);
}

Data::Data(const char* fileName, const char* name, Type typ, double normal) : Name(name) , FileName(fileName) , Normal(normal) {
  read(typ);
}

Data::Data(const char *name) : Name(name), FileName(0), Normal(1.0) {}


void Data::read(Type typ) { 
  if (FileName == 0)
    throw Bad_Error("Data::read() You used constructor without filename. No read() possible");
  if (typ==Newdat) {
    createFromNewdat();
    return;
  }


  ifstream file(FileName);

  if (!file) {
    cout << "FILE NOT FOUND: " << FileName << endl;
    throw Bad_Error("Data::Data() Data file not found");
  }
  char buffer[1000];
  double X,Y,dx=0,dy;
  double xMin, xMax;
  //ofstream o1;
  //ifstream i1;
  //  cout << "FILE: " << FileName << endl;
  while (file.getline(buffer,999)) {
    int pos = 0,count=0;
    double value[100];
    bool nocomment=true;
    while (buffer[pos] != 0) {
      if (buffer[pos++] == '#') nocomment=false; 
    }
    pos = 0;
    if (nocomment) {
      do {
        while (isspace(buffer[pos])) pos++; // eat beginning whitespaces
        string s;
        while (! isspace(buffer[pos]) && buffer[pos] != 0) s+= buffer[pos++]; // not 0 and not white space
        //	cout << "converting: ["<<count<<"] : " << s << "  || ";
        value[count++] = strtod(s.c_str(),'\0');
      } while (buffer[pos]);
      //      cout << endl;
      DataEntry tmpEntry;
      switch (typ) {
        case Yerror:
          X = value[0];
          Y = value[1];
          dy = value[2];
          points.push_back(DataEntry(X,Y*Normal,dx,dy*Normal));
          break;
        case XYerror:
          X = value[0];
          Y = value[1];
          dx = value[2];
          dy = value[3];
          points.push_back(DataEntry(X,Y*Normal,dx,dy*Normal));
          break;
        case AsymXYerror:
          X = value[0];
          xMin = value[1];
          xMax = value[2];
          Y = value[3];
          dy = value[4];
          tmpEntry = DataEntry(X,Y*Normal);
          tmpEntry.setYError(dy*Normal);
          tmpEntry.setXWithErrors(X, xMin, xMax);
          points.push_back(tmpEntry);
          break;
        default:
          throw Bad_Error("Data::read() unknown error type");
      }
    }
  }
}

void Data::createFromNewdat()
{
  NewdatReader reader(FileName);
  typedef std::vector<NewdatReader::Band> Bands;
  const Bands& ttBands = reader.TTBands();
  Bands::const_iterator it, end;
  end = ttBands.end();
  for ( it = ttBands.begin(); it != end; ++it) {
    NewdatReader::Band band = *it;
    DataEntry entry;
    entry.setXWithErrors(band.lEff, band.lMin, band.lMax);
    entry.setYWithErrors(band.power, band.power-band.minusError, band.power+band.plusError);
    points.push_back(entry);
  }
}



vector< vector< DataEntry> >* Data::fill2d() {
  unsigned int n = (unsigned int)sqrt((double)points.size());
  if (n*n != points.size()) return 0;
  
  
  vector< vector< DataEntry> > * v = new vector< vector<DataEntry> >(n);
  for (unsigned int i=0; i < n; i++) (*v)[i].resize(n);
  
  list<DataEntry>::iterator k = points.begin();
  for (unsigned int i=0; i <n; i++) {
    for (unsigned int j=0; j<n;j++) {
      (*v)[i][j] = *k;
      k++;
    }
  }
  return v;
} 

void Data::sort() {
  points.sort();
}
