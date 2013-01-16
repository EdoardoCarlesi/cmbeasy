#ifndef DATA_H
#define DATA_H

#include "global.h"

#include <list>
#include <vector>
using namespace std;

/*!
  Small class holding one single data point with x,y,
  and (possibly asymetric) x- and y-errors
*/
class DataEntry
{
  public:
    //! constructor for symmetric x (and y) errors
    DataEntry(const float X=0, const float Y=0, const float DX=0.0, const float DY=0.0)
      : mX(X), mY(Y) {
        mXmErr = X-DX;
        mXpErr = X+DX;
        mYmErr = Y-DY;
        mYpErr = Y+DY;
        mXErrIsSymmetric=mYErrIsSymmetric=true;
      }

    void setXWithErrors(const float X, const float xMin, const float xMax) {
      mXErrIsSymmetric = false;
      mX = X;
      mXmErr = xMin;
      mXpErr = xMax;
    }

    void setYWithErrors(const float Y, const float yMin, const float yMax) {
      mYErrIsSymmetric = false;
      mY = Y;
      mYmErr = yMin;
      mYpErr = yMax;
    }

    void setYError(const float dy) {
        mYErrIsSymmetric = true;
        mYmErr = mY-dy;
        mYpErr = mY+dy;
    }

    float x() const { return mX; }
    float xMiddle() const { return x(); }
    float xLeft() const { return mXmErr; }
    float xRight() const { return mXpErr; }
    float dx() const {
      if (!mXErrIsSymmetric)
        throw Bad_Error("data.h: DataEntry::dx() called for point with asymetric errors.");
      return mX-mXmErr;
    }

    float y() const { return mY; }
    float yMiddle() const { return y(); }
    float yLower() const { return mYmErr; }
    float yUpper() const { return mYpErr; }
    float dy() const {
      if (!mYErrIsSymmetric)
        throw Bad_Error("data.h: DataEntry::dy() called for point with asymetric errors.");
      return mY-mYmErr;
    }

    void operator*=(const float a) {
      mY *= a; mYmErr *= a; mYpErr *= a; }

    bool operator<(const DataEntry& e) const {
      return mX < e.x();
    }

  private:
    float mX, mY;
    float mXmErr;    //!< x value minus the xError
    float mXpErr;  //!< x value plus the xError
    float mYmErr;    //!< x value minus the xError
    float mYpErr;  //!< x value plus the xError

    bool mXErrIsSymmetric;
    bool mYErrIsSymmetric;
};

/*!
  Data provides basic read capabilities of ASCII data files.
  The data fromat needs to have one of the formats
  X Y DX DY   (XYerror)
  X Y DY   (Yerror)
  X Xmin Xmax Y DY   (AsymXYerror)
  and "#" marks comment lines.
  All data points are stored in the list<DataEntry> points, which you can
  publically access.
*/
class Data {
 public:
  enum Type { XYerror, Yerror, AsymXYerror, Newdat};
  const char* Name;
  list<DataEntry> points; //!< The data as DataEntry's
  const char *FileName; //!< The name of the file to read from
  double Normal;  //!< normalization, i.e. factor with which Y and dY will be multiplied on reading from file

  Data(const char* fileName, const char*name,double normal=1.0); //!< Create Data object and read in data from file 
  Data(const char* fileName, const char*name, Type t, double normal=1.0); //!< Create Data object and read in data from file
  Data(const char* name); //!< Rarely used. 
  void read(Type t); //!< read in the data from file 
  void createFromNewdat(); //!< use a NewdatReader object to create data from a file in .newdat format (e.g. Boomerang 03, Acbar '08 etc)

  const Data& operator *=(const double factor) {
    for (list<DataEntry>::iterator i = points.begin(); i != points.end(); i++) {
      (*i) *= factor;
    }
    return *this;
  }

  vector< vector< DataEntry> >* fill2d();
  void sort(); // sort the data points 

};

#endif
