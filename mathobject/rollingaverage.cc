#include "rollingaverage.h"
#include "global.h"


RollingAverage::RollingAverage(unsigned int size) : Size(size), Number(0), Sum(0), y(0), idx(0), PerformFullCounter(0) {
  y = new double[size];
  if (!y) throw Bad_Error("RollingAverage:: couldn't allocate memory");
  for (unsigned int i =0; i < size; i++) y[i] = 0;
}

RollingAverage::~RollingAverage() {
  delete[] y;
}


void RollingAverage::clear() {
  for (unsigned int i =0; i < Size; i++) y[i] = 0;
  Number = 0;
  Sum = 0;
  idx = 0;
  PerformFullCounter = 0;
}

void RollingAverage::push(double x) {
  Sum -= y[idx]; // subtract the value we will overwrite
  Sum += x;  // add the new value	       
  y[idx++] = x; // store the new value
  if (idx == Size) idx = 0; // roll-over
  Number++; // one more value
  if (Number >= Size) Number = Size; // at most we have Size numbers in the array
}

double RollingAverage::average() {
  if (!Number) throw Bad_Error("RollingAverage::average() zero number");
  if (PerformFullCounter++ < 10) return Sum/Number;
  PerformFullCounter = 0; // reset 
  Sum = 0;
  for (unsigned int i = 0; i < Number; i++) {  // at most to number, which might later be size
    Sum += y[i];
  }
  return Sum/Number;
}

