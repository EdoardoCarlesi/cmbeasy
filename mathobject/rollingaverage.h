#ifndef ROLLINGAVERAGE
#define ROLLINGAVERAGE


/*!
  Small class to compute a "rolling" average over values.
  The constructer takes as an argument the maximum number
  of values. By calling push(x), the value x will be stored and
  if the number of values exceeds the limit size, the one stored
  earliest will be erased. Mathematically speaking:

  average = (\sum_0^size x_i ) / min(size,#values stored)

  and at each push(), x_i = x_i+1 
  
  The average is updated at each call to average() by updating
  the sum of all values. Each 100 calls to average(), it is computed
  from scratch to limit elimination of significant figures
  

*/

struct RollingAverage {
  unsigned int Size; //!< the maximum number of values
  unsigned int Number;  //!< the number of values
  double Sum; //!< the sum of all values
  double *y; //!< the values
  unsigned int idx;  //!< current position for push()
  int PerformFullCounter;  //!< Every 100 calls to average(), we re-compute the average from scratch
  RollingAverage(unsigned int size);
  ~RollingAverage();
  
  void push(double x); //!< push x into the average-array
  void clear();  //!< clear the array
  double average(); //!< compute the average
};


#endif
