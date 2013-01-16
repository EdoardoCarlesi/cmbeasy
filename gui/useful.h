#ifndef USEFUL
#define USEFUL

#include <map>

/*!
  Return the element of the best matching
  index of this map
*/
template<class S> S* bestMatch2(map<double, S> &mp, double k, double delta) {  
  typedef typename map<double,S>::iterator bmIterator;
  if (mp.empty()) return 0;
  bmIterator i =  mp.lower_bound(k);  // either exact match or i is at the next higher value
  S* best = 0;
  double D = 2*delta + 1e99;  // at least bigger then delta :-)
  if (i != mp.end()) { best = &(i->second);   D = fabs(i->first - k); } 

  
  if (i != mp.begin()) {
    double D2 = fabs( (--i)->first - k);
    if (D2 < D)  { D = D2; best = &(i->second); }
  }

  if (D <= delta) return best; else return 0; 
}

/*!
  Return the best matching double value, i.e. the
  index of this map 
*/

template<class S> double bestMatch1(map<double, S> &mp, double k, double delta) {  
    typedef typename map<double,S>::iterator bmIterator;
  if (mp.empty()) return 0;
  bmIterator i =  mp.lower_bound(k);  // either exact match or i is at the next higher value
  double best = 0;
  double D = 2*delta + 1e99;  // at least bigger then delta :-)
  if (i != mp.end()) { best = (i->first);   D = fabs(i->first - k); } 
  
  if (i != mp.begin()) {
    double D2 = fabs( (--i)->first - k);
   
    if (D2 < D)  { D = D2; best = (i->first); }
  }

  if (D <= delta) return best; else return 0; 
}



struct SequentialOutOfRange {};  // error struct for function below
/*!
  transverse through map, counting from 0 to q 
  and return the index (i.e. typ S) corresponding
  to the element reached in this map
*/ 
template<class S,class T> S  sequential(const map<S,T> &m, unsigned int q) {
  typedef typename map<S,T>::const_iterator SequentialIterator;
  unsigned int c = 0;
  if (q >= m.size())  throw SequentialOutOfRange();
  SequentialIterator i = m.begin();
  while (c < q) { i++; c++;}
  return i->first;
}

#endif
