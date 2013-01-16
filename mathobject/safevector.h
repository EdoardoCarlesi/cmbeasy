#ifndef SAFEVECTOR
#define SAFEVECTOR

using namespace std;
#include <vector>
#include <iostream>
#include <list>

struct SafeVectorOutOfRange {};

/*!
  Safe vector class that does range checking for
  the [] operator.
  If necessary, it asks vector for resize() 
*/
template <class T> class SafeVector : public vector<T> {
 public:
   typedef typename vector<T>::reference SafeVectorReference;
   typedef typename vector<T>::const_reference Const_SafeVectorReference;
   typedef typename list<T>::const_reference SafeListIterator;

  SafeVector<T>() : vector<T>() {}
  SafeVector<T>(const int i) : vector<T>(i) {}
  template <class InputIterator>
  SafeVector<T>(InputIterator i, InputIterator j) : vector<T>(i, j) {}
  SafeVector<T>(const list<T>& li) : vector<T>(li.size()) {
    for (SafeListIterator i = li.begin(); i != li.end(); i++) 
      push_back(*i);
  }

  SafeVectorReference operator[](size_t __n) {
    if (__n >= this->size() ) this->resize(__n + 1);
    return *(this->begin() + __n); 
  }
  Const_SafeVectorReference  operator[](size_t __n) const  { 
    if (__n >= this->size() ) throw SafeVectorOutOfRange(); // resize(__n + 1);
    return *(this->begin() + __n); 
  }
};

#endif

