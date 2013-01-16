#ifndef CLEANVECTOR
#define CLEANVECTOR
// vim:ts=2:sw=2:et

#include <vector>
#include <iostream>

using namespace std;

/*!
  A slight modification of the standard vector. The added
  functionality is automatic deleting of all data. Hence, CleanVector
  only makes sense to store pointers.
  */
template <class T> class CleanVector : public vector<T> {
  public:
    typedef typename vector<T>::iterator VectorIterator;
    typedef typename vector<T>::const_iterator ConstVectorIterator;
    CleanVector<T>() : vector<T>() {}
    CleanVector<T>(const int i) : vector<T>(i) {}
    ~CleanVector() {
      //cout << "+++++++++++++++++++++ CLEANVECTOR DESTRUCTION +++++++"<<endl;
      clear();
      //cout << "++ACCOMPLISHED"<<endl;
    }

    void clear() {
      for (VectorIterator i=this->begin(); i!=this->end(); i++) {
        delete (*i);
        (*i)=0;
      }
      vector<T>::clear();
    }

    static void clear(vector<T> &v) {
      v.clear();
    }

    void resize(size_t neu) {
      if (neu < this->size() && neu >0) {
        for (unsigned int i=neu; i<this->size(); i++) {
          T now=((*this)[i]);
          delete now;
        }
      }
      vector<T>::resize(neu);
    }

    void printStatus() const {
      int count=0;
      for (ConstVectorIterator i=this->begin(); i !=this->end(); i++)
        cout << count++ <<  "   " << *i << endl;
    }

};

#endif

