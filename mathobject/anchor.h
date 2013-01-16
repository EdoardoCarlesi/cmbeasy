#ifndef ANCHOR_H
#define ANCHOR_H

using namespace std;
#include <map>
#include <iostream>
#include "global.h"
class AnchorEnabled;  // defined below

/*!
  Anchor provides basic functionalities to 
  keep track of objects (i.e. call delete for these
  objects at some point of time to prevent
  memory-leaks)
  
  You can:
  
  add an object to the anchor and 
  remove objects from it:

  \code
  Anchor a;
  Spline *s;  // Splines are publicly derived from AnchorEnabled
  a.add(s);  // add it 
  a.remove(s);  // remove it
  \endcode

  To get rid of all the objects in an anchor,
  you can either call explicitly
  \code 
  a.kill();
  \endcode
  Or the anchore does it on destruction.

  An anchored object will remove itself from
  the anchor on destruction. This prevents
  segfaults with double freeing.

  Only objects that are subclasses of AnchorEnabled
  can be added to an Anchor. This makes it very
  easy for an object to gain Anchor functionality. Just
  inherit, call the AnchorEnable() constructor with
  a given anchor as argument and whooops: The objects 
  is  looked after by an anchor. See for instance
  the Spline implementation for an example.
*/

class Anchor {
 private:
  map<AnchorEnabled*,AnchorEnabled*> obj;

 public:
  static int count;
  ~Anchor() { kill(); }
  //! add object to watchlist anchor
  void add(AnchorEnabled* s) {
#pragma omp critical (Anchor)
    {
      obj[s]=s; count++;
    }
  }
  //! remove spline from watchlist of anchor
  void remove(AnchorEnabled* s) {
#pragma omp critical (Anchor)
    {
      obj.erase(obj.find(s));
      count--;
    }
  }
  //! kill whole watchlist
  void kill();
  void printStatus(ostream &o=cout);

};


/*! 
  AnchoEnabled is just a struct that is abstract and
  needed for Anchor.
  Each Object that is derived from AnchorEnabled
  may be kept track in Anchor.
  Also, Constructing an AnchorEnabled Object with
  a none-null pointer to an Anchor will result in a call
  of add(this) to the Anchor in question
 */
class AnchorEnabled {
 public:
  Anchor *MyAnchor; 
  AnchorEnabled(Anchor* a) : MyAnchor(a) { 
    if (a) a->add(this);
  };   
  virtual ~AnchorEnabled() { 
    //cout << "{{{{{ ~AnchorEnabled" << endl;
    if (MyAnchor) { 
      //cout << "~AnchorEnabled() " << endl;
      MyAnchor->remove(this); 
    }
  }
};

#endif 
