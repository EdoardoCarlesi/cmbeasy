#include "anchor.h"

void Anchor::kill() {
  // cout << "---- ANCHOR HAS BEEN ASKED TO KILL --"<<endl;
  bool isEmpty;
#ifdef WITH_OMP
#pragma omp critical (Anchor)
#endif
  {
    isEmpty = obj.empty();
  }
  while (!isEmpty) {
    //cout << "---- ANCHOR IS killing --   "<< obj.begin()->second  << endl;
    delete obj.begin()->second;
    //cout << "finished " << endl;
#ifdef WITH_OMP
#pragma omp critical (Anchor)
#endif
    {
      isEmpty = obj.empty();
    }
  }
}

void Anchor::printStatus(ostream &o) {
  o << "Anchor content " << endl;
  for (map<AnchorEnabled*,AnchorEnabled*>::iterator i = obj.begin(); i != obj.end(); i++) 
    o << i->first << endl;
  o<< "===== end of anchor == (total count in class: "<< count << ")" << endl;
}

int Anchor::count = 0;
