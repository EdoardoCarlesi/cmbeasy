#ifndef GLOBAL_H
#define GLOBAL_H

using namespace std;
#include <string>
/*!
  Main Error class. Throughout the 
  package this exception is thrown 
  whenever there is an error. Usually
  main() or the gui catches it and
  displays the error message (and
  stops the program, if necessary). 

  It is also quite usefull for debugging
  or data-extraction (if you just need
  the first few k-values for example).

  Just insert a line "throw Bad_Error("stop");"
  and the program will terminate at the
  point you like.
*/
struct Bad_Error {
  string s;
  Bad_Error(string q="") {s  =q; };
};

#ifndef NO_WARNINGS
#define WARN(m) cerr << endl << "[WARNING] " <<   m  << endl << endl
#endif 

#ifdef NO_WARNINGS
#define WARN(m) 
#endif

#endif
