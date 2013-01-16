#include "chainshop.h"
#include <iostream>
#include "global.h" 
using namespace std; 
  
int main(int argc, char* argv[]) {
  if (argc < 3) { 
    cerr << "Usage: distill file file file (...) outputfile" << endl << endl;
    cerr << "Two files will be generated: outputfile.dat and outputfil.mcc" << endl;
  }

  const char* out(argv[argc-1]);
  vector<string> namen(argc-2);
  for (int i =0; i < argc-2; i++) namen[i] = argv[i+1];

  try {
    ChainShop mine;
    mine.distillChain(namen,out);
  } catch (Bad_Error x) {
    cout << endl << endl << "BAD ERROR: " << endl << endl;
    cout << x.s << endl;
  }
  return 0;
}
 
 
 
  
  
