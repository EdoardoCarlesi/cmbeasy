#include <string>
#include <fstream>
#include <math.h>
#include "controlpanel.h"
#include "cmbcalc.h"

int main() {
  int lmoin,kmax0;
  cout << "Maximum value of l, keta: (1500, 3000)" << endl;
  cout << "l: ";
  cin >> lmoin;
  cout << "keta: ";
  cin >> kmax0;
  CmbCalc::precalculateBesselFunctions(ControlPanel::cmbeasyDir("/resources/jlgen.dat"),lmoin,kmax0);
  return 0;
}

