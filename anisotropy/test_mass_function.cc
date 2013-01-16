/* mass_function.cc */
#include <iostream>
#include <math.h>
#include "../mathobject/spline.h"
#include "../mathobject/anchor.h"
#include "mass_function.h"

int main(){
	Anchor SplineAnchor;
	Spline* p_k = new Spline(1000, "power_spectrum", &SplineAnchor);
	// Fill spline
	MassFunction *mf = new MassFunction(p_k);
	return 0;
}


