#ifndef RECOMBINATIONFACTORY_H
#define RECOMBINATIONFACTORY_H

#include "recfast.h"
#include "recfastalpha.h"

#include "quintcosmos.h"
#include "controlpanel.h"

#ifndef PRERELEASE
#include "recfastquint.h"
#endif

class RecombinationFactory {
 public:
  static Recombination* recombination(Cosmos& c, ControlPanel& control, Anchor* a=0) {
    Recombination::recombination r=control.recombination;
    switch (r) {
    case Recombination::RecfastAlpha:
      return new RecfastAlpha(c, control, a);
      break;
#ifndef PRERELEASE
    case Recombination::RecfastQuint:
      return new RecfastQuint(c, control, a);
      break;
#endif
    default:
      return new Recfast(c, control, a);
    }
  }
};

#endif 
