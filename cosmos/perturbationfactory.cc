#include "perturbationfactory.h"
#include "synchronous.h"
#include "quintsynchronous.h"
#include "speedyinvariant.h"
#include "quintcosmos.h"
#include "speedydeinvariant.h"
#include "coupledquintcosmos.h"
  #include "coupledinvariant.h"
  #include "speedycoupledinvariant.h"
  #include "cninvariant.h"

Perturbation*  PerturbationFactory::perturbation(Gauge::gauge gauge, Cosmos *cosmos) {
  switch (gauge) {
  case Gauge::synchronous:
    return new Synchronous(cosmos);
  case Gauge::quintSynchronous:
    return new QuintSynchronous((QuintCosmos*)cosmos);
  case Gauge::speedyInvariant:
    return new SpeedyInvariant(cosmos);
  case Gauge::speedyDEInvariant:
    return new SpeedyDEInvariant((QuintCosmos*)cosmos);
  case Gauge::cnInvariant:
     return new CnInvariant((CnCosmos*)cosmos);
  case Gauge::speedyCoupledInvariant:
     return new SpeedyCoupledInvariant((CoupledQuintCosmos*)cosmos);
  default:
    throw Bad_Error("Perturbationfactory:: gauge unknown");
  }
}

