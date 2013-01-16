#ifndef PERTURBATIONFACTORY_H
#define PERTURBATIONFACTORY_H 

#include "gauge.h"

class Perturbation;
class Cosmos;
class QuintCosmos;

/*!
  Small class with one single static function and objective:

  create a perturbation object of a given gauge and return
  a pointer to it.

*/
class PerturbationFactory {
 public:
  PerturbationFactory() {};

  static Perturbation* perturbation(Gauge::gauge, Cosmos*);

};

#endif
