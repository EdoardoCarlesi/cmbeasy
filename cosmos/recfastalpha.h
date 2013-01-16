#ifndef RECFASTALPHA_H
#define RECFASTALPHA_H

#include "recombination.h"

/*!
  Recombination routine for arbitrary fine structure constant.
  Numerically stable for some 10-30% change at the epoch of 
  recombination. 

  The changes made resemble the ones detailed in Kaplinghat et. al. (2000).
  
  I've experimented quite a lot to get stability and independence from
  assuming physics according to some redshift. Therefore RecfastAlpha
  looks quite a bit different from Recfast.
  
  Still, the original Fortran Code has been written by Douglas Scott (dscott@astro.ubc.ca)
  based on calculations in the papers Seager, Sasselov & Scott
  (ApJ, 523, L1, 1999; ApJS, 128, 407, 2000).

WARNING: this class has not been updated to recfast 1.4.2 like the standard Recfast

*/


class RecfastAlpha : public Recombination {
  bool Crossed98;
  bool CrossedHelium;
 

  double C; //!< speed of light
  double k_B; //!< boltzmann constant
  double mu_T; // !< total molecular weight
  double m_H; //! < mass of one hydrogen atom
  double H; //!< hubble in seconds;
  double OmegaB, OmegaT; //!< omega_baryons, omega_matter, 
  double x_H0, x_He0, x0; //!< set by getInit()
  double Tnow; //!< is initialized to t_cmb() today
  double H_frac;
  double fu; //!< fudge factor
  double fHe; //!< n_He_tot / n_H_tot

  double Nnow;
  
  double CR;
  double CB1_He1, CB1_He2; // !< // ionization for HeI & HeII
  double CB1;
  double CDB, CDB_He;
  double CK, CK_He;
  double CT;
  double CL, CL_He;
  double Bfact;
  double SwitchTbEvolution;
  void getInit(double z);
 public:
  RecfastAlpha(Cosmos& c, ControlPanel& control, Anchor* a=0) : Recombination(c, control, a)  {};
  void compute();
  void ion(const double, const double*,double*);
  double BoltzmannOverWeightAndC2() { return  k_B / (C*C * mu_T*m_H);} // k_b / (c^2 * mean_molecular weight of matter); (prefactor of eqn(68) Ma & Bertschinger 
  virtual double fineStructureConstant(double z); //!< re-implement your varying alpha here
};

#endif
