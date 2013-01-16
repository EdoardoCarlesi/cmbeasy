#ifndef RECFAST_H
#define RECFAST_H

#include "recombination.h"

/*!
  Slightly modified version of Recfast.
  Original Fortran Code written by Douglas Scott (dscott@astro.ubc.ca)
  based on calculations in the papers Seager, Sasselov & Scott
  (ApJ, 523, L1, 1999; ApJS, 128, 407, 2000).

  Updated to recfast 1.4.2 Jan '09,
  which includes the updates in
  Wong, Moss & Scott (2008).
  
  Updated to recfast 1.5 Jan '10.

  Differences to the original recfast:
       - calculate c_s^2 for baryons
       - take the background from cosmos, so e.g.  Hz is correct.

       The copyright notice of the original recfast.for:

  Copyright 1999-20010 by University of British Columbia.  All rights reserved.

  THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO 
  REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  
  BY WAY OF EXAMPLE, BUT NOT LIMITATION,
  U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF 
  MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT 
  THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE 
  ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.   


 Here is a copy of the description and the nicely documented list
 of parameters from the original version:

  Name:    RECFAST
  Version: 1.5

Purpose:  Calculate ionised fraction as a function of redshift.
          Solves for H and He simultaneously, and includes
          H "fudge factor" for low z effect, as well as
          HeI fudge factor.

  Description: Solves for ionisation history since recombination
  using the equations in Seager, Sasselov & Scott (ApJ, 1999).
  The Cosmological model can be flat or open.
  The matter temperature is also followed, with an update from
  Scott & Scott (2009).
  The values for \alpha_B for H are from Hummer (1994).
  The singlet HeI coefficient is a fit from the full code.
  Additional He "fudge factors" are as described in Wong, Moss
  and Scott (2008).
  Extra fitting function included (in optical depth) to account
  for extra H physics described in Rubino-Martin et al. (2010).
  Care is taken to use the most accurate constants.
  Note that some parameters are fixed (e.g. N_nu=3, nu's are
  massless, w=-1, etc.) - some users may want to explictly
  imput their own H(z) to account for extra physics.
  This is provided as a PROGRAM, which can be easily converted
  to a SUBROUTINE for use in CMB Boltzmann codes.

  z is redshift - W is sqrt(1+z), like conformal time
  x is total ionised fraction, relative to H
  x_H is ionized fraction of H - y(1) in R-K routine
  x_He is ionized fraction of He - y(2) in R-K routine
  (note that x_He=n_He+/n_He here and not n_He+/n_H)
  Tmat is matter temperature - y(3) in R-K routine
  f's are the derivatives of the Y's
  alphaB is case B recombination rate
  alpHe is the singlet only HeII recombination rate
  a_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
  b_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
  c_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
  d_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
  a_VF is Verner and Ferland type fitting parameter for Helium
  b_VF is Verner and Ferland type fitting parameter for Helium
  T_0 is Verner and Ferland type fitting parameter for Helium
  T_1 is Verner and Ferland type fitting parameter for Helium
  Tnow is the observed CMB temperature today
  OmegaT is the total Omega_0
  OmegaL is the Omega_0 contribution from a Cosmological constant
  OmegaK is the Omega_0 contribution in curvature (1-O_T-O_L)
  OmegaB is Omega in baryons today
  OmegaC is the Omega_0 in (cold) dark matter: OmegaT=OmegaC+OmegaB
  Yp is the primordial helium abundace
  fHe is He/H number ratio = Yp/4(1-Yp)
  Trad and Tmat are radiation and matter temperatures
  epsilon is the approximate difference (=Trad-Tmat) at high z
  OmegaB is Omega in baryons today
  H is Hubble constant in units of 100 km/s/Mpc
  HOinp is input value of Hubble constant in units of 100 km/s/Mpc
  HO is Hubble constant in SI units
  bigH is 100 km/s/Mpc in SI units
  Hz is the value of H at the specific z (in ION)
  G is grvitational constant
  n is number density of hydrogen
  Nnow is number density today
  x0 is initial ionized fraction
  x_H0 is initial ionized fraction of Hydrogen
  x_He0 is initial ionized fraction of Helium
  rhs is dummy for calculating x0
  zinitial and zfinal are starting and ending redshifts
  fnu is the contribution of neutrinos to the radn. energy density
  zeq is the redshift of matter-radiation equality
  zstart and zend are for each pass to the integrator
  w0 and w1 are conformal-time-like initial and final zi and zf's
  Lw0 and Lw1 are logs of w0 and w1
  hw is the interval in W
  C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
  m_e,m_H: electron mass and H atomic mass in SI
  not4: ratio of 4He atomic mass to 1H atomic mass
  sigma: Thomson cross-section
  a: radiation constant for u=aT^4
  Pi: Pi
  Lambda: 2s-1s two photon rate for Hydrogen
  Lambda_He: 2s-1s two photon rate for Helium
  DeltaB: energy of first excited state from continuum = 3.4eV
  DeltaB_He: energy of first excited state from cont. for He = 3.4eV
  L_H_ion: level for H ionization in m^-1
  L_H_alpha: level for H Ly alpha in m^-1
  L_He1_ion: level for HeI ionization
  L_He2_ion: level for HeII ionization
  L_He_2s: level for HeI 2s
  L_He_2p: level for He 2p (21P1-11S0) in m^-1
  Lalpha: Ly alpha wavelength in SI
  Lalpha_He: Helium I 2p-1s wavelength in SI
  mu_H,mu_T: mass per H atom and mass per particle
  H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
  dHdz is the derivative of H at the specific z (in ION)
  CDB=DeltaB/k_B  Constants derived from B1,B2,R
  CDB_He=DeltaB_He/k_B         n=2-infinity for He in Kelvin
  CB1=CDB*4.      Lalpha and sigma_Th, calculated
  CB1_He1: CB1 for HeI ionization potential
  CB1_He2: CB1 for HeII ionization potential
  CR=2*Pi*(m_e/h_P)*(k_B/h_P)	once and passed in a common block
  CK=Lalpha**3/(8.*Pi)
  CK_He=Lalpha_He**3/(8.*Pi)
  CL=C*h_P/(k_B*Lalpha)
  CL_He=C*h_P/(k_B*Lalpha_He)
  CT=(8./3.)*(sigma/(m_e*C))*a
  Bfact=exp((E_2p-E_2s)/kT)	Extra Boltzmann factor
  fu is a "fudge factor" for H, to approximate low z behaviour
  b_He is a "fudge factor" for HeI, to approximate higher z behaviour
  Heswitch is an integer for modifying HeI recombination
           Parameters and quantities to describe the extra triplet states
           and also the continuum opacity of H, with a fitting function
           suggested by KIV, astro-ph/0703438
  a_trip: used to fit HeI triplet recombination rate
  b_trip: used to fit HeI triplet recombination rate
  L_He_2Pt: level for 23P012-11S0 in m^-1
  L_He_2St: level for 23S1-11S0 in m^-1
  L_He2St_ion: level for 23S1-continuum in m^-1
  A2P_s: Einstein A coefficient for He 21P1-11S0
  A2P_t: Einstein A coefficient for He 23P1-11S0    
  sigma_He_2Ps: H ionization x-section at HeI 21P1-11S0 freq. in m^2
  sigma_He_2Pt: H ionization x-section at HeI 23P1-11S0 freq. in m^2
  CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
  CfHe_t: triplet statistical correction
  Hswitch is an integer for modifying the H recombination
  AGauss1 is the amplitude of the 1st Gaussian for the H fudging
  AGauss2 is the amplitude of the 2nd Gaussian for the H fudging
  zGauss1 is the ln(1+z) central value of the 1st Gaussian
  zGauss2 is the ln(1+z) central value of the 2nd Gaussian
  wGauss1 is the width of the 1st Gaussian
  wGauss2 is the width of the 2nd Gaussian
  tol: tolerance for the integrator
  cw(24),w(3,9): work space for DVERK
  Ndim: number of d.e.'s to solve (integer)
  Nz: number of output redshitf (integer)
  I: loop index (integer)
     ind,nw: work-space for DVERK (integer)
*/


class Recfast : public Recombination
{
  double C; //!< speed of light
  double k_B; //!< boltzmann constant
  double mu_T; //!< total molecular weight
  double h_P; //!<  Planck's constants
  double m_H; //!< mass of one hydrogen atom
  double not4; //!< ratio of 4He atomic mass to 1H atomic mass
  double H; //!< hubble in seconds;
  double OmegaB, OmegaT; //!< omega_baryons, omega_matter, 
  double x_H0, x_He0, x0; //!< set by getInit()
  double Tnow; //!< is initialized to t_cmb() today
  double H_frac;
  double fu; //!< fudge factor
  double b_He; //!< "fudge factor" for HeI
  double fHe; //!< n_He_tot / n_H_tot
  double L_He_2p; //!< level for He 2p (21P1-11S0) in m^-1

  double Nnow;

  double CR;
  double CB1_He1, CB1_He2; //!< ionization for HeI & HeII
  double CB1;
  double CDB, CDB_He;
  double CK, CK_He;
  double CT;
  double CL, CL_He;
  double Bfact;
  double SwitchTbEvolution;

  Spline* mdHdz;   //!< spline to hold dH/dz

  void getInit(double z);

 public:
  Recfast(Cosmos& c, ControlPanel& control, Anchor* a=0) : Recombination(c, control, a), mdHdz(0)  {}
  void compute();
  void ion(const double, const double*,double*);

  // k_b / (c^2 * mean_molecular weight of matter);
  // (prefactor of eqn(68) Ma & Bertschinger 
  double BoltzmannOverWeightAndC2() { return  k_B / (C*C * mu_T*m_H); }

  /*! options for HeI recombnaion from RecFast (Heswitch)
   *
   *  0) no change from old Recfast'
   *  1) full expression for escape probability for singlet 1P-1S transition
   *  2) also including effect of contiuum opacity of H on HeI singlet
   *     (based in fitting formula suggested by Kholupenko, Ivanchik & Varshalovich, 2007)
   *  3) only including recombination through the triplets
   *  4) including 3 and the effect of the contiuum
   *     (although this is probably negligible)
   *  5) including only 1, 2 and 3
   *  6) including all of 1 to 4
   */
   typedef unsigned int HeliumRecombinationSwitch;
};
#endif
