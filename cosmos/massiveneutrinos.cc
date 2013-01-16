// vim: sw=2 ts=2 noet expandtab


#include "massiveneutrinos.h"

#include "miscmath.h"

#include <cmath>
#include <sstream>

#define QMAXI 30.0 //! From CMBFAST
#define NQ 1000 //! From CMBFAST

static bool initialized = false;
static unsigned int firstIndex;

static const double qmax=10000;

namespace MassiveNeutrinos {
  static double* qdn = 0;
  static Nu1d_ nu1d;
  static double* dlfdlq = 0;
  static double* denl;
//X   static const unsigned int mqGridSize = 15;
  static const unsigned int nqmax = 20;
  static const unsigned int mNuNRlMax = 15;
//X   static const unsigned int nqmax = 25;
//X   static const unsigned int mNuNRlMax = 50;
}

void MassiveNeutrinos::init()
{
  initialized = true;
  qdn = new double[nqmax];
  double dq = 1.;
  double q;
  for (unsigned int i = 0; i <  nqmax; ++i) {
    q = i + .5;
    qdn[i] = dq * q * q * q / (exp(q) + 1.);
  }

  // dlnf_0(q) / dlnq  where f(q) = 1 / (exp(q) + 1)
  dlfdlq = new double[nqmax];  //!< dlnf_0 / dlnq array up to nqmax-1
  for (unsigned int i = 0; i < nqmax; ++i) {
    double q = i * dq + .5;
    dlfdlq[i] = -q / (exp(-q)  + 1.);
  }

  denl = new double[mNuNRlMax+1];
  for (unsigned int j = 1; j <= mNuNRlMax; ++j) {
    denl[j] = 1. / (double) ((j << 1) + 1);
  }
}

void MassiveNeutrinos::setFirstIndex(unsigned int i)
{
  firstIndex = i;
}

unsigned int MassiveNeutrinos::pertArraySize()
{
  return nqmax*(mNuNRlMax+1);
}

unsigned int MassiveNeutrinos::qGridSize()
{
  return nqmax;
}


/*
MassiveNeutrinos::~MassiveNeutrinos()
{
  delete qdn;
  qdn=0;
}
*/

/*!
  Initialize interpolation tables for massive neutrinos.
  Use cubic splines interpolation of log rhonu and pnu vs. log a.
*/

void MassiveNeutrinos::initnu1(const double amnu)
{
    double a;
    int i;
    double rhonu, pnu;

    nu1d.amin = 1e-9;
    nu1d.amax=10000.;
    nu1d.dlna = -log(nu1d.amin/nu1d.amax) / 9999;

    for (i = 0; i <  NRHOPN; ++i) {
      a = nu1d.amin * exp(i * nu1d.dlna);
      ninu1(a, &rhonu, &pnu, amnu);
      nu1d.r1[i] = log(rhonu);
      nu1d.p1[i] = log(pnu);
    }

    Miscmath::splini();
    Miscmath::splder(nu1d.r1, nu1d.dr1, NRHOPN);
    Miscmath::splder(nu1d.p1, nu1d.dp1, NRHOPN);
    Miscmath::splder(nu1d.dr1, nu1d.ddr1, NRHOPN);
}


void MassiveNeutrinos::propagateMassiveNeutrinoMomentsLongitudinal(const double *y, double *yprime,
                                                                   const double a, const double tau,
                                                                   const double nuMass, const double k,
                                                                   const double PhiDot, const double Psi,
                                                                   const double beta, const double X)
{
  double *Nmprime = &yprime[firstIndex];
  const double *NR_0 = &y[firstIndex];
  const double *NR_1 = &NR_0[nqmax];
  const double *NR = &NR_1[nqmax];
  double *NRprime_0 = Nmprime;
  double *NRprime_1 = &NRprime_0[nqmax];
  double *NRprime = &NRprime_1[nqmax];
  double dq = 1.;
  double q,aq,v;
  double akv[nqmax]; // this is q k / epsilon in M&B
  double aq2table[nqmax];
  for (unsigned int i = 0; i < nqmax ; ++i) {
    q = i * dq + .5;
    aq = a *nuMass / q;
    aq2table[i] = aq*aq;
    v = 1. / sqrt(aq * aq + 1.);
    akv[i] = k * v;   // this is q*k / sqrt( a^2 m^2 + q^2)
  }

  double ell;
  //  l = 0, 1, 2,lmax_nu_NR.
  for (unsigned int i = 0; i < nqmax  ; ++i) {
    NRprime_0[i] = -akv[i] * NR_1[i] + PhiDot *dlfdlq[i];
    NRprime_1[i] = ( akv[i] * (NR_0[i] - 2.*NR[i])
        - k*k*Psi*dlfdlq[i]/akv[i] )/3.
        - 1./3.*dlfdlq[i]*X*beta*akv[i]*aq2table[i];
    NRprime[i] = akv[i]*denl[2]*( 2.*NR_1[i] - 3.*NR[nqmax + i] );
      for (unsigned int l = 3; l < mNuNRlMax ; ++l) {
        ell = (double ) l;
        NRprime[(l-2)*nqmax + i] =
          akv[i]*denl[l]*( ell*NR[(l-3)*nqmax + i] -
                          (ell+1.)*NR[(l-1)*nqmax + i]);
      }
    //  Truncate moment expansion.
      NRprime[(mNuNRlMax - 2)*nqmax + i] = akv[i] * NR[(mNuNRlMax-3)*nqmax + i]
                                  - ( (double)(mNuNRlMax + 1) ) / tau * NR[(mNuNRlMax-2)*nqmax + i];

  }

}

void MassiveNeutrinos::propagateMassiveNeutrinoMomentsSynchronous(const double *y, double *yprime,
                                                                  const double a, const double tau,
                                                                  const double nuMass, const double k,
                                                                  const double hdot, const double etadot)
{
 //keep old cmbfast notation
 unsigned int iq0 = firstIndex;
 unsigned int iq1 = iq0 + nqmax;
 unsigned int iq2 = iq1 + nqmax;
 double dq = 1.;
 double q,aq,v;
 double akv[nqmax];
  for (unsigned int i = 0; i < nqmax ; ++i) {
    q = i * dq + .5;
    aq = a * nuMass / q;
    v = 1. / sqrt(aq * aq + 1.);
    akv[i] = k * v;   // this is q*k / sqrt( a^2 m^2 + q^2)
  }
  //  l = 0, 1, 2,lmaxnu.
  for (unsigned int i = 0; i < nqmax  ; ++i) {
    unsigned int ind = iq0 + i;
    yprime[ind] = -akv[i] * y[ind + nqmax] + hdot * dlfdlq[i] / 6.;
    ind = iq1 + i;
    yprime[ind] = akv[i] * (y[ind - nqmax] - y[ind + nqmax] * 2) / 3;
    ind = iq2 + i;
    yprime[ind] = akv[i] * (y[ind - nqmax] * 2 - y[ind + nqmax] * 3.0) / 5.0 
      - (hdot / 15. + 0.4 * etadot) * dlfdlq[i];
    ind = firstIndex + i + mNuNRlMax * nqmax;
    //  Truncate moment expansion.
    yprime[ind] = akv[i] * y[ind - nqmax] - (mNuNRlMax + 1) / tau * y[ind];
  }
  for (unsigned int l = 3; l < mNuNRlMax ; ++l) {
    double ell = double(l);
    for (unsigned int i = 0; i < nqmax ; ++i) {
      unsigned int ind = firstIndex + i + l * nqmax;
      yprime[ind] = akv[i] * denl[l] * (ell * y[ind - nqmax] - (ell + 1.) * y[ind + nqmax]);
    }
  }
}

static std::list<double> debugScaleFactors = std::list<double>();

void MassiveNeutrinos::nu2(const double a, const double betanu,
                           const double rhonu, const double pnu, const double deltaphi,
                           double *drhonu, double *fnu, double *dpnu, double *shearnu,
                           const double *psi0, const double * psi1, const double *psi2)
{
    if (!initialized)
      init();

    double q, v, g0[4], g1[nqmax+2], g2[nqmax+2], g3[nqmax+2], g4[nqmax+2], aq;
    double dum3[nqmax+2];
    unsigned int iq;
    double gf1, gf2, gf3, gf4;


    //  Compute the perturbations of density, energy flux, pressure, and 
    //  shear stress of one flavor of massive neutrinos, in units of the mean 
    //
    //  density of one flavor of massless neutrinos, by integrating over 
    //  momentum. 

    //  const=7*pi**4/120.  = 5.68219698
    static const double Const = 7./120.*pow(M_PI, 4.);

    if (nqmax == 0) {
      *drhonu = 0.;
      *fnu = 0.;
      *dpnu = 0.;
      *shearnu = 0.;
      return;
    }

    //  q is the comoving momentum in units of k_B*T_nu0/c. 
    g1[0] = 0.;
    g2[0] = 0.;
    g3[0] = 0.;
    g4[0] = 0.;
    dum3[0] = 0.;

    for (iq = 2; iq <= nqmax+1; ++iq) {
      q = iq -  1.5;
      aq = a / q;
      v = 1. / sqrt(aq * aq + 1.);
      g1[iq-1] = qdn[iq - 2] * psi0[iq-2] / v;
      g2[iq-1] = qdn[iq - 2] * psi0[iq-2] * v;
      g3[iq-1] = qdn[iq - 2] * psi1[iq-2];
      g4[iq-1] = qdn[iq - 2] * psi2[iq-2] * v; 
      dum3[iq-1]=qdn[iq - 2]*v*v*v*aq*aq;
    }

    double tempant;
    Miscmath::splint(dum3, &tempant, nqmax+1);

    tempant=(tempant+dum3[nqmax]*2./(nqmax-0.5))/Const;
    tempant=tempant*deltaphi*betanu;

    Miscmath::splint(g1, g0, nqmax+1);
    Miscmath::splint(g2, &g0[1], nqmax+1);
    Miscmath::splint(g3, &g0[2], nqmax+1);
    Miscmath::splint(g4, &g0[3], nqmax+1);
    gf1 = g1[nqmax];
    gf2 = g2[nqmax];
    gf3 = g3[nqmax];
    gf4 = g4[nqmax];
    *drhonu = (g0[0] + gf1 * 2. / (nqmax-0.5)) / 5.68219698+(rhonu-3.*pnu)*deltaphi*betanu;
    *dpnu = (g0[1] + gf2 * 2. / (nqmax-0.5)) / 5.68219698 / 3.;
    *dpnu -= tempant/3.;
    *fnu = (g0[2] + gf3 * 2. / (nqmax-0.5)) / 5.68219698;
    *shearnu = (g0[3] + gf4 * 2. / (nqmax-0.5)) / 5.68219698 * 2. / 3.;

}

double MassiveNeutrinos::nuder(const double a, const double adotoa, const double beta,
                               const double rhonu, const double phidot,
                               const double *psi2, const double *psi2dot)
{
    if (!initialized)
      init();

    static double vdot;

    double q, v, aqdot, g0, g1[nqmax+1], aq;
    static double gf1;

    //  Compute the time derivative of the mean density in massive neutrinos 
    //
    //  and the shear perturbation. 


    if (nqmax == 0) {
      double shearnudot = 0.;
      return shearnudot;
    }

    //  q is the comoving momentum in units of k_B*T_nu0/c. 
    g1[0] = 0.;
    for (unsigned int iq = 1; iq <= nqmax; ++iq) {
      q = iq - .5;
      aq = a/q;
      aqdot = aq * (adotoa + beta*phidot);
      v = 1. / sqrt(aq * aq + 1.);
      vdot = -aq * aqdot / pow(aq * aq + 1. , 1.5);
      g1[iq] = qdn[iq - 1] * (psi2dot[iq-1] * v + psi2[iq-1] * vdot);
    }
    Miscmath::splint(g1, &g0, nqmax+1);
    gf1 = g1[nqmax];
    double shearnudot = (g0 + gf1 * 2. / (nqmax-0.5)) / 5.68219698 * 2. / 3.;

    return shearnudot;
}

void MassiveNeutrinos::ninu1(const double a, double *rhonu, double *pnu,  const double amnu)
{
  double q, v, aq, dq;
  static double qdn, dum1[NQ+1], dum2[NQ+1];

  //  const=7*pi**4/120.  = 5.68...

  //  q is the comoving momentum in units of k_B*T_nu0/c.
  //  Integrate up to qmax and then use asymptotic expansion for remainder.

  dq = QMAXI/NQ;
  dum1[0] = 0.;
  dum2[0] = 0.;

  // see Eqn. (52) of (Ma & Bertschinger)
  for (int i = 1; i <= NQ; ++i) {
    q = i * dq;
    // aq = a * amnu / q;
    // changed for coupled neutrinos
    aq = a / q;
    v = 1. / sqrt(aq * aq + 1.);
    qdn = dq * q * q * q / (exp(q) + 1.);
    dum1[i] = qdn / v;
    dum2[i] = qdn * v;
  }
  Miscmath::splint(dum1, rhonu, NQ + 1);
  Miscmath::splint(dum2, pnu, NQ + 1);

  //  Apply asymptotic corrrection for q>qmax and normalize by relativistic
  //  energy density.
  *rhonu = (*rhonu + dum1[1000] / dq) / 5.68219698;
  *pnu = (*pnu + dum2[1000] / dq) / 5.68219698 / 3.;
} // ninu1_

void MassiveNeutrinos::nu1(const double a, double *rhonu, double *pnu)
{
  double d;
  int i;

  // special case of a=0 (for coupled neutrinos, a=a*mass_nu(phi)
  // which can be 0
  if ( a == 0 ) {
    *rhonu = 1.;
    *pnu = 1./3.;
    return;
  }


  d = log(a / nu1d.amin) / nu1d.dlna + 1.;
  i = (int) d;
  d -= i;
  if (i < 1) {
    //  Use linear interpolation, bounded by results for massless neutrinos. 
    *rhonu = nu1d.r1[0] + (d - 1) * nu1d.dr1[0];
    *pnu = nu1d.p1[0] + (d - 1) * nu1d.dp1[0];
    *rhonu = std::min( exp(*rhonu) ,1.);
    *pnu = std::min( exp(*pnu) ,.3333333333);
  } else if (i >= NRHOPN) {
    //  This should not happen, unless the user evolves to z<0! 
    *rhonu = nu1d.r1[NRHOPN-1] + (d + i - 10000) * nu1d.dr1[NRHOPN-1];
    *pnu = nu1d.p1[NRHOPN-1] + (d + i - 10000) * nu1d.dp1[NRHOPN-1];
    *rhonu = exp(*rhonu);
    *pnu = exp(*pnu);
  } else {
    //  Cubic spline interpolation.
    *rhonu = nu1d.r1[i - 1] + d * (nu1d.dr1[i - 1] + d * ((nu1d.r1[i] - nu1d.r1[i - 1]) * 3. - nu1d.dr1[i - 1] * 2. - nu1d.dr1[i] + d * (nu1d.dr1[i - 1] + nu1d.dr1[i] + (nu1d.r1[i - 1] - nu1d.r1[i]) * 2.)));
    
    *pnu = nu1d.p1[i - 1] + d * (nu1d.dp1[i - 1] + d * ((nu1d.p1[i] - nu1d.p1[i - 1]) * 3. - nu1d.dp1[i - 1] * 2. - nu1d.dp1[i] + d * (nu1d.dp1[i - 1] + nu1d.dp1[i] + (nu1d.p1[i - 1] - nu1d.p1[i]) * 2.)));
    
    *rhonu = exp(*rhonu);
    *pnu = exp(*pnu);
  }
} // nu1_


void MassiveNeutrinos::setIninitalLongitudinalScalarPerturbations(double *y, const double a,
                                        const double mass_nu, const double deltan,
                                        const double Vn, const double Pi_nu)
{
  double *NR_0 = &y[firstIndex]; // massive neutrinos
  double *NR_1 = &NR_0[nqmax];
  double *NR_2 = &NR_1[nqmax];
  double *NR = &NR_2[nqmax];
  double dq = 1.;
  double aq,q,dlfdlq;
  for (unsigned int i = 0; i < nqmax; ++i) {
    q = (i+1) * dq - .5;
    aq = a * mass_nu / q;
    dlfdlq = -q / ( exp(-q) + 1.);
    NR_0[i] = -0.25*dlfdlq * deltan;
    NR_1[i] = -dlfdlq * Vn* sqrt(aq*aq+1.) / 3.;
    NR_2[i] = -dlfdlq * Pi_nu/12.;
    for (unsigned int l = 3; l <= mNuNRlMax; ++l) {
      NR[(l-3)*nqmax + i ] = 0.;
    }
  }
}
