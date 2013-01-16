// vim: sw=2 ts=2 noet expandtab


#include "massiveneutrinosnew.h"

#include "miscmath.h"
#include "cncosmos.h"
#include <iostream>
#include <cmath>

#warning for debug
#include <sstream>

#define QMAXI 30.0 //! From CMBFAST
#define NQ 1000 //! From CMBFAST

static bool initialized = false;
static unsigned int firstIndex;

static const double qmax=100; // original value: 10000, changed to speed up calculation

namespace MassiveNeutrinosNew {
  static double* qdn = 0;
  static Nu1d_ nu1d;
  static double* dlfdlq = 0;
  static double* denl;
//X   static const unsigned int mqGridSize = 15;
  static const unsigned int nqmax = 20;
  static const unsigned int mNuNRlMax = 15; //50
}

void MassiveNeutrinosNew::init()
{
#warning check dq
  initialized = true;
  qdn = new double[nqmax];
  double dq = 1.;
  double q;
  for (int i = 0; i <  nqmax; ++i) {
    q = i + .5;
    qdn[i] = dq * q * q * q / (exp(q) + 1.);
  }

  // dlnf_0(q) / dlnq  where f(q) = 1 / (exp(q) + 1)
  dlfdlq = new double[nqmax];  //!< dlnf_0 / dlnq array up to nqmax-1
  for (int i = 0; i < nqmax; ++i) {
    double q = i * dq + .5;
    dlfdlq[i] = -q / (exp(-q)  + 1.);
  }

  denl = new double[mNuNRlMax+1];
  for (int j = 1; j <= mNuNRlMax; ++j) {
    denl[j] = 1. / (double) ((j << 1) + 1);
  }
}

void MassiveNeutrinosNew::setFirstIndex(unsigned int i)
{
  firstIndex = i;
  //cout << "First Index set to: " << firstIndex << endl;
}

unsigned int MassiveNeutrinosNew::pertArraySize()
{
  return nqmax*(mNuNRlMax+1);
}

unsigned int MassiveNeutrinosNew::qGridSize()
{
  return nqmax;
}


/*
MassiveNeutrinosNew::~MassiveNeutrinosNew()
{
  delete qdn;
  qdn=0;
}
*/

/*!
  Initialize interpolation tables for massive neutrinos.
  Use cubic splines interpolation of log rhonu and pnu vs. log a.

  In principle, reading the data from disk is not a speed improvement, as
  reading the data is as costy as caculating it, as it seem.

  As it could lead to inconsitencies, we just skip it.

 */

void MassiveNeutrinosNew::initnu1(const double amnu)
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


void MassiveNeutrinosNew::propagateMassiveNeutrinoMoments(const double *y, double *yprime,
                                             const double a, const double tau,
                                             const double nuMass, const double k,
                                             const double PhiDot, const double Psi,
                                             const double beta, const double X, bool stop)
{
  double *Nmprime = &yprime[firstIndex];
  const double *NR_0 = &y[firstIndex];
  const double *NR_1 = &NR_0[nqmax];
  const double *NR = &NR_1[nqmax];
  double *NRprime_0 = Nmprime;
  double *NRprime_1 = &NRprime_0[nqmax];
  double *NRprime = &NRprime_1[nqmax];
#warning qspace
//X   double dq=qmax/nqmax;
  double dq = 1.;
  double q,aq,v;
  double akv[nqmax]; // this is q k / epsilon in M&B
  double aq2table[nqmax];
  for (int i = 0; i < nqmax ; ++i) {
    q = i * dq + .5;
    aq = a *nuMass / q;
    aq2table[i] = aq*aq;
    v = 1. / sqrt(aq * aq + 1.);
    akv[i] = k * v;   // this is q*k / sqrt( a^2 m^2 + q^2)
  }

  //  cout<< "beta: " << beta << endl;

  double ell;
  //  l = 0, 1, 2,lmax_nu_NR.
  ofstream *psi0, *psi1, *psi2;
  psi0=psi1=psi2=0;
  static bool wroteInt = false;
  static const double zStop = 1000;
  if ( false
    && !wroteInt && !( (1./a-1.)>zStop))
  {
    wroteInt=true;
    psi0 = new ofstream("psi0");
    psi1 = new ofstream("psi1");
    psi2 = new ofstream("psi2");
  }
  for (int i = 0; i < nqmax  ; ++i) {
    NRprime_0[i] = -akv[i] * NR_1[i] + PhiDot *dlfdlq[i];
    NRprime_1[i] = ( akv[i] * (NR_0[i] - 2.*NR[i])
        - k*k*Psi*dlfdlq[i]/akv[i] )/3.
        - 1./3.*dlfdlq[i]*X*beta*akv[i]*aq2table[i];
    NRprime[i] = akv[i]*denl[2]*( 2.*NR_1[i] - 3.*NR[nqmax + i] );
      for (int l = 3; l < mNuNRlMax ; ++l) {
        ell = (double ) l;
        NRprime[(l-2)*nqmax + i] =
          akv[i]*denl[l]*( ell*NR[(l-3)*nqmax + i] -
                          (ell+1.)*NR[(l-1)*nqmax + i]);
      }
    //  Truncate moment expansion.
      NRprime[(mNuNRlMax - 2)*nqmax + i] = akv[i] * NR[(mNuNRlMax-3)*nqmax + i]
                                  - ( (double)(mNuNRlMax + 1) ) / tau * NR[(mNuNRlMax-2)*nqmax + i];
    if (psi0 && psi1 &&psi2) {
      q = i * dq + .5;
      (*psi0) << q << "     " << NR_0[i] << endl;
      (*psi1) << q << "     " << NR_1[i] << endl;
      (*psi2) << q << "     " << NR[i] << endl;
    }
  }
  if (psi0 && psi1 &&psi2) {
//X     cout << "wrote psi" << endl;
    psi0->close();
    psi1->close();
    psi2->close();
  }

  if (stop) { // control higher neutrino moments cutoff
    for (int i = 0; i < nqmax  ; ++i) {
      NRprime_0[i] = 0;
      NRprime_1[i] = 0;
      NRprime[i] = 0;
      for (int l = 3; l < mNuNRlMax ; ++l) {
        NRprime[(l-2)*nqmax + i] = 0;
      }
      NRprime[(mNuNRlMax - 2)*nqmax + i] = 0;
    }
  }
}

static std::list<double> debugScaleFactors = std::list<double>();

void MassiveNeutrinosNew::nu2(const double a, const double betanu,
                           const double rhonu, const double pnu, const double deltaphi,
                           double *drhonu, double *fnu, double *dpnu, double *shearnu,
                           const double *psi0, const double * psi1, const double *psi2)
{
    if (!initialized)
      init();

    double q, v, g0[4], g1[nqmax+2], g2[nqmax+2], g3[nqmax+2], g4[nqmax+2], aq;
    double dum3[nqmax+2];
    int iq;
    double gf1, gf2, gf3, gf4;

#warning remove debug
    if (debugScaleFactors.empty()) {
      debugScaleFactors.push_back(1e-6);
      debugScaleFactors.push_back(1e-3);
      debugScaleFactors.push_back(30);
      debugScaleFactors.push_back(0.01);
    }

  /*
    //debug:
    // Nu_Integrate(am=  9.438214661763289E-063 , beta=   52.0000000000000      , rhonu =
    // 1.00000000000000      , pnu=  0.333333333333333      ,·
    //               dphi=  2.376670423144260E-020
    //
    if (a != 9.438214661763289E-063 ) {
      ifstream debin("debug-input.dat");
          char dummy[1024];
          double qdnDummy[16];
          debin.getline(dummy, 1024); debin.getline(dummy, 1024); debin.getline(dummy, 1024);
          for (int psiNo = -1; psiNo <= 2; ++psiNo) {
            double* psix;
            switch (psiNo) {
              case -1: psix = qdnDummy; break;
              case 0: psix = const_cast<double*>(psi0); break;
              case 1: psix = const_cast<double*>(psi1); break;
              case 2: psix = const_cast<double*>(psi2); break;
            }
            for (int i = 1; i <= 15; ++i) {
              double qin, psi;
              debin >> qin >> psi;
              psix[i-1] = psi;
            }
            debin.getline(dummy, 1024);
            debin.getline(dummy, 1024);
          }
          cout << "read psi" << endl;
            for (int i = 0; i <= 14; ++i) {
              cout << "qdnRead: " << qdnDummy[i]  << " and mine: "  << qdn[i] << endl;
            }

          for (int psiNo = 0; psiNo <= 2; ++psiNo) {
            double* psix;
            switch (psiNo) {
              case 0: psix = const_cast<double*>(psi0); break;
              case 1: psix = const_cast<double*>(psi1); break;
              case 2: psix = const_cast<double*>(psi2); break;
            }
            cout << "================= psi " << psiNo << endl;
            for (int i = 0; i <= 14; ++i) {
              cout << i << "       " << psix[i] << endl;
            }
          }

      nu2(9.438214661763289E-063 , 52., 0.01, 1.0, 10000., //1.0, 0.333333333333333, 2.376670423144260E-020, 
          drhonu, fnu, dpnu, shearnu, psi0, psi1, psi2);
    }

    //end debug
  */

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


//X     throw Bad_Error("fix massive nu");

#warning debug
    ofstream *g1o, *g2o, *g3o, *g4o, *dum3o;
    g1o=g2o=g3o=g4o=dum3o=0;
    if (false && a>debugScaleFactors.front()) {
      double scale = debugScaleFactors.front();
      debugScaleFactors.pop_front();
      std::stringstream conv;
      conv << scale;
      std::string name = "propdebug-";
      g1o = new ofstream((name+"g1at"+conv.str()+".dat").c_str());
      g2o = new ofstream((name+"g2at"+conv.str()+".dat").c_str());
      g3o = new ofstream((name+"g3at"+conv.str()+".dat").c_str());
      g4o = new ofstream((name+"g4at"+conv.str()+".dat").c_str());
      dum3o = new ofstream((name+"dum3t"+conv.str()+".dat").c_str());
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
#warning debug
      if (g1o) {
//X         (*g1o) << qdn[iq - 2] << "   " <<  (psi0[iq-2] / v/qdn[iq - 2]) << endl;
        (*g1o) << iq << "   " << (qdn[iq - 2] * psi0[iq-2] / v) << endl;
        (*g2o) << iq << "   " <<  (psi0[iq-2] * v) << endl;;
        (*g3o) << qdn[iq - 2] << "   " <<  (psi1[iq-2]/qdn[iq - 2]) << endl;;
        (*g4o) << qdn[iq - 2] << "   " <<  (psi2[iq-2]*v/qdn[iq - 2]) << endl;;
      }
    }

//X     for (iq = 1; iq <= nqmax+1; ++iq) {
//X       q = iq;
//X       aq=a/q;
//X       v = 1. / sqrt(aq * aq + 1.);
//X       double aqdn=q*q*q/(exp(q)+1.);
//X #warning check index, used to be iq only
//X #warning dbug
//X       if  (dum3o){
//X         (*dum3o) << aqdn << "   " <<  (v*v*v*a*a/(q*q)) << endl;;
//X       }
//X     }
#warning dbug
    if (g3o) {
      (*g1o).close();
      (*g2o).close();
      (*g3o).close();
      (*g4o).close();
      (*dum3o).close();
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
//X #warning deb
//X     cout << pnu << endl;
    *drhonu = (g0[0] + gf1 * 2. / (nqmax-0.5)) / 5.68219698+(rhonu-3.*pnu)*deltaphi*betanu;
    *dpnu = (g0[1] + gf2 * 2. / (nqmax-0.5)) / 5.68219698 / 3.;
//X #warning deb
    *dpnu -= tempant/3.;
//X #warning deb
    *fnu = (g0[2] + gf3 * 2. / (nqmax-0.5)) / 5.68219698;
    *shearnu = (g0[3] + gf4 * 2. / (nqmax-0.5)) / 5.68219698 * 2. / 3.;

  /*
    cout << "drhonu: " << *drhonu << endl;
    cout << "fnu:  " << *fnu << endl;
    cout << "dpnu: " << *dpnu << endl;
    cout << "shearnu: " << *shearnu << endl;
  exit(0);
   
    cout << gf1 << endl;
    cout << "-" << endl;
    cout << gf2 << endl;
    cout << "+" << endl;
    cout << gf3 << endl;
    cout <<  gf4 << endl;
    cout << g0[3] << endl;
*/

} // nu2_ 

double MassiveNeutrinosNew::nuder(const double a, const double adotoa, const double beta,
                               const double rhonu, const double phidot,
                               const double *psi2, const double *psi2dot)
{
    if (!initialized)
      init();

//X     throw Bad_Error("fix massive nu");

    static double vdot, d;

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
    for (int iq = 1; iq <= nqmax; ++iq) {
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

void MassiveNeutrinosNew::ninu1(const double a, double *rhonu, double *pnu,  const double amnu)
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

void MassiveNeutrinosNew::nu1(const double a, double *rhonu, double *pnu)
{
  double d;
  int i;

  // throw Bad_Error("nu1");

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
    *rhonu = min( exp(*rhonu) ,1.);
    *pnu = min( exp(*pnu) ,.3333333333);
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


void MassiveNeutrinosNew::setIninitalLongitudinalScalarPerturbations(double *y, const double a,
                                        const double mass_nu, const double deltan,
                                        const double Vn, const double Pi_nu)
{
  double *NR_0 = &y[firstIndex]; // massive neutrinos
  double *NR_1 = &NR_0[nqmax];
  double *NR_2 = &NR_1[nqmax];
  double *NR = &NR_2[nqmax];
#warning qspace
//X   double dq=qmax/nqmax;
  double dq = 1.;
  double aq,q,dlfdlq;
  for (int i = 0; i < nqmax; ++i) {
    q = (i+1) * dq - .5;
    aq = a * mass_nu / q;
    dlfdlq = -q / ( exp(-q) + 1.);
    NR_0[i] = -0.25*dlfdlq * deltan;
    NR_1[i] = -dlfdlq * Vn* sqrt(aq*aq+1.) / 3.;
    NR_2[i] = -dlfdlq * Pi_nu/12.;
    //   cout << " -0- " << NR_0[i] << endl;
    //      cout << " -1- " << NR_1[i] << endl;
    //      cout << " -2- " << NR_2[i] << endl;
    //      cout << "address one: " << &NR_2[i] << endl;
    for (int l = 3; l <= mNuNRlMax; ++l) {
      NR[(l-3)*nqmax + i ] = 0.;
      //   cout << "address two: " << &NR[(l-3)*nqmax + i ] << endl;
    }
  }
}
