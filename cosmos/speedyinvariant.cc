#include "speedyinvariant.h"

#include "cosmos.h"

#include "massiveneutrinos.h"


// uncomment this, if you would like to write to a file during fderivs
//X #define OUTPUTFILE
// #define OUTPUTFILETENSOR

SpeedyInvariant::SpeedyInvariant(Cosmos *cosmos)
              :Perturbation(cosmos), s5(sqrt(5.0)), s6(sqrt(6.0))
{
  /*
  cout << "at ls: " << cosmos->Visibility->integrate(cosmos->Visibility->start(), 1000) << endl;
  cout << "after ls: " << cosmos->Visibility->integrate(1000, cosmos->Visibility->stop()   ) << endl;
  throw Bad_Error("nix");
  */
  // To speed up, some fractions
  // that are needed in the Photon 
  // Distribution propagation 

  // denom1 :   l / (2l -1)
  // denom2 :   (l + 1) / (2l + 3)
  for (int l = 0; l <= LMX0; l++) {
    denom1[l] = ((double)l) / (double)(2*l -1);
    denom2[l] = ((double)l+1) / (double) (2*l +3);

    if (l > 2) {
      denomE1[l] =  sklm(2,l,0) / (double) (2*l -1);
      denomE2[l] = sklm(2,l+1,0) / (double) (2*l +3);
    }
  }

  double dilute = 5*cosmos->a_equ();
  TauDilute= cosmos->a2tau(dilute);
  //  cout << "TauDilute: " << TauDilute << endl;
}

SpeedyInvariant::~SpeedyInvariant()
{
}

void SpeedyInvariant::insertThreshold(Species species) {	
  pair<double,double> TauJump = tauTightCoupling(cosmos->Opac->start()*1.01,2000,species);
  //cout << "inserting: " << TauJump.first << " -> " << TauJump.second << " for species: " << species << endl;
  //cout << "           " << cosmos->tau2z(TauJump.first) << ", z=" << cosmos->tau2z(TauJump.second) << " for species: " << species << endl;
  bool ok=true;
  for (;;) {
    ok =true;
    for (map<double, pair<double, Species> >::iterator i = Thresholds.begin(); i!=Thresholds.end(); i++){
      if (TauJump.first >= i->first && TauJump.first <=  i->second.first) ok = false;
      if (i->first >= TauJump.first && i->first <= TauJump.second) ok =false; 
    }
    if (ok) break;
    double TauJumpGap = TauJump.second - TauJump.first;
    //  cout << "TauJumpfirst: " << TauJump.first << " TauJumpGap: " <<TauJumpGap <<endl;
    TauJump.first *= 1.00000001;
    do {
      TauJump.second = TauJump.first + TauJumpGap;
      TauJumpGap *= 1.3;
    } while (TauJump.second == TauJump.first);
  }
  Thresholds[TauJump.first] = make_pair(TauJump.second,species);
}

/*!
  This are the coefficients kappa with indices s, l and m.  Hence
  the rather strange name.
*/
double SpeedyInvariant::sklm(double s, double l, double m) {
  return sqrt( (l*l - m*m)*(l*l - s*s) / (l*l));
}

/*!
  In contrast to the other perturbation classes, isTightCoupling is 
  only used once to find the times at which we cross. These are later
  used in PropagateScalar
*/
bool SpeedyInvariant::isTightCoupling(const double tau, Species species) {
  if (tau > cosmos->tau_0())
    return false;

  double TauCompton = 1.0/opac(tau);
  double threshold = 2e-1;

  double compare=0,R=0;
  double adota = tau2adot(tau)/tau2a(tau);
  double earliestTime = cosmos->TimeStep[0];
  earliestTime -= 30;
#ifndef PRERELEASE
#warning allowing early sources
  if (earliestTime < 0) {
    earliestTime += 30;
  }
  if (warnOnEarlySource && earliestTime < 130) {
    warnOnEarlySource = false;
    cout << flush << "  Warning: using very early earliest source at tau=" << earliestTime << " . Probably a mistake." << endl << "   ";
  }
#endif
#ifdef PRERELEASE
  if (earliestTime < 100) throw Bad_Error("SpeedyInvarinat::isTightCoupling() earliest source is very early. No mistake ?");
#endif

  switch (species) {
  case photon:
    if (tau > earliestTime-2) return false; //!< Slightly different from baryon to get different times (better to code, no physics)
    if (k*TauCompton > threshold) return false;
    return true;
    break;
  case baryon:
    if (tau > earliestTime) return false;
    compare = max(k,adota);
    R = cosmos->tau2R(tau);  // 4/3 * rho_photon / rho_baryon
    if (compare*TauCompton/R > 0.2*threshold) return false;
    return true;
    break;
  case octopole:
    if (k*TauCompton > 1e-1 || TauCompton > 5e-2) return false; 
    return true;
    break;
  case multipole:
    if (k*TauCompton > 5e-1 || TauCompton > 2e-1) return false;
    return true;
    break;
  default:
    throw Bad_Error("Perturbation::isTightCoupling(): Neither photons nor baryons requested");
  }
}

/*!
  Here, the fadeout factor Gamma is calculated. 
  Essentially, it is a tanh() centerted at x=k*tau = 1000 with
  width 50.
  If you don't want to fadeout, simply return 1.0.
  This switches off the speed up due to neglectance of 
  relativistic oscillations at late time 
*/
double SpeedyInvariant::fadeOut(double tau) {
  double x = k*tau;
  double FadeOut = 1.0;
  if (x > 100) {
    double center = 1000;
    double width = 50;
    center = max(center,k*TauDilute);
    double arg = (k*tau - center)/width;   // the argument for the tanh
    // the next few lines carefully deal with the argument to prevent overflows
    double eplus=1,eminus=1;
    arg = min(arg,50.0);
    arg = max(arg,-50.0);
    eplus = exp( arg); 
    eminus = 1.0/eplus; // exp(-arg);
    FadeOut = 0.5*(  1 - (eplus - eminus)/(eplus + eminus) );
  }
  if (isnan(FadeOut)) throw Bad_Error("o mei");
  //return 1.0;
  return FadeOut;
}

/*!
  tau-derivatives of Gauge-Invariant Variables. y and yprime come
  from odeint(), meaning that they usually will not coincide with our
  objects y and yprime. However, fderivsWrapper() will at the end
  of the integration fill y and yprime of our object with the result. 
  In other words: scalarSources() etc will have a fully valid y and 
  yprime arrary and when scalarSources() calls fderivs(), it does so
  with our objects y and yprime

  This is a very high efficiency version.

*/
void SpeedyInvariant::fderivs(const double tau, const double *y, double *yprime) {
  yprime[0] = 0;
  called++;
  if (tau < 1000) called3++;

 // first some abbreviations for often used formulae

  double a= tau2a(tau);   // we take it from cosmos :-), as this is also mother, children will be fast
  double a2 = a * a;
  double adotoa = tau2adot(tau)/a;

  double Gpi4 = Gpi8()*0.5;
  double Gpi4a2 = Gpi4*a2;

  double Opac = opac(tau); // Opacity
  double TauCompton = 1.0/Opac;  //!< the compton free path length

  // Now, we calculate the suppression (if any)
  double FadeOut =1.0;
  FadeOut = fadeOut(tau);
  double K = k; // when fading out, we rescale k  for some equations
  K *= FadeOut; 
  double FadeOutOpac = FadeOut*Opac;
  double DoNothing = 1e-10; //1e-10; // if FadeOut < DoNothing the r.h.s of the photon and neutrino evol. eqn are set to 0

  // initialize E,N,M ... etc makes bookkeeping of Multipoles quite easy 
  // note:
  // y and yprime come from Miscmath::odeint() and are not the same as our
  // object's y and yprime (our object's y will however later contain 
  // the result of the odeint calculation)
  setPointers(y,yprime); 
  double Dgb = y[1];
  double Vb = y[2];
  double Dgc = y[3];
  double Vc = y[4]; 

  // Dg and V for photons
  double Dgp = 4*M[0];
  double Vp = M[1];

  double Dgn = 0;
  double Vn = 0;
  double Pin = 0;

  if (cosmos->nuR() > 0) {
    Dgn = 4*N[0];
    Vn = N[1];
    Pin = 12.0/5.0*N[2];
  }

  double w_nuNR = 0;
  rho_nuNR = P_nu_NR = delta_nu_nr_longit = onepwNRVn_nr_longit = Pi_nu_NR = 0;
  del_rho_nu_NR = f_nu_NR = del_P_nu_NR = sigma_nu_nr = 0;
  double neutrinoMass_kbT = 0;
  if (cosmos->haveMassiveNeutrinos()) {
    neutrinoMass_kbT = cosmos->mass_nu_kbT();
    w_nuNR = cosmos->tau2w_nuNR(tau);
    MassiveNeutrinos::nu1(neutrinoMass_kbT*a, &rho_nuNR, &P_nu_NR);
    MassiveNeutrinos::nu2(neutrinoMass_kbT*a, 0. /*no coupling to DE*/,
        rho_nuNR, P_nu_NR, 0. /* field perturbation, only needed when coupled to dark energy*/,
        &del_rho_nu_NR, &f_nu_NR, &del_P_nu_NR, &sigma_nu_nr, NR_0, NR_1, NR);

    delta_nu_nr_longit =del_rho_nu_NR/rho_nuNR;
    onepwNRVn_nr_longit = (1.+w_nuNR)*f_nu_NR/(rho_nuNR+P_nu_NR);
    Pi_nu_NR = 1.5*sigma_nu_nr/P_nu_NR;
  }

  // For the following, we want to be able to overwrite M2 and E2... 
  // this may not really be necessary, but still...
  double M2 = M[2], E2 = E[2];
  // if we have tight coupling and not the full quadrupole evolution yet then we estimate

  if (CoupledOctopole) { // we are overwriting the numerics (which mimics the formula anyhow) plus some corrections later on
    M2 = 8.0/9.0*k*Vp*TauCompton;
    E2 = -s6*0.25 * M2;      
    Mprime[2] = 0; // for the time being
  } else {  // if we propagate the quadrupole fully, we go ahead and get \dot M2 and \dot E2 already here
    if (FadeOut > DoNothing) {
      Mprime[2] = -FadeOutOpac*(9.0/10.0 * M[2] + s6*0.1*E[2]) + K*(2.0/3.0 * Vp - 3.0/7.0*M[3]);
      Eprime[2] = -K*( s5 / 7.0 *E[3] ) - FadeOutOpac* (E[2] + s6*0.1 * (M[2] - s6*E[2]));
    } else {
      Mprime[2] = 0;  // do nothing if fadeout has reached such a level
      Eprime[2] = 0;
    }
  }

  // first, we  calculate Phi
  double a2ak = adotoa / k;
  double dgrho = tau2rho_cdm(tau) *(Dgc + 3 *a2ak*Vc);  //we collect all species, first: cdm
  dgrho += tau2rho_b(tau) * (Dgb + 3*a2ak*Vb);              //baryons
  dgrho += tau2rho_g(tau) * (Dgp + 4*a2ak*Vp);              //photons
  dgrho += tau2rho_nu(tau) * (Dgn + 4*a2ak*Vn);           // massless neutrinos
  // massive neutrinos are in longitudinal gauge
  dgrho += tau2rho_nuNR(tau)*(delta_nu_nr_longit + 3.*a2ak*onepwNRVn_nr_longit);
  dgrho *= Gpi4a2;
  // no phicontrib for massive neutrinos, since they are in longitudinal gauge
  double phicontrib = 3*(tau2rho_b(tau) + tau2rho_cdm(tau)) + 4*(tau2rho_g(tau) + tau2rho_nu(tau));
  phicontrib *= Gpi4a2;
  Phi = dgrho / (k2 + phicontrib);

  // then we collect all the shear sources there are
  double Pip = 12.0/5.0*M2;  // photons

  double neutrinoShearContrib=0;
  neutrinoShearContrib += tau2p_nu(tau)*Pin;   // massless neutrinos
  const double p_nm = tau2p_nuNR(tau);
  neutrinoShearContrib += Pi_nu_NR*p_nm; // massive neutrinos

  // then we get Psi, first we get the total shear (times perturbation)
  double dgpres = cosmos->tau2p_g(tau)*Pip;
  dgpres *= 2.0*Gpi4a2;
  dgpres += neutrinoShearContrib * 2.*Gpi4a2;
  Psi = -Phi - dgpres / k2;

  // (rho + p)*V of all species and from this  d/dtau Phi
  double RhoPV = tau2rho_cdm(tau) * Vc + tau2rho_b(tau)*Vb;
  RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
  RhoPV += tau2rho_nuNR(tau)*onepwNRVn_nr_longit;
  PhiDot  = adotoa*Psi - Gpi4a2*(RhoPV)/k;

  double R = cosmos->tau2R(tau);  // 4/3 * rho_photon / rho_baryon
  double OneOverRPlusUnity = 1.0/(1+R);
  double cs2 = soundSpeed(tau);  // baryon soundspeed (almost 0, for tau > 1e-2) 

  // d/dtau for Dgb
  yprime[1] = -k*Vb - 3*cs2*adotoa*Dgb ; 
  // and d/dtau for M[0], which is 1/4 of Dgp 		  
  Mprime[0] = -K/3.0 * Vp;

  //now we calculate d/dtau of Dg and V for cdm
  yprime[3] = -k*Vc;
  yprime[4] =  k*Psi - adotoa * Vc;

  double slip = 0;
  // get the slip for as long as either baryons or photons can be treated as tightly coupled
  double Vb_prime=0, Vp_prime=0;  // Photon and baryon velocity prime. Init with 0, but doesnt matter
  if (CoupledBaryon || CoupledPhoton) {
    // First we need the slip just as usual 
    // the difference is that we don't use it directly but plug it into some other equation
    // the slip is d/dtau ( V_b - V_photon)  the gauge-invariant formulation of eqn 74 in Ma & Bertschinger
    // plus a Quadrupole term that is neglected there...
    // in principle we need two additional things, (d^2/dtau^2 a)/a  **and** d/dtau phi
    double adotdota = (adotoa * adotoa -  Gpi8()*a2*tau2p(tau)) * .5;

    slip = -adotdota*Vb;
    slip += - adotoa*k*(0.5*Dgp + Psi - 4.0/5.0*M2 - 2*Phi);
    slip +=  k*(cs2*yprime[1] - Mprime[0] + 2.0/5.0*Mprime[2] + PhiDot*(1-3*cs2) );
    slip *= TauCompton;  // thats the tau_compton part

    if (a < 1e-4 || true) {
      slip += adotoa*2*R*(Vb - Vp);   // main contribution 
    } else {
      slip -= (dopac(tau)/(Opac*Opac)*(1+R)  +  2)*adotoa*(Vb - Vp);
    }
    slip *= OneOverRPlusUnity; //   slip /= (1+R);
    slip /= 1.0 + 2*adotoa*TauCompton* OneOverRPlusUnity;
  }


  // get the photon evolution (minus a tine M2 correction which comes later)
  if (CoupledPhoton) {
    Vp_prime = R*k*(M[0] - 2.0/5.0*M2 - Phi) - adotoa*Vb + cs2*k*(Dgb -  3*Phi) - slip;
    Vp_prime *= OneOverRPlusUnity; //  Vp_prime /= 1+R;
    Vp_prime += k*Psi;
  } else  {
    Vp_prime = Opac*(Vb-Vp) + K*(Psi - Phi) + K*(M[0]-2.0/5.0*M[2]);  
  }


  // Here, we go to higher order for M2, if we still use the approximation formulae for M2
  // we also add a tiny correction to Vp_prime then
  double M2primeLeadingOrder=0;
  if (CoupledOctopole) {
    M2primeLeadingOrder = 8.0/9.0*k*( Vp_prime + 2*adotoa*Vp)*TauCompton;
    // the delta m2 is the leading order effect... 
    double DeltaM2 = -(11.0/6.0*M2primeLeadingOrder + k2*TauCompton*29.0/70.0*M2)*TauCompton;
    // The delta slip is really very tiny and could be omitted
    double DeltaSlip = TauCompton*OneOverRPlusUnity * k* (   2*adotoa*DeltaM2 - M2primeLeadingOrder);
    M2 += DeltaM2;
    slip += DeltaSlip;
    if (CoupledPhoton) {
      // in principle, the line below is the one to be used, however the criterion for CoupledOctopole
      // is such that the second term in the bracket is negligible, so we use the line below
      //      Vp_prime -= 2.0/5.0*k*(R + 2*adotoa*TauCompton/(1+R))*DeltaM2/(1+R); // correct in the coupled case
      Vp_prime -= 2.0/5.0*k*R*DeltaM2*OneOverRPlusUnity;
    } else Vp_prime -= 2.0/5.0*K*DeltaM2; // correct in the uncoupled case 
  }

  // evolve the baryons
  if (CoupledBaryon) {
    Vb_prime = cs2*k*(Dgb - 3*Phi) - adotoa*Vb;
    Vb_prime +=  R*(k*(-2.0/5.0*M2 +  M[0] - Phi) + slip);
    Vb_prime *= OneOverRPlusUnity;  // Vb_prime /= 1+R;
    Vb_prime += k*Psi;
  } else   Vb_prime = k*(Psi -3*cs2*Phi) + Dgb*k*cs2 - adotoa*Vb + Opac*R*(Vp-Vb);


  if (CoupledOctopole) {
    Mprime[2] = 8.0/9.0*k*( Vp_prime + 2*adotoa*Vp)*TauCompton;
    Mprime[3] =3.0/5.0*k*TauCompton*( M2primeLeadingOrder + 2*adotoa*M2);
    Eprime[2] = -s6*0.25*Mprime[2];
    Eprime[3] = k*TauCompton*(Eprime[2] + 2*adotoa*E2) / s5; 
    for (int l = 4; l <= lmaxg; l++) {  // first let's clear all from octupole on
      Mprime[l] = 0;
      Eprime[l]  = 0;
    }
  } else {
    if (CoupledMultipole) {
      Mprime[3]  = -Opac * M[3] + K*denom1[3] * M[2];
      Eprime[3]  = K*denomE1[3] * E[2]  - Opac*E[3];
      for (int l = 4; l <= lmaxg; l++) {  // first let's clear all from octupole on
        Mprime[l] = 0;
        Eprime[l]  = 0;
      }
    } else {
      if (FadeOut > DoNothing) {
        for (int l = 3; l < lmaxg; l++) Mprime[l]  = - FadeOutOpac * M[l] + K*( denom1[l] * M[l-1] - denom2[l]*M[l+1]);
        // truncation of the scheme similar to  Ma & Bertschinger 
        Mprime[lmaxg] = Combi*K*M[lmaxg - 1] - M[lmaxg]*FadeOut*( (lmaxg+1)/tau + Opac);
        for (int l = 3; l < lmaxg; l++) Eprime[l]  = K*( denomE1[l] * E[l-1] - denomE2[l]*E[l+1] ) - FadeOutOpac*E[l];
        // truncation of the scheme exactly like in Ma & Bertschinger 
        Eprime[lmaxg] =  Combi*K*E[lmaxg - 1] - E[lmaxg]*FadeOut*( (lmaxg+1)/tau + Opac);
      } else {
        for (int l = 3; l <= lmaxg; l++) {Mprime[l] =0; Eprime[l]=0;}
      }
    }
  }

  if (FadeOut < 0.5) called2++;



  Mprime[1] = Vp_prime;
  yprime[2] = Vb_prime;

  // massless neutrinos

  if (cosmos->nuR() > 0) {
    if (FadeOut > DoNothing) {
      Nprime[0] = -K/3.0 * N[1];
      Nprime[1] = K*(N[0] - 2.0/5.0*N[2] + Psi - Phi);
      for (int l = 2; l < lmaxnr; l++)
        Nprime[l] = K*(denom1[l]*N[l-1] - denom2[l]*N[l+1]);
      // truncation just like photons
      Nprime[lmaxnr] = K*((double) 2*lmaxnr +1)/((double) 2*lmaxnr-1) *N[lmaxnr - 1] - FadeOut*N[lmaxnr]* (lmaxnr+1)/tau;
    } else {
      for (int l = 0; l <= lmaxnr; l++) Nprime[l] = 0;
    }
  }

  if (cosmos->haveMassiveNeutrinos()) {
    MassiveNeutrinos::propagateMassiveNeutrinoMomentsLongitudinal(y, yprime, a, tau,
                                   neutrinoMass_kbT, k, PhiDot, Psi, 0. /* no coupling */,
                                    0. /* field perturbation, only needed for coupling */);
  }
}

void SpeedyInvariant::initialScalarPerturbations(const ControlPanel &control, const double tau)
{

  /*
  double aeq = cosmos->a_equ();
  double dilute = 5*aeq;
  double tau_dilute = cosmos->a2tau(dilute);
  */
  // cout << "k: " << k <<" tau_dilute: "<< tau_dilute << endl;

  // unfortunatly, we need to define our M,N,E ourselves, cause
  // they are non -constant (we have to initialize them)
  // see setPointers() for more
  double *M = &y[5];                    // Temperature
  double *E = &M[lmaxg-1];  // E[2] is first after M[lmaxg]
  double *N = 0;
  if (cosmos->nuR()>0) {
   N = &E[lmaxg+1];  //  massles neutrinos
  }

  //  cout << "N[lmaxnr] " << &N[lmaxnr] << "  y[nvar]: " << &y[nvar] << endl;
  double Dgp, Dgc, Dgb, Dgn, Vb, Vc, Vp, Vn, Psi, Pi_nu;
  k2 = k*k;

  double x = k*tau;

  double Omega_nu = cosmos->rho2omega(cosmos->tau2rho_nu(tau),tau) + cosmos->rho2omega(cosmos->tau2rho_nuNR(tau), tau);
  double Omega_gamma = cosmos->rho2omega(cosmos->tau2rho_g(tau),tau);
  double Omega_cdm = cosmos->rho2omega(cosmos->tau2rho_cdm(tau),tau);
  double Omega_b = cosmos->rho2omega(cosmos->tau2rho_b(tau),tau);

  double Q = 1.0/(3.0*Omega_cdm+8.0*(Omega_nu + Omega_gamma));

  double P = 1.0/(4.0*Omega_nu + 15.0);
  double T=1./4.;
  double U =1.0/(30.0+4.0*Omega_nu);

  double ad=control.adiabaticContribution; // adiabatic contribution to mixed
  double icdm=control.isoCDMContribution; // isoCDM contribution to mixed
  double ibaryon=control.isoBaryonContribution; //isoBaryon contribution to mixed
  double ineutrino=control.isoNeutrinoContribution; //isoNeutrino contribution to mixed

  double normalize = 1.0; // see comment below on rescaling

  switch (control.initialCond) {
  case ControlPanel::adiabatic: 
    Dgp = 1;
    Dgn = Dgp;
    Dgb = 3.0/4.0*Dgp;
    Dgc = Dgb;
    Vp = -5.0/4.0*P*x*Dgp;
    Vb = Vp;
    Vn = Vp;
    Vc = -5.0/4.0*P*x*Dgp;
    Pi_nu = -Dgp*P*x*x;
    Phi = (5.0+2.0*Omega_nu)/(30.0+8.0*Omega_nu)*Dgp;
    Psi= -Phi-Omega_nu*Pi_nu/(x*x);
    normalize = Psi * (4.0*Omega_nu + 15.0)/10;
    normalize= -1/normalize;
    break;
  case ControlPanel::isoCDM:
    Dgc = 1.0;
    Dgp = 0*Dgc;
    Dgn = Dgp;
    Dgb = 3.0/4.0*Dgp;

    Vc = (4.0*Omega_nu-15.0)*Omega_cdm*U*x*Dgc/12.0;
    Vp = -15.0/4.0 * Omega_cdm*U*x;
    Vn = Vp;
    Vb = Vp;
    Pi_nu = -2*Omega_cdm*U*x*x;

    Phi = 2.0*Omega_cdm*(5.0+4.0*Omega_nu)*P*Q*Dgc;  
    Psi= -Phi-Omega_nu*Gpi8()*Pi_nu/(x*x);
   break;
 case ControlPanel::isoBaryon:
    //first order
    Dgc = 0;
    Dgp = 0;
    Dgn = 0;
    Dgb = 1;
    Vc = Omega_b*(4*Omega_nu-15)/3*T*U*Dgb*x;
    Vp = -15*Omega_b*T*U*Dgb*x;
    Vn = Vp;
    Vb = Vp;
    Pi_nu = -8*Omega_b*T*U*Dgb*x*x;
    Phi = Omega_b*(15+4*Omega_nu)*U/4*Dgb;  
    Psi= -Phi-Omega_nu*Gpi8()*Pi_nu/(x*x);
    break;
 case ControlPanel::isoNeutrino:
    Dgp= 1.0;
    Dgn = -Omega_gamma/Omega_nu*Dgp;
    Dgc = 3.0/4.0*Dgp;
    Dgb = 3.0/4.0*Dgp;
    Vc = Omega_gamma*P*x*Dgp;
    Vp = (15.0+4.0*Omega_gamma+4.0*Omega_nu)/4.0*P*x*Dgp;
    Vn = -15.0/4.0*Omega_gamma/Omega_nu*P*x*Dgp;
    Vb = Vp;
    Pi_nu = - 3.0*Omega_gamma/Omega_nu*P*x*x*Dgp;
    Phi = Omega_gamma/(15.0+4.0*Omega_nu)*Dgp;  
    Psi= -Phi-Omega_nu*Gpi8()*Pi_nu/(x*x);

    break;
 case ControlPanel::mixed:
    // linear combination of the four fundamental modes
    Dgp = 1*ad + 0*icdm +0*ibaryon+ 1*ineutrino;
    Dgn = 1*ad + 0*icdm +0*ibaryon -Omega_gamma/Omega_nu*ineutrino ;
    Dgb = 3.0/4.0*ad + 0*icdm +1*ibaryon + 3.0/4.0*ineutrino ;
    Dgc = 3.0/4.0*ad +  1*icdm + 0*ibaryon+ 3.0/4.0*ineutrino;
    Vp = -5.0/4.0*P*x*ad  -15*Omega_cdm*T*U*x*icdm -15*Omega_b*T*U*x*ibaryon + (15.0+4.0*Omega_gamma+4.0*Omega_nu)/4.0*P*x*ineutrino ;
    Vb = -5.0/4.0*P*x*ad  -15.0*Omega_cdm*T*U*x*icdm -15.0*Omega_b*T*U*x*ibaryon + (15.0+4.0*Omega_gamma+4.0*Omega_nu)/4.0*P*x*ineutrino;
    Vn = -5.0/4.0*P*x*ad  -15.0*Omega_cdm*T*U*x*icdm -15.0*Omega_b*T*U*x*ibaryon  -15.0/4.0*Omega_gamma/Omega_nu*P*x*ineutrino ;
    Vc = -5.0/4.0*P*x*ad + Omega_cdm*(4*Omega_nu-15)/3*T*U*x*icdm+  Omega_b*(4*Omega_nu-15)/3*T*U*x*ibaryon + Omega_gamma*P*x*ineutrino ;
    Pi_nu = -ad*P*x*x -Omega_cdm*8*T*U*x*x*icdm -Omega_b*8*T*U*x*x*ibaryon + 3.0*Omega_gamma/Omega_nu*P*x*x*ineutrino ;
    Phi = (5.0+2.0*Omega_nu)/(30.0+8.0*Omega_nu)*ad + Omega_cdm*(15.0+4.0*Omega_nu)*T*U*icdm  + Omega_b*(15.0+4.0*Omega_nu)*T*U*ibaryon+ Omega_gamma/(15.0+4.0*Omega_nu)*ineutrino ; 
    Psi= -Phi-Omega_nu*Gpi8()*Pi_nu/(x*x);
    break;
  default:
    throw Bad_Error("SpeedyInvariant::initialScalarPerturbations():  This initial condition is notsupported.");
  }

  // Normalize all variables such that the curvature chi is =-1
  // Needed for correct initial Powerspectrum Amplitude
  Dgp *= normalize;
  Dgn *= normalize;
  Dgb *= normalize;
  Dgc *= normalize;
  Vp *= normalize;
  Vb *= normalize;
  Vn *= normalize;
  Vc *= normalize;
  Pi_nu *= normalize;
  Phi *= normalize;
  Psi *= normalize;
  // Hm. Well. Done. Sorry for the inconvenience caused ...

  // if we want to crudely change one of the perturbation variables independently
  // from what adiabaticity or isocurvature tell us, we can do so using InitialConditionFactors
  // in the ControlPanel (so it can be used e.g. from the montecarlo driver/xdriver etc.)
  Dgp *= control.initialConditionFactors.DgpFactor;
  Dgn *= control.initialConditionFactors.DgnFactor;
  Dgb *= control.initialConditionFactors.DgbFactor;
  Dgc *= control.initialConditionFactors.DgcFactor;
  Vp *= control.initialConditionFactors.VpFactor;
  Vb *= control.initialConditionFactors.VpFactor;
  Vn *= control.initialConditionFactors.VnFactor;

  //  cout << "Psi: "  << Psi << "  Psi*norm: " << (4*Omega_nu + 15)*Psi/10*normalize << endl;

  // fill matter
  y[1] = Dgb;
  y[2] = Vb;
  y[3] = Dgc;
  y[4] = Vc;

  // fill photons and neutrinos
  M[0] = Dgp / 4.0;
  M[1] = Vp;
  // set all higher moments to 0, obs: E[2] is the first one needed
  for (int l = 2; l <= lmaxg; l++) {
    M[l] =0; E[l] = 0;
  }
  if (cosmos->nuR()>0) {
    N[0] = Dgn / 4.0;
    N[1] = Vn;
    N[2] = 5.0/12.0*Pi_nu;
    for (int l = 3; l <= lmaxnr; l++) {
      N[l] =0;
    }
  }


  if (cosmos->haveMassiveNeutrinos() && (control.initialCond != ControlPanel::adiabatic)) {
    throw Bad_Error("For massive neutrinos, only adiabatic initial conditions are implemented.");
  } else if (cosmos->haveMassiveNeutrinos()) {
    const double mass_nu = cosmos->mass_nu_kbT();
    double deltan = Dgn - 3.*(1.+cosmos->tau2w_nuNR(tau))*Phi;
    if (cosmos->nuR()>0) {
      MassiveNeutrinos::setFirstIndex((N-y)+lmaxnr+1);
    } else {
      MassiveNeutrinos::setFirstIndex((E-y)+lmaxg+1);
    }
    MassiveNeutrinos::setIninitalLongitudinalScalarPerturbations(y, cosmos->tau2a(tau), mass_nu, deltan, Vn, Pi_nu);
  }
}

/*
void SpeedyInvariant::initialTensorPerturbations() {	

#ifdef OUTPUTFILETENSOR   
  ofs = new ofstream("colddark.dat");
  ofs2 = new ofstream("colddark2.dat");
#endif

  yt[1] = 1;
  yt[2] = 0.;
  int ind1 = 3;
  int ind2 = ind1 + lmaxt +1;
  int ind3 = ind2 + lmaxt +1;
  for (int l = 0; l <= lmaxt; ++l) {
    yt[ind1 + l] = 0.;
    yt[ind2 + l] = 0.;
  }
  for (int l = 0; l <= lmaxnr; l++) {
    yt[ind3 + l] = 0;
  }
}
*/


#define LM2 7 // 7
#define LM3 4
#define NQ1 15

void SpeedyInvariant::getReady(const ControlPanel& control) {
  lasttau = 0;
  called=0; called2=0;called3=0;
  if (control.scalar) {
    lmaxg = LMAX0;          // number of photon momenta
    //    if (k > 0.5*0.7) lmaxg = 2500;
    if (control.highPrecisionTransfer)
      lmaxnr = LMAXNR0; else lmaxnr = LM2;
  }


  Combi =  ((double) 2*lmaxg +1)/((double) 2*lmaxg-1);

  //  cout << " k:  " << k << "  got ready: " << lmaxg << "  "  << lmaxnr << endl;

  if (control.tensor)
    lmaxt = LMAXT0; else lmaxt = 0;

  insertThreshold(photon);
  insertThreshold(baryon);
  insertThreshold(octopole);
  insertThreshold(multipole);

  CoupledPhoton = true;
  CoupledBaryon = true;
  CoupledOctopole = true;
  CoupledMultipole = true;

  Perturbation::getReady(control);  // now call parent objects get ready
}

/*!
  For the sources, to be precise the ISW effect, one needs d/d tau
  of Phi and Psi. In the QuintSpeedyInvariant class, these are needed anyhow
  within the fderivs() function. However, here they are only calculated
  for the sources, which is why we need to do it.

  Please note that fillPhiDot() assumes that it is called after fderivsWrapper()
  has been called, i.e. Phi and Psi and all other stuff is up to date.
*/
void SpeedyInvariant::fillPhiDot(double tau)
{
  double Vb = y[2];
  double Vc = y[4];
  double Vp = M[1];
  double Vn = 0;
  if (cosmos->nuR()>0) {
    Vn = N[1];
  }
  double a = tau2a(tau);   // we take it from cosmos :-), as this is also mother, children will be fast
  double a2 = a * a;
  double adotoa = tau2adot(tau)/a;

  double Gpi4 = Gpi8()*0.5;
  double Gpi4a2 = Gpi4*a2;

  double Pip = 12.0/5.0*M[2];  // photons

  double RhoPV = tau2rho_cdm(tau) * Vc + tau2rho_b(tau)*Vb;
  RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
  RhoPV += tau2rho_nuNR(tau)*onepwNRVn_nr_longit;
  PhiDot  = adotoa*Psi - Gpi4a2*(RhoPV)/k;

  double Pin = 0;
  if (cosmos->nuR()>0) {
    Pin = 12.0/5.0*N[2];   // massless neutrinos
  }

  double Psi2nd = 12.0/5.0*a2*(tau2rhodot_g(tau)*M[2] + tau2rho_g(tau)*Mprime[2])/3.0;
  if (cosmos->nuR()>0) {
    Psi2nd     += 12.0/5.0*a2*(tau2rhodot_nu(tau)*N[2] + tau2rho_nu(tau)*Nprime[2])/3.0;
    Psi2nd += 2.0*a*tau2adot(tau)*tau2p_nu(tau)*Pin;
  }
  Psi2nd   += 2.0*a*tau2adot(tau)*tau2p_g(tau)*Pip;

  if (cosmos->haveMassiveNeutrinos()) {
    double shearnudot = MassiveNeutrinos::nuder(a*cosmos->mass_nu_kbT(), adotoa,
                                                0 /*no coupling*/, tau2rho_nuNR(tau),
                                                0 /*field pert, only needed for coupling*/,
                                                NR, NRprime);
    double w_nuNR = cosmos->tau2w_nuNR(tau);
    double wdot_nuNR = cosmos->Tau2wdot_nuNR->fastY(tau);
    double pNuNR = tau2p_nuNR(tau);
    double rhodot_nuNR = cosmos->Tau2RhoDot_nuNR->fastY(tau);
    double pdotNuNR = tau2rho_nuNR(tau)*wdot_nuNR+rhodot_nuNR*w_nuNR;
    double PiDot_nu_NR = 1.5/P_nu_NR*shearnudot;
    PiDot_nu_NR += -Pi_nu_NR*(pdotNuNR/pNuNR+4.*adotoa);
    Psi2nd += 2.*tau2adot(tau)*a*pNuNR*Pi_nu_NR + a2*(pdotNuNR*Pi_nu_NR+pNuNR*PiDot_nu_NR);
  }

  PsiDot = -PhiDot - Gpi8()/k2 * Psi2nd;
}


void SpeedyInvariant::scalarSources(double tau, double *d, double *dp, double *dk) {  
  fderivsWrapper(tau,y,yprime);  // get Phi and all derivatives of M[l] etc
  fillPhiDot(tau);   // get PhiDot and PsiDot right

  double Dgp = 4*M[0];
  double Vb = y[2];

  double FadeOut = fadeOut(tau);
  Dgp = Dgp*FadeOut - (1-FadeOut)*8*Psi;  // Dgp = -8*Psi is the particular solution after fadeout

  double C = (M[2] - s6*E[2])/10.0; 
  double Cdot = (Mprime[2] - s6*Eprime[2])/10.0;
  double M2ddot = -dopac(tau)*(9*M[2] + s6*E[2])/10.0 - opac(tau)*(9*Mprime[2] + s6*Eprime[2])/10.0;
  M2ddot += k*(2.0/3.0*Mprime[1] - 3.0/7.0*Mprime[3]);
  double E2ddot = -k*s5/7.0 *Eprime[3] - dopac(tau)*(E[2] + s6*C) - opac(tau)*(Eprime[2] + s6*Cdot);
  double Cddot =  (M2ddot - s6*E2ddot)/10.0;
  double Vbdot = yprime[2];

  C = C*FadeOut + (1-FadeOut)*0;
  Cdot = Cdot*FadeOut + (1-FadeOut)*0;
  Cddot = Cddot*FadeOut + (1-FadeOut)*0;

  *d =  -expmmu(tau)*(PhiDot - PsiDot);  // ISW
  *d += visibility(tau)*( -Phi + Psi  + 0.25*Dgp + Vbdot/k + C/2.0 + 3.0/(2.0*k2)*Cddot );
  *d += dvisibility(tau)*(Vb / k + 3.0*Cdot/(k2));
  *d += ddvisibility(tau)*3.0/(2.0*k2)*C;

  double x = k * (tau_0() - tau);
  if (x > 0.) *dp = visibility(tau) * 3. / 2. * C / (x * x); else *dp = 0.;
}

/*
void SpeedyInvariant::fderivsTensor(const double tau,   double *y, double *yprime) {
    static double  psie, htpr, a;
    static int l;
    static double a2, htdpr, rhonu, ep, ht, adotoa;

    static double cs2;
    static double deltap0, deltat0, tcp, pnu;
    static int ind1, ind2,ind3;
    static double tcp1, tcp2;

    //  Evaluate the time derivatives of the perturbations. 
    // ep is used to stop the tight coupling approximation. 
    if (k > epsw * .06) ep = .01; else ep = .0117;


    a = tau2a(tau);
    a2 = a * a;
    cs2 =soundSpeed(tau);

    // Tight Coupling parameters 
    tcp = 0.;
    tcp1 = k / opac(tau);
    tcp2 = 1. / (opac(tau) * tau);

    if (tcp1 > ep || tcp2 > ep) tcp = 1.;

    //  Compute expansion rate. 
    if (cosmos->amnu == 0.) {
      rhonu = 1.;
      pnu = .33333333333333331;
    } else cosmos->nu1(a, &rhonu, &pnu);

//  8*pi*G*rho*a**2 and 8*pi*G*P*a**2. 
// Tensors 
ind1 = 3;
ind2 = ind1 + lmaxt + 1;
    ind3 = ind2 + lmaxt + 1;
    double shearg=y[ind1]/15+y[ind1+2]/21+y[ind1+4]/35;
    double shearr=y[ind3]/15+y[ind3+2]/21+y[ind3+4]/35;

    double Pi=  Gpi8()*a2*(tau2p_g(tau)*shearg+tau2p_nu(tau)*shearr);

    adotoa = tau2adot(tau)/a;
    ht = y[1];
    htpr = y[2];
    yprime[1] = htpr;
    htdpr = adotoa * -2 * htpr - k2 * ht + 24*Pi;
    yprime[2] = htdpr;

    // Photon perturbations 
    psie = y[ind1] / 10. + y[ind1 + 2] / 7. + y[ind1 + 4] * 3. / 70. - y[ind2]
                 * 3. / 5. + y[ind2 + 2] * 6. / 7. - y[ind2 + 4] * 3. / 70.;

    if (tcp == 1.) {
      // no tight coupling approx 
	yprime[ind1] = -k * y[ind1 + 1] - opac(tau) * y[ind1] + opac(tau) * psie - htpr;
 	yprime[ind2] = -k * y[ind2 + 1] - opac(tau) * y[ind2] - opac(tau) * psie;
	// l=1...lmaxt 
	for (l = 1; l < lmaxt ; ++l) {       
	    yprime[ind1 + l] = k * denl[l] * (l * y[
		    ind1 - 1 + l] - (l + 1) * y[ind1 + 1 + l]) - opac(tau) * y[ind1 + l];
	    yprime[ind2 + l] = k * denl[l] * (l * y[
		    ind2 - 1 + l] - (l + 1) * y[ind2 + 1 + l]) - opac(tau) * y[ind2 + l];
	}

	// Truncate moment expansion  
	yprime[ind1 + lmaxt] = k * y[ind1 - 1 + 
		lmaxt] - (lmaxt + 1) / tau * y[ind1 + 
		lmaxt] - opac(tau) * y[ind1 + lmaxt];
	yprime[ind2 + lmaxt] = k * y[ind2 - 1 + 
		lmaxt] - (lmaxt + 1) / tau * y[ind2 + 
		lmaxt] - opac(tau) * y[ind2 + lmaxt];
    } else {
      deltat0 = htpr * -4. / opac(tau) / 3.;
      deltap0 = -deltat0 / 4.;
      y[ind1] = deltat0;
      y[ind2] = deltap0;
      for (l = 0; l <= lmaxt ; ++l) {
	    yprime[ind1 + l] = 0.;
	    yprime[ind2 + l] = 0.;   
      }
    }

    // massless neutrinos
    yprime[ind3] = -k*y[ind3+1]-htpr;
    for (l= 1; l < lmaxnr; l++) {  // l=1,lmaxnr-1
      yprime[ind3+l]=k*denl[l]*(l*y[ind3-1+l]-(l+1)*y[ind3+1+l]);
    }
    // Truncate moment expansion
    yprime[ind3+lmaxnr]=k*y[ind3-1+lmaxnr]-(lmaxnr+1)/tau     *y[ind3+lmaxnr];
}  

void SpeedyInvariant::tensorSources(double tau, double *dt, double *dte, double *dtb) {
  double psie, htpr, psieddot, x, htdpr, x2, psiedot;
  int ind1, ind2;
  
  fderivsTensor(tau,yt,ytprime);
  double *ypr = ytprime;  // shorthand 
  
  
  x = k * (tau_0() - tau);
  x2 = x * x;
  ind1 = 3;
  ind2 = ind1 + lmaxt + 1;
  
  htpr = yt[2];
  htdpr = ypr[2];
  psie = yt[ind1] / 10. + yt[ind1 + 2] / 7. + yt[ind1 + 4] * 3. / 70. - yt[ind2]
    * 3. / 5. + yt[ind2 + 2] * 6. / 7. - yt[ind2 + 4] * 3. / 70.;
  psiedot = ypr[ind1] / 10. + ypr[ind1 + 2] / 7. + ypr[ind1 + 4] * 3. / 70. 
    - ypr[ind2] * 3. / 5. + ypr[ind2 + 2] * 6. / 7. - ypr[ind2 + 4] * 3. / 70.;
  psieddot = (opac(tau) * psiedot + dopac(tau) * psie)
    * -.3 - htdpr * .1 - k * (ypr[ind1 + 1] * 3. / 70. 
	    + ypr[ind1 + 3] / 15. + ypr[ind1 + 5] / 42. - ypr[ind2 + 1] * 33. 
			      / 35. + ypr[ind2 + 3] * 8. / 15. - ypr[ind2 + 5] / 42.);
  
  if (x > 0.) {
    //    if (tau > 200) htpr = 0;
    *dt =  (- expmmu(tau) * htpr + visibility(tau)* psie) / x2;  

#ifdef OUTPUTFILETENSOR
    if (k > 0.173 && k < 0.178) {
      cout << ":OUTT" <<   "  " << ofs << endl;
      (*ofs) << tau << "  "<< *dt << endl;
     
    } 
#endif 

    *dte = visibility(tau) * (psie - psieddot / k2 - 
			      psie * 6. / x2 - psiedot * 4. / k / x) - 
      dvisibility(tau)* (psie * 4. / x / k + psiedot * 2. / k2) 
      - ddvisibility(tau) * psie / k2;
    *dtb = (visibility(tau) * (psie * 2. / x + psiedot / k) 
	    + dvisibility(tau) * psie / k) * 2.;
  } else {
      *dt = 0.;
      *dte = 0.;
      *dtb = 0.;
  }    
  *dte = -(*dte);
  *dtb = -(*dtb);  
} 
*/


void SpeedyInvariant::calcPerturbationNr(const ControlPanel &control) {
  //   Calculate number of equations
  if (control.scalar) {
    nvar = 4;  // CDM and Baryon Delta and Velocity
    nvar += lmaxg +1; // Photons  0.. lmaxg
    nvar += lmaxg -1; // Polarization 2 .. lmaxg
    if (cosmos->nuR() > 0) {
      nvar += lmaxnr +1; // massless neutrinos
    }
    if (cosmos->haveMassiveNeutrinos()) {
      nvar += MassiveNeutrinos::pertArraySize();
    }
    //cout << "Allocated " << nvar << " entries." << endl;
  } else nvar = 0;

  if (control.tensor) nvart = 2*lmaxt+4+lmaxnr+1+lmaxnu+1; else nvart = 0;

  y = new double[nvar+1];
  yprime = new double[nvar+1];
  yt = new double[nvart+1];
  ytprime = new double[nvart+1];

  std::fill(y, y+nvar+1, 0.);
  std::fill(yprime, yprime+nvar+1, 0.);
  std::fill(yt, yt+nvart+1, 0.);
  std::fill(ytprime, ytprime+nvart+1, 0.);
}

/*!
  Set M, E, N, Nm, Mprime, Eprime, Nprime, Nmprime such that
  they have the right position within y and yprime.
  This just centralises a bit the bookkeeping for the
  fderivs() and initialScalarPerturbation()
*/
void SpeedyInvariant::setPointers(const double *y, double *yprime) {
  // set pointers M, E, N for Temperature, Polarization and 
  // massless neutrinos. Just convenience, but diff eqn. much more
  // readable with this
  M = &y[5];                    // Temperature
  // Polarization, E[2] is the first one we need and we 
  // therefore readjust the beginning such that E[2] points to the first
  // free variable after the M's
  E = &M[lmaxg-1];  // E[2] is first after M[lmaxg]

  Mprime = &yprime[5];
  Eprime = &Mprime[lmaxg-1];

  N = Nm = 0;
  Nprime = Nmprime = 0;
  if (cosmos->nuR()>0) {
    N = &E[lmaxg+1]; // N[0] is first after E[lmaxg]
    Nprime = &Eprime[lmaxg+1];
    if (cosmos->haveMassiveNeutrinos()) {
      Nm = &N[lmaxnr+1];  //  massive neutrino momenta
      Nmprime = &Nprime[lmaxnr+1];
    }
  } else if (cosmos->haveMassiveNeutrinos()) {
      Nm = &E[lmaxg+1];
      Nmprime = &Eprime[lmaxg+1];
  }

  int temp_nqmax = MassiveNeutrinos::qGridSize();
  NR_0 = NR_1 = NR = NRprime_0 = NRprime_1 = NRprime = 0;

  if (cosmos->haveMassiveNeutrinos()) {
    NR_0 = Nm;
    NR_1 = &NR_0[temp_nqmax];
    NR = &NR_1[temp_nqmax];
    NRprime_0 = Nmprime;
    NRprime_1 = &NRprime_0[temp_nqmax];
    NRprime = &NRprime_1[temp_nqmax];
  }
}


void SpeedyInvariant::propagateScalar(double *tau, const double tauend, const double precision)
{
  bool original = false;
  k2 = k*k;
  if (hnext == 0) hnext = (*tau)*1e-1;

  /*
     for (map<double, pair<double, Species> >::iterator i = Thresholds.begin(); i!=Thresholds.end(); i++) {
     cout << "Jump at: " << i->first << " width: " << i->second.first - i->first <<  " species: " << i->second.second << endl;
     }
     */ 

  do {
    bool YourTurn = true;
    if (!Thresholds.empty()) {
      if (Thresholds.begin()->first < tauend) {
        YourTurn = false;
        double Jump = Thresholds.begin()->first; // until the jump
        // cout << "jump propagating: from " << *tau << " to: " << Jump << endl;
        if (Jump > *tau) {
          hnext = Miscmath::odeint(y,nvar,*tau,Jump,precision,hnext,0, (moDerivs)&Perturbation::fderivsWrapper, *this,original);
          *tau = Thresholds.begin()->second.first;  // set it to right after the jump
        }
        double Vp=0,Vp_prime=0,adotoa=0,M2=0,M2primeLeadingOrder=0; // init for switch (c++ thing)
        double* yprime=0 ; 
        switch (Thresholds.begin()->second.second) {  // which species
          case photon:
            CoupledPhoton = false;
            break;
          case baryon:
            CoupledBaryon = false;
            break;
          case octopole:
            CoupledOctopole = false;
            hnext *= 1e-2;
            yprime =  new double[nvar+1];
            if (yprime) {
              setPointers(y,yprime); 
              fderivs(*tau,y,yprime);
              Vp = M[1];
              Vp_prime = Mprime[1];
              adotoa =   tau2adot(*tau)/tau2a(*tau);
              M2 = 8.0/9.0*k*Vp/opac(*tau);
              //cout << "M2 leading: " << M2 ;
              M2primeLeadingOrder = 8.0/9.0*k*( Vp_prime + 2*adotoa*Vp)/opac(*tau);
              M2 -= (11.0/6.0*M2primeLeadingOrder + k2/opac(*tau)*29.0/70.0*M2)/opac(*tau);
              //cout << " now: " << M2 << endl;
              y[7] = M2;
              delete[] yprime;
            }
            break;
          case multipole:
            CoupledMultipole = false;
            break;
          default:
            throw Bad_Error("SpeedyInvariant::propagateScalar() no such species");
        }
        Thresholds.erase(Jump);
      }
    }
    if (YourTurn) {
      // cout << "ordinary propagation: "  << *tau << " to: " << tauend<< endl;
      hnext = Miscmath::odeint(y,nvar,*tau,tauend,precision,hnext,0, (moDerivs)&Perturbation::fderivsWrapper, *this,original);
      *tau = tauend;
    }
  } while (*tau < tauend);

  for (int i=1; i<=nvar; i++) {
    if (isnan(y[i])) {
      cout << endl << "y["<<i<<"] isnan at time: " << *tau << "  [ " << cosmos->tau2a(*tau) << " ]"<<endl;
      throw Bad_Error("isnan !");
    }
  }
}
