#include "speedydeinvariant.h"
#include "quintcosmos.h"

SpeedyDEInvariant::SpeedyDEInvariant(QuintCosmos *cosmos) : SpeedyInvariant(cosmos), quintcosmos(cosmos) {
}

static bool noDarkEnergyPerturbations = true;

/*!
  tau-derivatives of Gauge-SpeedyDEInvariant Variables. y and yprime come
  from odeint(), meaning that they usually will not coincide with our
  objects y and yprime. However, fderivsWrapper() will at the end
  of the integration fill y and yprime of our object with the result. 
  In other words: scalarSources() etc will have a fully valid y and 
  yprime arrary and when scalarSources() calls fderivs(), it does so
  with our objects y and yprime
  
  This is a very high efficiency version.

*/
void SpeedyDEInvariant::fderivs(const double tau, const double *y, double *yprime) {
  //cout << "fderivs " << k << " " << tau2a(tau) << endl;
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

  // dark energy
  double Dq = y[qidx]; // please note that this is in longitudinal gauge, not gauge invariant like the rest
  double Vq = y[qidx+1]; // gauge invariant = longitud for V, however we have (1+w)*V here (better conditioned)

  //double X = y[xidx];
  //double Xprime = y[xidx+1];

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

    delta_nu_nr_longit = del_rho_nu_NR/rho_nuNR;
    onepwNRVn_nr_longit = (1.+w_nuNR)*f_nu_NR/(rho_nuNR+P_nu_NR);
    Pi_nu_NR =1.5*sigma_nu_nr/P_nu_NR;
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
  // the dark energy  fluid has a special formula, because we propagate not D_g, but delta_longit
  // please do also note that V_dark has the (1+w) built in already
  dgrho += quintcosmos->tau2rho_q(tau)*(Dq + 3*a2ak*Vq);
  // massive neutrinos are in longitudinal gauge, too
  dgrho += tau2rho_nuNR(tau)*(delta_nu_nr_longit + 3.*a2ak*onepwNRVn_nr_longit);
  dgrho *= Gpi4a2;  // note: Gpi4a2 = a^2/(2 Mpl^2)
  // there is no Phi Contrib from the DE fluid nor from massive nu, because we did not use D_g, but delta_longit...
  // hence no need to subtract the Phi that came in to make delta gauge invariant
  double phicontrib = 3*(tau2rho_b(tau) + tau2rho_cdm(tau)) + 4*(tau2rho_g(tau) + tau2rho_nu(tau));
  phicontrib *= Gpi4a2;
  Phi = dgrho / (k*k + phicontrib);

  //const double qdot = quintcosmos->Tau2PhiDot->fastY(tau);
  //double wq = quintcosmos->tau2w_q(tau);

  // then we collect all the shear sources there are
  double Pip = 12.0/5.0*M2;  // photons

  double neutrinoShearContrib=0;
  if (cosmos->nuR()>0) {
    neutrinoShearContrib += cosmos->tau2p_nu(tau)*Pin;   // massless neutrinos
  }
  if (cosmos->haveMassiveNeutrinos()) {
    const double p_nm = cosmos->tau2p_nuNR(tau);
    neutrinoShearContrib += Pi_nu_NR*p_nm; // massive neutrinos
  }

  // then we get Psi, first we get the total shear (times perturbation)
  double dgpres = cosmos->tau2p_g(tau)*Pip;
  dgpres += neutrinoShearContrib;
  dgpres *= 2.0*Gpi4a2;
  Psi = -Phi - dgpres / k2;

  double RhoPV = tau2rho_cdm(tau) * Vc + tau2rho_b(tau)*Vb;
  RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
  RhoPV += tau2rho_nuNR(tau)*onepwNRVn_nr_longit;
  RhoPV += quintcosmos->tau2rho_q(tau)* Vq; // simple :-)
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
                                    0. /* field pert, only needed for coupling*/);
  }

  // dark energy

  double w = quintcosmos->tau2w_q(tau);
  double wdot =  quintcosmos->tau2wdot_q(tau);

  double mu2 =  quintcosmos->quintessence()->speedOfSound2(); // speed of sound in DE rest frame
  double P_Anisotropic = 0; // no anisotropic pressure for this dark energy
  double P = 0; // Pressure
  if (IsPhantomCrossing && fabs(1+w) < 1e-2) 
    P = mu2*Dq; // cs^2 always for crossing dark energy models
  else  P =  mu2*Dq + Vq/k*( 3.0*adotoa*(mu2 - w) +  wdot/(1+w));


  double Psi2nd = 12.0/5.0*a2*(tau2rhodot_g(tau)*M[2]
                  + tau2rho_g(tau)*Mprime[2])/3.0;
  if (cosmos->nuR()>0) {
    Psi2nd += 12.0/5.0*a2*(tau2rhodot_nu(tau)*N[2] + tau2rho_nu(tau)*Nprime[2])/3.0;
  }
  Psi2nd += 2.0*a*tau2adot(tau)*(tau2p_g(tau)*Pip + tau2p_nu(tau)*Pin);

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

  yprime[qidx] = -3*(1+w)*PhiDot  -  k*Vq - 3*adotoa*(P - w*Dq);
  yprime[qidx+1] = (3*w - 1)*adotoa*Vq + k*( (1+w)*Psi + P  -  2.0/3.0 * P_Anisotropic);

  // CAREFUL with FADEOUT for models with mu2  closer to 0 !!! Only allowed
  // if perturbations decay relativistically within the horizon

  if (mu2 > 1e-2) {
    yprime[qidx] *= FadeOut;
    yprime[qidx+1] *= FadeOut;
  }


#if 0
   Quintessence* quint = quintcosmos->quintessence();
   double q = quintcosmos->tau2phi(tau);
   double Xdoubleprime = -2.*adotoa*Xprime - (k*k + a2 * quint->Vprime2(q, a, adotoa, tau))*X
                        + (PsiDot - 3.*PhiDot)*qdot - 2.*a2*quint->Vprime(q, a, adotoa) * Psi;
   yprime[xidx] = Xprime;
   yprime[xidx+1] = Xdoubleprime;
#else
   yprime[xidx] = 0;
   yprime[xidx+1] = 0;
#endif


#ifdef OUTPUTFILE
  (*ofs) << tau << " " << Dq << " " << Vq << " " << endl;
#endif
}

void SpeedyDEInvariant::initialScalarPerturbations(const ControlPanel &control, const double tau)
{
  // copy the phantomcrossing flag for fderivs()
  IsPhantomCrossing = control.isPhantomCrossing();

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
  double Dgp, Dgc, Dgb, Dgn, Dgq;
  double Vb,   Vc,  Vp,  Vn,  Vq;
  double Pi_nu;
  k2 = k*k;

  double x = k*tau;

  double Omega_nu = cosmos->rho2omega(cosmos->tau2rho_nu(tau),tau) + cosmos->rho2omega(cosmos->tau2rho_nuNR(tau), tau);
  double Omega_gamma = cosmos->rho2omega(cosmos->tau2rho_g(tau),tau);
  double Omega_cdm = cosmos->rho2omega(cosmos->tau2rho_cdm(tau),tau);
  double Omega_b = cosmos->rho2omega(cosmos->tau2rho_b(tau),tau);
  //double Omega_q = quintcosmos->rho2omega(quintcosmos->tau2rho_q(tau),tau);
  double w=quintcosmos->tau2w_q(tau);
  double Q = 1.0/(3.0*Omega_cdm+8.0*(Omega_nu + Omega_gamma));
  double V = (105-45*w+4*Omega_nu*(3*w-1))/(36*(w-1));

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
    Dgq = 3.0/4.0*(1.0+quintcosmos->tau2w_q(tau))*Dgp;
    Vp = -5.0/4.0*P*x*Dgp;
    Vb = Vp;
    Vn = Vp;
    Vc = -5.0/4.0*P*x*Dgp;
    Vq = -5.0/4.0*P*x*Dgp;
    Pi_nu = -Dgp*P*x*x;
    Phi = (5.0+2.0*Omega_nu)/(30.0+8.0*Omega_nu)*Dgp;
    Psi= -Phi-Omega_nu*Pi_nu/(x*x);
    normalize = Psi * (4.0*Omega_nu + 15.0)/10.;
    normalize = -1./normalize;



    break;
  case ControlPanel::isoCDM:   // CDM isocurvature initial conditions

    Dgc = 1.0;
    Dgp = 0*Dgc;
    Dgn = Dgp;
    Dgb = 3.0/4.0*Dgp;
    Dgq = Omega_cdm*(15 + 2*Omega_nu)*(1+w)*U;

    Vc = (4.0*Omega_nu-15.0)*Omega_cdm*U*x*Dgc/12.0;
    Vp = -15.0/4.0 * Omega_cdm*U*x;
    Vn = Vp;
    Vb = Vp;
    Vq = Omega_cdm*U*V*x;
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
    Dgq = 4*Omega_b*(15+2*Omega_nu)*(1+quintcosmos->tau2w_q(tau))*T*U*Dgb;
    Vc = Omega_b*(4*Omega_nu-15)/3*T*U*Dgb*x;
    Vp = -15*Omega_b*T*U*Dgb*x;
    Vn = Vp;
    Vb = Vp;
    Vq = Omega_b*(105-45*quintcosmos->tau2w_q(tau)+4*Omega_nu*(3*quintcosmos->tau2w_q(tau)-1))/(9*(quintcosmos->tau2w_q(tau)-1))*T*U*Dgb*x;
    Pi_nu = -8*Omega_b*T*U*Dgb*x*x;
    Phi = Omega_b*(15+4*Omega_nu)*U/4*Dgb;  
    Psi= -Phi-Omega_nu*Gpi8()*Pi_nu/(x*x);
    break;

  case ControlPanel::isoNeutrino:
    Dgp= 1.0;
    Dgn = -Omega_gamma/Omega_nu*Dgp;
    Dgc = 3.0/4.0*Dgp;
    Dgb = 3.0/4.0*Dgp;
    Dgq = 0;
    Vc = Omega_gamma*P*x*Dgp;
    Vp = (15.0+4.0*Omega_gamma+4.0*Omega_nu)/4.0*P*x*Dgp;
    Vn = -15.0/4.0*Omega_gamma/Omega_nu*P*x*Dgp;
    Vb = Vp;
    Vq = Omega_gamma*P*x*Dgp;
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
    Dgq = 3.0/4.0*(1.0+quintcosmos->tau2w_q(tau))*ad +  4*Omega_cdm*(15+2*Omega_nu)*(1+quintcosmos->tau2w_q(tau))*T*U*icdm+ 4*Omega_b*(15+2*Omega_nu)*(1+quintcosmos->tau2w_q(tau))*T*U*ibaryon+ 0*ineutrino;
    Vp = -5.0/4.0*P*x*ad  -15*Omega_cdm*T*U*x*icdm -15*Omega_b*T*U*x*ibaryon + (15.0+4.0*Omega_gamma+4.0*Omega_nu)/4.0*P*x*ineutrino ;
    Vb = -5.0/4.0*P*x*ad  -15.0*Omega_cdm*T*U*x*icdm -15.0*Omega_b*T*U*x*ibaryon + (15.0+4.0*Omega_gamma+4.0*Omega_nu)/4.0*P*x*ineutrino;
    Vn = -5.0/4.0*P*x*ad  -15.0*Omega_cdm*T*U*x*icdm -15.0*Omega_b*T*U*x*ibaryon  -15.0/4.0*Omega_gamma/Omega_nu*P*x*ineutrino ;
    Vc = -5.0/4.0*P*x*ad + Omega_cdm*(4*Omega_nu-15)/3*T*U*x*icdm+  Omega_b*(4*Omega_nu-15)/3*T*U*x*ibaryon + Omega_gamma*P*x*ineutrino ;
    Vq = -5.0/4.0*P*x*ad + (105.0-45.0*quintcosmos->tau2w_q(tau)+4*Omega_nu*(3*quintcosmos->tau2w_q(tau)-1))/(9.0*(quintcosmos->tau2w_q(tau)-1))*Omega_cdm*T*U*x*icdm + (105.0-45.0*quintcosmos->tau2w_q(tau)+4*Omega_nu*(3*quintcosmos->tau2w_q(tau)-1))/(9.0*(quintcosmos->tau2w_q(tau)-1))*Omega_b*T*U*x*ibaryon+  Omega_gamma*P*x*ineutrino;
    Pi_nu = -ad*P*x*x -Omega_cdm*8*T*U*x*x*icdm -Omega_b*8*T*U*x*x*ibaryon + 3.0*Omega_gamma/Omega_nu*P*x*x*ineutrino ;
    Phi = (5.0+2.0*Omega_nu)/(30.0+8.0*Omega_nu)*ad + Omega_cdm*(15.0+4.0*Omega_nu)*T*U*icdm  + Omega_b*(15.0+4.0*Omega_nu)*T*U*ibaryon+ Omega_gamma/(15.0+4.0*Omega_nu)*ineutrino ; 
    Psi= -Phi-Omega_nu*Gpi8()*Pi_nu/(x*x);
    break;
  default:
    cout << "Initial condition not yet implemented" << endl;
    throw exception();
    break;
  }


  // Normalize all variables such that the curvature chi is =-1
  // Needed for correct initial Powerspectrum Amplitude

  Dgp *= normalize;
  Dgn *= normalize;
  Dgb *= normalize;
  Dgc *= normalize;
  Dgq *= normalize;
  Vp *= normalize;
  Vb *= normalize;
  Vn *= normalize;
  Vc *= normalize;
  Vq *= normalize;
  Pi_nu *= normalize;
  Phi *= normalize;
  Psi *= normalize;

  // Hm. Well. Done. Sorry for the inconvenience caused ...

  //cout << "finitial: " << Psi << "   " << Dgb << endl;

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

  // fill dark energy. Please note that for DE, we use longitudinal gauge 
  y[qidx]= Dgq - 3.*(1.+w)*Phi;
  y[qidx+1]= Vq * (1.0 + w);
  if (noDarkEnergyPerturbations) {
    y[qidx] = y[qidx+1] = 0;
  }

  double mu2 =  quintcosmos->quintessence()->speedOfSound2(); // speed of sound in DE rest frame
  double wdot =  quintcosmos->tau2wdot_q(tau);
  double Pressure = 0; // Pressure
  double Dq = y[qidx];
  double phi = quintcosmos->Tau2Phi->fastY(tau);
  double phiDot = quintcosmos->Tau2PhiDot->fastY(tau);
  double rhoq = quintcosmos->tau2rho_q(tau);
  double X = Vq/k*phiDot;
  double adotoa=tau2adotoa(tau);
  if (IsPhantomCrossing && fabs(1+w) < 1e-2) {
    Pressure = mu2*Dq; // cs^2 always for crossing dark energy models
  } else {
    Pressure =  mu2*Dq + Vq*(1+w)/k*( 3.0*adotoa*(mu2 - w) +  wdot/(1+w));
  }
  const double Ucommaq = quintcosmos->quintessence()->Vprime(phi, tau2a(tau), adotoa);
  double a2 = tau2a(tau)*tau2a(tau);
  double Xprime = a2/phiDot*(Pressure*rhoq+Psi/a2*phiDot*phiDot+Ucommaq*X);
  Xprime = phiDot/(1.+w)*(Dgq-X*Ucommaq/rhoq)+phiDot*(Psi-3.*Phi);
  //double Xprime_inicheck = (y[qidx]/(1.+w)+Psi); //-Vq/k*(a2*Ucommaq/phiDot);
  //double Xprime_inicheckb = 1./(1.+w)*(Dgq-X*Ucommaq/rhoq)+(Psi-3.*Phi);

  // field perturbations:
#if 0
  y[xidx] = X;
  y[xidx+1] = Xprime;
#else
  y[xidx] = 0;
  y[xidx+1] = 0;
#endif

}


void SpeedyDEInvariant::getReady(const ControlPanel& control) {
 SpeedyInvariant::getReady(control);  // now call parent objects get ready
 qidx = nvar-1;
 xidx = nvar-3;
}


void SpeedyDEInvariant::fillPhiDot(double tau)
{
  return; // PhiDot and PsiDot are already calculated by fderivs
}


void SpeedyDEInvariant::calcPerturbationNr(const ControlPanel &control) {
  //   Calculate number of equations
  if (control.scalar) {
    //  iq0 = (lmaxg << 1) + 9 + lmaxnr;
    //iq1 = iq0 + nqmax;
    //iq2 = iq1 + nqmax;

    nvar = 4;  // CDM and Baryon Delta and Velocity
    nvar += lmaxg +1; // Photons  0.. lmaxg
    nvar += lmaxg -1; // Polarization 2 .. lmaxg
    if (cosmos->nuR() > 0) {
      nvar += lmaxnr +1; // massless neutrinos
    }
    if (cosmos->haveMassiveNeutrinos()) {
      nvar += MassiveNeutrinos::pertArraySize();
    }
  } else nvar = 0;

  nvar += 2; // 2 more for dark energy
  nvar += 2; // for derivatives of field perturbation

  if (control.tensor) nvart = 2*lmaxt + 4 + lmaxnr+1; else nvart = 0;

  y = new double[nvar+1];
  yprime = new double[nvar+1];
  yt = new double[nvart+1];
  ytprime = new double[nvart+1];
}
