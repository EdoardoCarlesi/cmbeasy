#include "speedycoupledinvariant.h"
#include "coupledquintcosmos.h"
#include "coupledleaping.h"
#include "coupledinvariant.h"
#include "exponentialcoupling.h"

#include <algorithm>
#include <limits>
#include <vector>
// uncomment this, if you would like to write to a file during fderivs
//X #define OUTPUTFILE

SpeedyCoupledInvariant::SpeedyCoupledInvariant(QuintCosmos *cosmos) 
: SpeedyDEInvariant(cosmos), quintcosmos(cosmos)
{
}

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

void SpeedyCoupledInvariant::fderivs(const double tau, const double *y, double *yprime) {
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
#warning fadeout disabled for all coupling models
  FadeOut = 1.0;
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
//X   double Dc = y[3];
  double Vc = y[4];
//X   double delta_c_longit = Dc-3.*adotoa*Vc/k;

  quintcosmos->quintessence()->touch(tau);
  // dark energy - speedydeinv original
//X   double delta_q_longit = y[qidx]; // please note that this is in longitudinal gauge, not gauge invariant like the rest
//X   double onepluswVq = y[qidx+1]; // gauge invariant = longitud for V, however we have (1+w)*V here (better conditioned)

  double onepw = 1.+quintcosmos->tau2w_q(tau);

//X   double Dq = y[qidx];
  double Dgq = y[qidx];
  double Vq = y[qidx+1];

//X   double delta_q_longit = Dq-3.*onepw*adotoa*Vq/k;
  double onepwVq = Vq*onepw;


  // Dg and V for photons
  double Dgp = 4*M[0];
  double Vp = M[1];
  double Dgn = 4*N[0];
  double Vn = N[1];

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
  double dgrho = tau2rho_cdm(tau) *(Dgc + 3.*a2ak*Vc);  //we collect all species, first: cdm
//X   double dgrho = tau2rho_cdm(tau)*(delta_c_longit + 3.*a2ak*Vc);  //we collect all species, first: cdm
  dgrho += tau2rho_b(tau) * (Dgb + 3.*a2ak*Vb);              //baryons
  dgrho += tau2rho_g(tau) * (Dgp + 4.*a2ak*Vp);              //photons
  dgrho += tau2rho_nu(tau) * (Dgn + 4.*a2ak*Vn);           // massless neutrinos
#warning cleanup needed
  // the dark energy  fluid has a special formula, because we propagate not D_g, but delta_longit
  // please do also note that V_dark has the (1+w) built in already
//X    dgrho += quintcosmos->tau2rho_q(tau)*(delta_q_longit + 3.*a2ak*onepwVq); 
   dgrho += quintcosmos->tau2rho_q(tau)*(Dgq + 3.*a2ak*onepwVq); 
  // there is not Phi Contrib from the DE fluid, because we did not use D_g, but delta_longit...
  // hence no need to subtract the Phi that came in to make delta gauge invariant
  double phicontrib = 3.*(tau2rho_b(tau) + tau2rho_cdm(tau)) + 4.*(tau2rho_g(tau) + tau2rho_nu(tau));
  // and none from cdm either, if we also use delta_c_longit
//X   double phicontrib = 3.*tau2rho_b(tau) + 4.*(tau2rho_g(tau) + tau2rho_nu(tau));
  phicontrib += 3.*onepw*quintcosmos->tau2rho_q(tau);
  Phi = dgrho / (k*k/Gpi4a2 + phicontrib);  // note: Gpi4a2 = a^2/(2 Mpl^2) 

  // then we collect all the shear sources there are
  double Pip = 12.0/5.0*M2;  // photons
  double Pin = 12.0/5.0*N[2];   // neutrinos

  // then we get Psi, first we get the total shear (times perturbation)
  double dgpres = cosmos->tau2p_g(tau)*Pip + cosmos->tau2p_nu(tau)*Pin;
  dgpres *= 2.0*Gpi4a2;
  Psi = -Phi - dgpres / k2;

  double delta_c_longit = Dgc - 3.* Phi;
  double delta_q_longit = Dgq - 3.*onepw* Phi;
  double delta_g_longit = Dgp - 4.* Phi;
  double delta_n_longit = Dgn - 4.* Phi;

  double R = cosmos->tau2R(tau);  // 4/3 * rho_photon / rho_baryon
  double OneOverRPlusUnity = 1.0/(1+R);
  double cs2 = soundSpeed(tau);  // baryon soundspeed (almost 0, for tau > 1e-2) 

  // d/dtau for Dgb
  yprime[1] = -k*Vb - 3*cs2*adotoa*Dgb ; 
//X #warning coupledinvariant has the following:
//X   yprime[1] = -k*Vb ; 

  // and d/dtau for M[0], which is 1/4 of Dgp 		  
  Mprime[0] = -K/3.0 * Vp;

//X   cout << "Coupled: " << boolalpha <<  CoupledPhoton << " " << CoupledBaryon << " - " << CoupledMultipole << " - " << CoupledOctopole << endl;

  double slip = 0;
  // get the slip for as long as either baryons or photons can be treated as tightly coupled
  double Vb_prime=0, Vp_prime=0;  // Photon and baryon velocity prime. Init with 0, but doesnt matter
  if (CoupledBaryon || CoupledPhoton) {
    // First we need the slip just as usual 
    // the difference is that we don't use it directly but plug it into some other equation
    // the slip is d/dtau ( V_b - V_photon)  the gauge-invariant formulation of eqn 74 in Ma & Bertschinger
    // plus a Quadrupole term that is neglected there...
    // in principle we need two additional things, (d^2/dtau^2 a)/a  **and** d/dtau phi
    // (rho + p)*V of all species and from this  d/dtau Phi
    double RhoPV = tau2rho_cdm(tau) * Vc + tau2rho_b(tau)*Vb;
    RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
    RhoPV += quintcosmos->tau2rho_q(tau)* onepwVq; // simple :-)
    double PhiDot  = adotoa*Psi - Gpi4a2*(RhoPV)/k;
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

  // introduce some quantities of use for nonzero coupling (from Kodama & Sasaki modified to include field-
  // dependent couplings)
  CoupledQuintCosmos* coupledquintcosmos = dynamic_cast<CoupledQuintCosmos*>(quintcosmos);
  if (!coupledquintcosmos) throw Bad_Error("coupledinv speedy - error");

  Coupling* coupling = coupledquintcosmos->coupling();

  Quintessence* quint = dynamic_cast<Quintessence*>(quintcosmos->quintessence());
  if (!quint) throw Bad_Error("coupledinv speedy casting error");

  const double q = quint->q(a);
  const double qdot = quint->qDot(a);  // d/dtau of background field value
  const double rho_c = tau2rho_cdm(tau);
  const double rho_q = coupledquintcosmos->tau2rho_q(tau);

  const double Q_c_0 = coupling->phi2Q_c_0(q, qdot, rho_c);
  const double Q_q_0 = coupling->phi2Q_phi_0(q, qdot, rho_c);
  double adot= tau2adot(tau);
  double Mpl = cosmos->M_p();
  double wq= coupledquintcosmos->tau2w_q(tau);
  double S = coupling->phi2S(quint->q(a), rho_c);
  double Sdot = coupling->phi2Sprime(q, qdot, rho_c, adotoa);

  double X = Vq/k * qdot;
  const double DS = coupling->phi2DS(q, qdot, rho_c, Dgc, adotoa, Phi, X);

  double hq = (1.+wq)*rho_q;
  double Vqprime = k/onepw*Dgq+(2.*adotoa + Q_q_0 / hq)*Vq-3.*k*(1.+Q_q_0/3./adotoa/hq)*Phi+k*Psi;
  const double Ucommaq = quint->Vprime(quint->q(a),a,adotoa);
  double Xprime = -2.*adotoa*qdot-a*a*Ucommaq-a*a*S;
  Xprime *= Vq/k;
  Xprime += qdot*Vqprime/k;
#warning check again with Valeria that this A is correct
  double A = a*a/(2.*k*adotoa*Mpl*Mpl)*(tau2rho_b(tau)*Vb+tau2rho_cdm(tau)*Vc+4./3.*tau2rho_g(tau)*Vp+
                    4./3.*tau2rho_nu(tau)*Vn+(1.+wq)*coupledquintcosmos->tau2rho_q(tau)*Vq)-a*a/(2.*adotoa*adotoa*Mpl*Mpl)*Phi*phicontrib/3.0;

  const double Q_c_0prime = coupling->phi2Q_c_0prime(q, qdot, rho_c, a, adotoa);
  const double Qprime_c = Q_c_0prime;

  yprime[3] = Q_c_0/rho_c*Dgc-k*Vc-Q_c_0/rho_c*A + qdot/rho_c*DS + S*Xprime/rho_c;
  yprime[3] += -a*Qprime_c/rho_c /adotoa*Phi +Q_c_0/rho_c*Psi;
  yprime[4] = -(adotoa+S*qdot/rho_c)*Vc+S*qdot/rho_c*Vq+k*Psi;

  // dark energy
  double w = quintcosmos->tau2w_q(tau);
  double wdot =  quintcosmos->tau2wdot_q(tau);
  double mu2 =  quintcosmos->quintessence()->speedOfSound2(); // speed of sound in DE rest frame
  double P_Anisotropic = 0; // no anisotropic pressure for this dark energy
  double P = 0; // Pressure
  if (IsPhantomCrossing && fabs(1+w) < 1e-2) 
    P = mu2*delta_q_longit; // cs^2 always equal 1 for crossing dark energy models close to the crossing
  else  P =  mu2*delta_q_longit + onepwVq/k*( 3.0*adotoa*(mu2 - w) +  wdot/(1.+w));

  double RhoPV = tau2rho_cdm(tau) * Vc + tau2rho_b(tau)*Vb;
  RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
  RhoPV += quintcosmos->tau2rho_q(tau)* onepwVq; // simple :-)
  PhiDot  = adotoa*Psi - Gpi4a2*(RhoPV)/k;

  const double pqdot = quintcosmos->tau2p_qdot(tau);
  const double rhoqdot = quintcosmos->tau2rho_qdot(tau);
  double cs_q = pqdot/rhoqdot;

  const double Q_q_0prime = coupling->phi2Q_phi_0prime(q, qdot, rho_c, a, adotoa);
  double Qprime_q = adotoa/a*Q_q_0-Q_q_0prime/a;

  Qprime_q = Q_q_0prime;

  yprime[qidx] = -(3.*adotoa*(1.-wq)-Q_q_0/rho_q)*Dgq
                 -k*(1.+wq)*Vq-3.*adotoa/k*(1.-cs_q)*(3.*adotoa+Q_q_0/rho_q)*Vq
                 +a/rho_q*(-qdot/a*DS-1./a*S*Xprime-1./a*Q_q_0*(A-Psi)
                           -Qprime_q*Phi/adotoa);

  yprime[qidx+1] = Vqprime;

  // CAREFUL with FADEOUT for models with mu2  closer to 0 !!! Only allowed
  // if perturbations decay relativistically within the horizon

  if (mu2 > 1e-2) { 
    yprime[qidx] *= FadeOut;
    yprime[qidx+1] *= FadeOut;
  }

 static bool overone = false;
 static bool over1e4 = false;

 if (a>1 && !overone)
	 /* commented to speed up */
 //  cout << "Dgc/a(a=1) = " << (Dgc/a) << endl;
 if (/*a>1e-4 &&*/ !over1e4)
 {
	 /* commented to speed up */
  // cout << Psi << " " << Phi << " " << quintcosmos->tau2phi(tau) << endl;
   //cout <<  "rho_r " << (quintcosmos->tau2rho_relativistic(tau)/cosmos->M_p()/cosmos->M_p()) << endl;
   //cout << Dgc << " is Dgc" << endl;
   over1e4 = true;
 }

 static int points = 0;
 if (!(0 == points++ % 10))
   return;

  double conversion_rho_c = rho_c;
  double convQ_c_0 = Q_c_0;
  double convQ_q_0 = Q_q_0;
  static bool fixPhi = false;
  if (false && !fixPhi) {
    conversion_rho_c=35014.4889;
    convQ_c_0 = -0.05*conversion_rho_c*0.000471703521;
     /* commented to speed up */
    //cout << "rho_c is: " << (rho_c/Mpl/Mpl) << endl;
    //cout << "q_c_0 is: " << (Q_c_0/Mpl/Mpl) << endl;
    //cout << "convQ_c_0 is: " << (convQ_c_0) << endl;
    //cout << "qdot is: " << (qdot/Mpl) << endl;
    //cout << "frac: " << (Q_c_0/(adotoa*3.*rho_c)) << endl;
  }
  const double Dq = Dgq - 3.*(1.+wq)*(1.+convQ_q_0/(3.*adotoa*rho_q*(1.+wq)))*(Phi-adotoa*Vq/k);
  const double Dc = Dgc - 3.*(1.+convQ_c_0/(3.*adotoa*conversion_rho_c))*(Phi-adotoa*Vc/k);
  double Dp = delta_g_longit + 4.*adotoa*Vp/k;
  double Dn = delta_n_longit + 4.*adotoa*Vn/k;

  if (!fixPhi) {
    fixPhi = true;
     /* commented to speed up */
    //cout << "Dc" << Dc << " but we put in: " <<0.947036614 << endl;
    //cout << "Q_q_0: " << (Q_q_0/Mpl/Mpl) << endl;
    //cout << "yprime[qidx]: " << (yprime[qidx]) << endl;
    //cout << "Xprime: " << (Xprime/Mpl) << endl;
    //cout << "qprime from klein-gordon: " << ((-a*a*Ucommaq-a*a*S)/-2./adotoa/cosmos->M_p()) << endl;
  }
 
  double Psi2nd = 12.0/5.0*a2*(tau2rhodot_g(tau)*M[2] + tau2rhodot_nu(tau)*N[2] 
			      + tau2rho_g(tau)*Mprime[2] + tau2rho_nu(tau)*Nprime[2])/3.0;
 // add to this d/dtau a^2 term
 Psi2nd += 2.0*a*tau2adot(tau)*(tau2p_g(tau)*Pip + tau2p_nu(tau)*Pin);
 double PsiDot_n = -PhiDot - Gpi8()/k2 * Psi2nd; 

}

void SpeedyCoupledInvariant::initialScalarPerturbations(const ControlPanel &control, const double tau)
{
//X   cicheck.initialScalarPerturbations(control, tau);
//X 
  // copy the phantomcrossing flag for fderivs()
  IsPhantomCrossing = control.isPhantomCrossing();

#ifdef OUTPUTFILE   
  ofs = new ofstream("/tmp/robbers/valeria/speedycolddark.dat");
  ofs2 = new ofstream("/tmp/robbers/valeria/georg-bkg.dat");
#endif

  // unfortunatly, we need to define our M,N,E ourselves, cause
  // they are non -constant (we have to initialize them)
  // see setPointers() for more
  double *M = &y[5];                    // Temperature
  double *E = &M[lmaxg-1];  // E[2] is first after M[lmaxg]
  double *N = &E[lmaxg+1];  //  massles neutrinos 

  //  cout << "N[lmaxnr] " << &N[lmaxnr] << "  y[nvar]: " << &y[nvar] << endl;
  double Dgp, Dgc, Dgb, Dgn, Dgq, Vb,Vc,Vp,Vn,Vq,Psi,Pi_nu;
  k2 = k*k;

  double x = k*tau;

  double Omega_nu = cosmos->rho2omega(cosmos->tau2rho_nu(tau),tau);
  double Omega_gamma = cosmos->rho2omega(cosmos->tau2rho_g(tau),tau);
  double Omega_cdm = cosmos->rho2omega(cosmos->tau2rho_cdm(tau),tau);
  double Omega_b = cosmos->rho2omega(cosmos->tau2rho_b(tau),tau);
  double Omega_q = cosmos->rho2omega(quintcosmos->tau2rho_q(tau),tau);
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

  CoupledQuintCosmos* cqCosmos = dynamic_cast<CoupledQuintCosmos*>(quintcosmos);
  if (!cqCosmos || (control.initialCond != ControlPanel::adiabatic))
   throw Bad_Error("SpeedyCoupledInvariant::initialScalarPerturbations - only combination supported is adiabatic+coupling");

  Coupling* coupling = cqCosmos->coupling();

  Quintessence* quint = dynamic_cast<Quintessence*>(quintcosmos->quintessence());
  if (!quint) throw Bad_Error("coupledinv speedy casting error");
  quint->touch(tau);

  double a= tau2a(tau);
  double adotoa = tau2adot(tau)/a;

  const double q = quint->q(a);
  const double qdot = quint->qDot(a);  // d/dtau of background field value
  const double rho_c = tau2rho_cdm(tau);
  const double rho_q = cqCosmos->tau2rho_q(tau);
  const double wq = cqCosmos->tau2w_q(tau);
  const double hq = (1.+wq)*rho_q;

/*
  if (wq < (0.9*1./3.)) {
    cout << endl << "k=" << k << ": at z=" << cosmos->tau2z(tau) << " w_q=" << wq << endl;
    throw Bad_Error("SpeedyCoupledInvariant::initialScalarPerturbations - initial conditions only valid for w_q=1/3");
  }
*/
  const double Q_c_0 = coupling->phi2Q_c_0(q, qdot, rho_c);
  const double Q_q_0 = coupling->phi2Q_phi_0(q, qdot, rho_c);
  const double S = coupling->phi2S(quint->q(a), rho_c);
  const double Shat = S / rho_c;

#warning make initial conditions more general (w_x from background)
  const double w_eff_p = 1./3.;
  const double w_eff_nu = 1./3.;
  const double w_eff_b = 0;
  const double w_eff_c = Q_c_0/(3.*adotoa*rho_c);
  const double w_eff_q = wq + Q_q_0/(3.*adotoa*rho_q);

  double A1 = std::numeric_limits<double>::quiet_NaN();
  const double A2 = 1.+ Shat * Omega_cdm/Omega_q * qdot/(1.+wq);
  const double B1 = -(2.+Shat*qdot);

  // for Phi and Psi
  double Gpi4 = Gpi8()*0.5;
  double Gpi4a2 = Gpi4*a*a;
  double a2ak = adotoa / k;
  double dgrho = 0;
  double phicontrib = 0;
  double PhiDenom = 0;
  double dgpres = 0;
  double Pip, Pin;
  double Pi_nu_hat = 0;


  switch (control.initialCond) {
  case ControlPanel::adiabatic:
    Dgc = 1.;
    Dgp = (1.+w_eff_p)/(1.+w_eff_c)*Dgc;
    Dgn = (1.+w_eff_nu)/(1.+w_eff_p)*Dgp;
    Dgb = (1.+w_eff_b)/(1.+w_eff_c)*Dgc;
    Dgq = (1.+w_eff_q)/(1.+w_eff_c)*Dgc;
    A1 = Dgq/(1.+wq);
    Pi_nu_hat = -Dgp*P;
    Pi_nu = Pi_nu_hat*x*x;
    Vn = 5./4.*Pi_nu_hat;
    Vn *= x;
    Vp = (Dgp-Dgn)/4.+5./4.*Pi_nu_hat;
    Vp *= x;
    Vb = Vp;
    Vq = -A1+(3.+A2)*Dgn/8.-1./8.*(15.+4.*Omega_nu+A2*(5.+4.*Omega_nu))*Pi_nu_hat;
    Vq *= x/A2;
    Vc = (B1+2.)/B1*Vq+1./8.*Dgp/B1*x-1./(2.*B1)*Vp+1./(2.*B1)*Omega_nu*Pi_nu_hat*x;

    control.initialConditionFactors.apply(Dgp, Dgn, Dgb, Dgc, Vp, Vn, Vb, Vc);

    Phi = Omega_cdm *(Dgc + 3.*(1.+0.)*Vc/x);
    Phi += Omega_b *(Dgb + 3.*(1.+0.)*Vb/x);
    Phi += Omega_gamma *(Dgp + 3.*(1.+1./3.)*Vp/x);
    Phi += Omega_nu *(Dgn + 3.*(1.+1./3.)*Vn/x);
    Phi += Omega_q*(Dgq + 3.*(1.+wq)*Vq/x);

    PhiDenom = Omega_cdm*(1.+0.);
    PhiDenom += Omega_b*(1.+0.);
    PhiDenom += Omega_gamma*(1.+1./3.);
    PhiDenom += Omega_nu*(1.+1./3.);
    PhiDenom += Omega_q*(1.+wq);

    Phi /= 3.*PhiDenom+2./3.*x*x;
    Psi = -Phi-Omega_nu*Pi_nu_hat;

    normalize = Psi * (4.0*Omega_nu + 15.0)/10.;
    normalize = -1./normalize;

#if 0
    normalize = 1./Dgp;
#endif

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
  // Hm. Well. Done. Sorry for the inconvenience caused ...
  
  // fill matter
  y[1] = Dgb;
  y[2] = Vb;
  y[3] = Dgc;
  y[4] = Vc;

  // fill photons and neutrinos
  M[0] = Dgp / 4.0;
  M[1] = Vp;
  N[0] = Dgn / 4.0;
  N[1] = Vn;
  N[2] = 5.0/12.0*Pi_nu;
  // set all higher moments to 0, obs: E[2] is the first one needed
  for (int l = 2; l <= lmaxg; l++) {M[l] =0; E[l] = 0;}
  for (int l = 3; l <= lmaxnr; l++) N[l] =0; 
  
  // fill dark energy. Please note that for DE, we use longitudinal gauge 
  // for coupled cosmologies, we're propagating Dgq again (not delta_q_longit, as in the parent
  // class)
  y[qidx]= Dgq; // - 3*(1+w)*Phi;
//X   y[qidx] = Dgq-3.*(1.+w)*Phi+3.*(1.+w)*adotoa*Vq/k;
  y[qidx+1]=Vq; // * (1.0 + w);
} 

void SpeedyCoupledInvariant::getReady(const ControlPanel& control) {
  SpeedyDEInvariant::getReady(control);  // now call parent objects get ready

  qidx = nvar-1;
}

/*!
  For the sources, to be precise the ISW effect, one needs d/d tau
  of Phi and Psi. In the QuintSpeedyDEInvariant class, these are needed anyhow
  within the fderivs() function. However, here they are only calculated
  for the sources, which is why we need to do it.

  Please note that fillPhiDot() assumes that it is called after fderivsWrapper() 
  has been called, i.e. Phi and Psi and all other stuff is up to date.
*/
void SpeedyCoupledInvariant::fillPhiDot(double tau) {
 double Vb = y[2];
 double Vc = y[4]; 
 double Vp = M[1];
 double Vn = N[1];
 double Vq = y[qidx+1];
 double a= tau2a(tau);   // we take it from cosmos :-), as this is also mother, children will be fast
 double a2 = a * a;
 double adotoa = tau2adot(tau)/a;
 
 double Gpi4 = Gpi8()*0.5;
 double Gpi4a2 = Gpi4*a2;
 
 double Pip = 12.0/5.0*M[2];  // photons
 double Pin = 12.0/5.0*N[2];   // neutrinos
 
 // (rho + p)*V of all species and from this and eqn 2.50 d/dtau Phi
 double RhoPV = tau2rho_cdm(tau) * Vc + tau2rho_b(tau)*Vb;
 RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
 double onepwVq = Vq*(1.+quintcosmos->tau2w_q(tau));
 RhoPV += quintcosmos->tau2rho_q(tau)* onepwVq;
 PhiDot  = adotoa*Psi - Gpi4a2*(RhoPV)/k;
 // now d/dtau Psi from d/dtau of eqn 2.52, first d/dtau (p*Pi) p = rho/3
 double Psi2nd = 12.0/5.0*a2*(tau2rhodot_g(tau)*M[2] + tau2rhodot_nu(tau)*N[2] 
			      + tau2rho_g(tau)*Mprime[2] + tau2rho_nu(tau)*Nprime[2])/3.0;
 // add to this d/dtau a^2 term
 Psi2nd += 2.0*a*tau2adot(tau)*(tau2p_g(tau)*Pip + tau2p_nu(tau)*Pin);
 PsiDot = -PhiDot - Gpi8()/k2 * Psi2nd; 
}

void SpeedyCoupledInvariant::calcPerturbationNr(const ControlPanel &control) {
  //   Calculate number of equations 
  if (control.scalar) {
    //  iq0 = (lmaxg << 1) + 9 + lmaxnr;
    //iq1 = iq0 + nqmax;
    //iq2 = iq1 + nqmax;
    
    nvar = 4;  // CDM and Baryon Delta and Velocity
    nvar += lmaxg +1; // Photons  0.. lmaxg
    nvar += lmaxg -1; // Polarization 2 .. lmaxg
    nvar += lmaxnr +1; // massless neutrinos 
  } else nvar = 0;
 
  if (control.tensor) nvart = 2*lmaxt + 4 + lmaxnr+1; else nvart = 0;

  nvar += 2; //  2 more for dark energy
 
  y = new double[nvar+1];
  yprime = new double[nvar+1];
  yt = new double[nvart+1];
  ytprime = new double[nvart+1];
}
