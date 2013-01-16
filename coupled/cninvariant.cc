#include "cninvariant.h"
#include "cncosmos.h"
#include "massiveneutrinos.h"
#include "massiveneutrinosnew.h"
#include "perturbationtracker.h"
#include "stdio.h"
#include <iostream>
#include <sstream>
#include <limits>

// uncomment this, if you would like to write to a file during fderivs
//#define OUTPUTFILE

#define qmax 1000

string conversion(double d) {
std::ostringstream oss;
oss << d;
std::string value = oss.str();
return value;
}

double CnInvariant::maxDelta = std::numeric_limits<double>::infinity();
double CnInvariant::maxEta = std::numeric_limits<double>::infinity();

/* Unused MCMC cutoff variables */
double CnInvariant::maxGravPotential=std::numeric_limits<double>::infinity();
double CnInvariant::earliestStopZ=std::numeric_limits<double>::infinity();
double CnInvariant::cutoffIndex = std::numeric_limits<double>::infinity(); // default value - no scale dependence

int printVariable=0;
int printVariable2 =0;
static bool stop=false;	

CnInvariant::CnInvariant(CnCosmos *cosmos)
  :SpeedyDEInvariant(cosmos), quintcosmos(cosmos)
{
mDelta_c_longit = mDelta_b_longit = mV_b_longit = mV_c_longit = 0;
PhiDot = PsiDot = Psi = Phi = 0;
}

/*!
tau-derivatives of Gauge-CnInvariant Variables. y and yprime come
from odeint(), meaning that they usually will not coincide with our
objects y and yprime. However, fderivsWrapper() will at the end
of the integration fill y and yprime of our object with the result. 
In other words: scalarSources() etc will have a fully valid y and 
yprime arrary and when scalarSources() calls fderivs(), it does so
with our objects y and yprime

This is a very high efficiency version.
*/


void CnInvariant::fderivs(const double tau, const double *y, double *yprime)
{
  PerturbationTracker* mPertTracker = PerturbationTracker::self(this);
	//cout  << "called neutrinodriver.fderivs() " << endl; 
 std::fill(yprime, yprime+nvar+1, 0.);  // be safe (and not much slower)

called++;

if (tau < 1000) called3++;

double a = tau2a(tau);   // we take it from cosmos :-), as this is also mother, children will be fast
double a2 = a * a;
double adotoa = tau2adot(tau)/a;

const double w = quintcosmos->tau2w_q(tau);
const double wdot =  quintcosmos->tau2wdot_q(tau);

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
double Dq = y[qidx];   // please note that this is in longitudinal (Newtonian) gauge, not gauge invariant like the rest
double Vq = y[qidx+1]; // gauge invariant = longitud for V, however we have (1+w)*V here (better conditioned)

//const double onepwVq_longit = Vq;
//const double delta_q_longit = Dq;

double onepwVq_longit = Vq;
double delta_q_longit = Dq;

const double phi=quintcosmos->tau2phi(tau);
static const double Mp = quintcosmos->M_p();
static const double Mp2 = Mp*Mp;

double X = y[xidx];
double Xprime = y[xidx+1];
double X_true = X;
double Xprime_true = Xprime;

// Dg and V for photons
double Dgp = 4*M[0];
double Vp = M[1];

// massless neutrinos
double Dgn = 0;
double Vn = 0;
double Pin = 0;

if (cosmos->nuR() > 0) {
  Dgn = 4*N[0];
  Vn = N[1]; 
  Pin = 12.0/5.0*N[2];
}

double delta_nuNR_longit = 0;
double onepwNRVn_nr_longit = 0;
double w_nuNR = 0;

if (cosmos->haveMassiveNeutrinos()) {
// we are propagating massive neutrinos in longitudinal gauge
  w_nuNR = quintcosmos->tau2w_nuNR(tau);
}

double Vq_longit;
#warning (ehm, should this be so...)
static bool switched = false;
if (IsPhantomCrossing && fabs(w+1) < 1e-2) {
//X     if (!switched)
//X       cout << "switched at tau=" << tau << endl;
  switched = true;
  Vq_longit = 0.;
} else {
  switched = false;
  Vq_longit = Vq/(1.+w);
}

const double qdot = quintcosmos->Tau2PhiDot->fastY(tau);

double lna = log(a);
double Gpi8a2 = Gpi8()*a2;
double rho_mass_nu_0 = rho_nu0()*nuNR();
double plain_nu = /* Gpi8a2* */ rho_mass_nu_0*exp(-4.*lna);

  rho_nuNR = P_nu_NR = Pi_nu_NR = 0;
  del_rho_nu_NR = f_nu_NR = del_P_nu_NR = sigma_nu_nr = 0;
  double neutrinoMass_kbT = 0;
  double dlnmdphi=0;
  double csNuNR2 =0;
  double delPorho_nuNR_longit =0;
  double beta =0;
if (cosmos->haveMassiveNeutrinos()) {
  neutrinoMass_kbT = quintcosmos->mass_nu_kbT(phi);
  beta = dlnmdphi = neutrinoMass_kbT?quintcosmos->dMass_nu_kbTdphi(phi)/neutrinoMass_kbT:0.;

#warning check del_phi units
 MassiveNeutrinosNew::nu1(a*neutrinoMass_kbT, &rho_nuNR, &P_nu_NR);
 MassiveNeutrinosNew::nu2(a*neutrinoMass_kbT, dlnmdphi, rho_nuNR, P_nu_NR,
                        X/Mp/* /plain_nu */, &del_rho_nu_NR, &f_nu_NR,
                        &del_P_nu_NR, &sigma_nu_nr, NR_0, NR_1, NR);
  csNuNR2 = del_P_nu_NR/del_rho_nu_NR;
  delPorho_nuNR_longit = del_P_nu_NR/rho_nuNR;
  //mPertTracker->track("P_nu_NR", tau, P_nu_NR);
  //mPertTracker->track("del_P_nu_NR", tau, del_P_nu_NR);

  delta_nuNR_longit = del_rho_nu_NR/rho_nuNR;
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
  dgrho += tau2rho_g(tau) * (Dgp + 4*a2ak*Vp);              //photonsGamesPolitics
  dgrho += tau2rho_nu(tau) * (Dgn + 4*a2ak*Vn);           // massless neutrinos 
  // the dark energy  fluid has a special formula, because we propagate not D_g, but delta_longit
  // please do also note that V_dark has the (1+w) built in already
  dgrho += quintcosmos->tau2rho_q(tau)*(Dq + 3*a2ak*Vq); 
  // massive neutrinos are in longitudinal gauge, too
  
  dgrho += tau2rho_nuNR(tau)*(delta_nuNR_longit + 3.*a2ak*onepwNRVn_nr_longit);
  
  // there is not Phi Contrib from the DE fluid, because we did not use D_g, but delta_longit...(same for nu_NR)
  // hence no need to subtract the Phi that came in to make delta gauge invariant
  double phicontrib = 3*(tau2rho_b(tau) + tau2rho_cdm(tau)) + 4*(tau2rho_g(tau) + tau2rho_nu(tau));
  Phi = dgrho / (k*k/Gpi4a2 + phicontrib);  // note: Gpi4a2 = a^2/(2 Mpl^2) 
 
  double nuContribToPhi= tau2rho_nuNR(tau)*(delta_nuNR_longit + 3.*a2ak*onepwNRVn_nr_longit);
  nuContribToPhi /= (k*k/Gpi4a2 + phicontrib);
  
  //Phi += nuContribToPhi;

  // then we collect all the shear sources there are
  double Pip = 12.0/5.0*M2;  // photons

  double neutrinoShearContrib=0;
  if (cosmos->nuR()>0) {
    Pin = 12.0/5.0*N[2];   // massless neutrinos
    neutrinoShearContrib += cosmos->tau2p_nu(tau)*Pin;
  }

  if (cosmos->haveMassiveNeutrinos()) {
    const double p_nm = cosmos->tau2p_nuNR(tau);
    neutrinoShearContrib = Pi_nu_NR*p_nm;
  }

  // then we get Psi, first we get the total shear (times perturbation)
  double dgpres = cosmos->tau2p_g(tau)*Pip; // + cosmos->tau2p_nu(tau)*Pin;
  dgpres *= 2.0*Gpi4a2;
  dgpres += neutrinoShearContrib * 2.*Gpi4a2;
  Psi = -Phi - dgpres / k2;

  double RhoPV = tau2rho_cdm(tau) * Vc + tau2rho_b(tau)*Vb;
  RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
  RhoPV += quintcosmos->tau2rho_q(tau)* Vq; // simple :-)
  RhoPV += tau2rho_nuNR(tau)*onepwNRVn_nr_longit;
  PhiDot = adotoa*Psi - Gpi4a2*(RhoPV)/k;

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
    // (rho + p)*V of all species and from this  d/dtau Phi
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
  } else {
    Vb_prime = k*(Psi -3*cs2*Phi) + Dgb*k*cs2 - adotoa*Vb + Opac*R*(Vp-Vb);
  }

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

	// STOP CONDITIONS
  Quintessence* quint = quintcosmos->quintessence();
	double mQ = cosmos->M_p();
	bool stopCondition = delta_nuNR_longit >= maxDelta; // maxDelta is already normalized to int. units
	//cout << " delta_phi:  " << Dgq*(quint->q(a)/cosmos->M_p())*pow(10,-4.5) << endl;
	bool stopConditionEta = -beta*X*pow(Mp,-1) >= maxEta;  
	//if(stopCondition) cout << " Stop Condition " <<endl;
	//if(stopConditionEta) cout << " Stop Condition Eta " << endl; 
	if(stopCondition || stopConditionEta){
		// fix neutrino overdensity
	        //nuContribToPhi = (Gpi4a2/(k*k))*nuContribToPhi;
		//if (stopCondition) {
		//if (!stop) {		
		//maxEta = -beta*Dq*(quint->q(a)/cosmos->M_p());
		//stop = true;
		//}
		//delta_nuNR_longit = maxDelta; 
	//	y[nuNRidx] = maxDelta;	
	//yprime[nuNRidx] = mos->M_p()/(beta*quint->q(a));axDelta;
		// cout << " del_rho_nu_NR: " << del_rho_nu_NR <<  endl;
		// set velocities to zero
		onepwNRVn_nr_longit=0;
		// stop neutrino growth
		yprime[nuNRidx] = 0;
		// stop neutrino velocity growth
		//yprime[nuNRidx+1]=0;
		//mV_nu_nr_longit = 0;
		//}

		//if (stopConditionEta) {
		//X= -maxEta*pow(Mp,1)/beta;
		// stop velocities
			Xprime =0;
			//if (!stop) {
		//maxDelta = delta_nuNR_longit;
		//stop = true;
		//}
		// fix DE overdensity
		//yprime[qidx+1] = 0;
		//delta_q_longit = -maxEta / (beta*mQ);
		// set DE velocities to zero
		//Vq_longit = 0;
		// Stop overdensity growth
		//yprime[qidx] = 0;
		//yprime[qidx+1] = 0;
  		//mDelta_q_longit = Dq;
		//}
	}
//TODO add control on stop variable to cutoff higher neutrino moments
//	if(delta_nuNR_longit > maxDelta) stop=true;
  if (cosmos->haveMassiveNeutrinos()) {
      	 MassiveNeutrinosNew::propagateMassiveNeutrinoMoments(y, yprime, a, tau,
   neutrinoMass_kbT, k, PhiDot, Psi, beta, X/Mp, stop); // added boolean variable stop to control NU higher-moments cutoff
  }

  // dark energy
  const double dbetadphi = quintcosmos->dBetadphi(phi);
  const double rhonunr = quintcosmos->tau2rho_nuNR(tau);
  const double rhoq = quintcosmos->tau2rho_q(tau);
  const double wnunrdot = quintcosmos->tau2wDot_nuNR(tau);
  const double wqdot = wdot;
  const double wq = w;

  double mu2 =  quintcosmos->quintessence()->speedOfSound2(); // speed of sound in DE rest frame
  if (mu2 != 1) {
    throw Bad_Error("CnInvariant::fderivs() - mu2 != 1 not supported for this gauge.");
  }
  double P_Anisotropic = 0; // no anisotropic pressure for this dark energy
  double P = 0; // Pressure
  static bool pressureSwitched = false;

  const double Ucommaq = quintcosmos->quintessence()->Vprime(phi, a, adotoa);
  double pNuNR_ = tau2p_nuNR(tau);
  double qddot= -2.*adotoa*qdot-a*a*Ucommaq-a*a*beta*(rho_nuNR-3.*pNuNR_)/Mp;

  double altX = onepwVq_longit/(1.+w)/k*qdot;
  double t0 = qdot*Xprime/a2-Psi/a2*qdot*qdot-Ucommaq*altX;
  double t1  = a2/qdot*(P*rhoq);
  double t2  = a2/qdot*(-Psi/a2*qdot*qdot);
  double t3  = a2/qdot*(-Ucommaq*altX);
  double P_orig;
  P_orig = (rhoq==0?0:(qdot*Xprime/a2-Psi/a2*qdot*qdot-Ucommaq*X)/rhoq);
//#warning  not orig P;
  P = P_orig;

  double Xprime_from_P = a2/qdot*(P*rhoq+Psi/a2*qdot*qdot+Ucommaq*altX);

  double Psi2nd = 12.0/5.0*a2*(tau2rhodot_g(tau)*M2
      + tau2rho_g(tau)*Mprime[2])/3.0;
  if (cosmos->nuR()>0) {
    Psi2nd += 12.0/5.0*a2*(tau2rhodot_nu(tau)*N[2] + tau2rho_nu(tau)*Nprime[2])/3.0;
  }
  Psi2nd += 2.0*a*tau2adot(tau)*(tau2p_g(tau)*Pip + tau2p_nu(tau)*Pin);

  if (cosmos->haveMassiveNeutrinos()) {
#warning check nuder for varying beta
    double shearnudot =MassiveNeutrinosNew::nuder(a*neutrinoMass_kbT, adotoa, beta, rho_nuNR, qdot/Mp,
                                                NR, NRprime);
    double pNuNR = tau2p_nuNR(tau);
    double massless_neutrino_p = pNuNR/P_nu_NR;
    double wdot_nuNR = cosmos->Tau2wdot_nuNR->fastY(tau);
    double rhodot_nuNR = cosmos->Tau2RhoDot_nuNR->fastY(tau);
    double pdotNuNR = tau2rho_nuNR(tau)*wdot_nuNR+rhodot_nuNR*w_nuNR;
#warning figure out why the zero should be correct...
    double PDot_nu_NR = 0.*pdotNuNR/massless_neutrino_p;
    double PiDot_nu_NR = 1.5/P_nu_NR*shearnudot; //-pdotNuNR/pNuNR*Pi_nu_NR;
    Psi2nd -= 2.*a*tau2adot(tau)*pNuNR*Pi_nu_NR - a2*(0*pdotNuNR*Pi_nu_NR+pNuNR*PiDot_nu_NR);
  }
  PsiDot = -PhiDot - Gpi8()/k2 * Psi2nd;

  const double piTh = 0;
  double piTnu = Pi_nu_NR;

  double onepwVqprime  = 0;
  if (rhoq!=0) {
     onepwVqprime = (3.*wq-1.)*adotoa*onepwVq_longit
                      + k*P
                      + k * Psi * (1.+wq)
                      - 2./3.*wq*k*piTh
                      + rhonunr/rhoq*(1.-3.*w_nuNR)*beta*(qdot/Mp*onepwVq_longit-k*X/Mp);
  }
  double oldval = yprime[qidx];
	
  if (rhoq == 0) {
    yprime[qidx] = 0;
    yprime[qidx+1] = 0;
 } else { // Scalarfield overdensity evolution - see "Neutrino Clustering in Growing Neutrino Quintessence"

     	 yprime[qidx] = 3.*adotoa*(wq*delta_q_longit-P)
                   -beta*qdot/Mp*rhonunr/rhoq*((1.-3.*w_nuNR)*delta_q_longit-(delta_nuNR_longit-3.*delPorho_nuNR_longit))
                   -k*onepwVq_longit-3.*(1.+wq)*PhiDot
                   +rhonunr/rhoq*(1.-3.*w_nuNR)*(beta*Xprime + dbetadphi*qdot*X)/Mp;
    // since we propagate (1+w)Vq, this is the eq. for [(1+w)Vq'] = (1+w)V'+wqdot*V
    yprime[qidx+1] = onepwVqprime;
  }

  if (isinf(yprime[qidx])) {
      cout << "yprime[qidx]; " << yprime[qidx] << endl;
     throw Bad_Error("yprime[qidx]: isinf");
  }

  if (isnan(yprime[qidx]))
     throw Bad_Error("yprime[qidx]: isnan");

#warning fix phi dependent beta
   static const double d2lnmdp2 = 0;

   yprime[xidx] = Xprime;

   double Xdoubleprime = -2.*adotoa*Xprime - (k*k + a2 * quint->Vprime2(phi, a, adotoa, tau))*X
  + (PsiDot - 3.*PhiDot)*qdot - 2.*a2*quint->Vprime(phi, a, adotoa) * Psi
  -a2*(+dlnmdphi*rhonunr/Mp*(delta_nuNR_longit-3.*delPorho_nuNR_longit) //delta_nuNR_longit*(1.0-3.0*csNuNR2)
#warning check sign, (beta convention)
      +d2lnmdp2/Mp/*=dbetadphi*/*X*rhonunr*(1.-3.*w_nuNR)
                                 +2.0*dlnmdphi/Mp*rhonunr*(1.0-3.0*w_nuNR)*Psi);

  yprime[xidx+1]= Xdoubleprime;
	
	// Neutrino overdensity evolution 
	// 19.10.10 Changed signs to some parts according to the original article
  yprime[nuNRidx] = 3.*(adotoa + beta*qdot/Mp)*(w_nuNR*delta_nuNR_longit-delPorho_nuNR_longit)
                    - k*onepwNRVn_nr_longit
                    - 3.*(1.+w_nuNR)*PhiDot
                    +beta*(1.-3.*w_nuNR)*Xprime/Mp
                    +dbetadphi*qdot/Mp*X/Mp*(1.-3.*w_nuNR);

  // see the comment above for Vq'
  yprime[nuNRidx+1] = (1.-3.*w_nuNR)*(-beta*qdot/Mp-adotoa)*onepwNRVn_nr_longit
                      + k*delPorho_nuNR_longit
                      + k*(1.+w_nuNR)*Psi
                      -2./3.*k*w_nuNR*piTnu
                      + k*beta*X/Mp*(1.-3.*w_nuNR);

  const double nval1 = yprime[nuNRidx+1];

#warning not using fluid for neutrinos

  if (isinf(yprime[qidx])) {
    throw Bad_Error("isinf yprime[qidx]");
  }
  // CAREFUL with FADEOUT for models with mu2  closer to 0 !!! Only allowed
  // if perturbations decay relativistically within the horizon

  if (mu2 > 1e-2) {
    yprime[qidx] *= FadeOut;
    yprime[qidx+1] *= FadeOut;
  }

	// STOP CONDITIONS
	//double mQ = quint->q(a) / cosmos->M_p();
	//bool stopCondition = delta_nuNR_longit >= maxDelta; // maxDelta is already normalized to int. units
	//cout << " delta_phi:  " << Dgq*(quint->q(a)/cosmos->M_p())*pow(10,-4.5) << endl;
	//bool stopConditionEta = -beta*Dq*mQ >= maxEta;  
	//if(stopCondition) cout << " Stop Condition " <<endl;
	//if(stopConditionEta) cout << " Stop Condition Eta " << endl; 
		if(stopCondition || stopConditionEta){
		
		// stop field peturbation growth
		//
		//
		yprime[xidx] = 0;

		// stop neutrino growth
		yprime[nuNRidx] = 0;
		// stop neutrino velocity growth
		//yprime[nuNRidx+1]=0;

		// fix DE overdensity
		//yprime[qidx+1] = 0;
	//	Dq = -maxEta / (beta*mQ);
		// set DE velocities to zero
	//	Vq_longit = 0;
		// Stop overdensity growth
		yprime[xidx] = 0;
		//yprime[qidx+1] = 0;
  		//mDelta_q_longit = Dq;

	} 

  mDelta_c_longit = Dgc - 3.*Phi;
  mDelta_b_longit = Dgb - 3.*Phi;
  mV_c_longit = Vc;
  mV_b_longit = Vb;
  mDelta_q_longit = delta_q_longit;
  mDelta_nu_nr_longit = delta_nuNR_longit;
  mV_q_longit = Vq_longit;
  mV_nu_nr_longit = onepwNRVn_nr_longit/(1.+w_nuNR);
  double Dgq = Dq + 3.*(1.+w)*Phi;

//bool printCondition = false;
double kPrint=k;
double kMax=0.180;
double kMin=0.175;
double units=pow(10,-4.5)*pow(2,0.5);
double zz = 1./a; // z + 1

bool printCondition = kPrint >= kMin && kPrint < kMax; 

	if(printCondition){
printOSF("delta_nu", zz, units*mDelta_nu_nr_longit);
printOSF("delta_yNR", zz, units*y[nuNRidx]);
printOSF("delta_yQ", zz, units*y[qidx]);
printOSF("delta_dm", zz, units*mDelta_c_longit);
printOSF("delta_phi", zz, -X/Mp);
printOSF("Phi", zz, Phi);
printOSF("Phi", zz, nuContribToPhi);
//printOSF("delta_q-beta", zz, units*beta*delta_q_longit*mQ);
//printOSF("delta_q_big", zz, units*delta_q_longit*mQ);

}
/*
  PerturbationTracker* t= PerturbationTracker::self(this);
  t->track("n_delta", tau, delta_nuNR_longit);
  t->track("Phi", tau, Phi);
  t->track("PhiDot", tau, PhiDot);
  t->track("q_delta", tau, Dgq);
  t->track("c_delta", tau, mDelta_c_longit);
  t->track("b_delta", tau, mDelta_b_longit);
*/
  }

void CnInvariant::printOSF(string name, double xValue, double yValue) {

name+=".dat";
ofstream osf;
osf.open(name.c_str(), std::ios_base::app);
osf << xValue << "  "  << yValue <<endl;
osf.close();
}

void CnInvariant::printOSF5(string name, double xValue, double yValue, double zValue, double wValue, double kValue) {

name+=".dat";
ofstream osf;
osf.open(name.c_str(), std::ios_base::app);
osf << xValue << "  "  << yValue << "  " <<  zValue << "  " <<  wValue << "  " <<  kValue <<endl;
osf.close();
}

void CnInvariant::initialScalarPerturbations(const ControlPanel &control, const double tau)
{
  // copy the phantomcrossing flag for fderivs()
  IsPhantomCrossing = control.isPhantomCrossing();
	//cout << " ::: CnInvariant::initialScalarPerturbation() ::: " << endl;

//X   SpeedyDEInvariant::initialScalarPerturbations(control, tau);
//X   y[nuNRidx] = y[nuNRidx+1] = 0;
//X   return;

#ifdef OUTPUTFILE
  ofs = new ofstream("colddark.dat");
  ofs2 = new ofstream("colddark2.dat");
#endif

  // unfortunatly, we need to define our M,N,E ourselves, cause
  // they are non -constant (we have to initialize them)
  // see setPointers() for more
  double *M = &y[5];                    // Temperature
  double *E = &M[lmaxg-1];  // E[2] is first after M[lmaxg]

  double *N = 0;
  if (cosmos->nuR()>0) {
    N = &E[lmaxg+1];  //  massless neutrinos
  }

  //  cout << "N[lmaxnr] " << &N[lmaxnr] << "  y[nvar]: " << &y[nvar] << endl;
  double Dgp, Dgc, Dgb, Dgn, Dgq, Vb, Vc, Vp, Vn, Vq, Pi_nu;

  double Vqprime, hq, Q_q_0, S;

  const double a = tau2a(tau);
  const double adot = tau2adot(tau);
  const double adotoa = adot/a;
  static const double Mp = quintcosmos->M_p();
  const double q = quintcosmos->Tau2Phi->fastY(tau);
  const double qdot = quintcosmos->Tau2PhiDot->fastY(tau);
  const double massNu_kbT = quintcosmos->mass_nu_kbT(q);
  const double beta = massNu_kbT?quintcosmos->dMass_nu_kbTdphi(q)/massNu_kbT:0.;
  const double rhoq = quintcosmos->tau2rho_q(tau);
  const double pq = quintcosmos->tau2p_q(tau);
  const double rho_nuNR = quintcosmos->tau2rho_nuNR(tau);
  const double p_nuNR = quintcosmos->tau2p_nuNR(tau);
  const double Ucommaq = quintcosmos->quintessence()->Vprime(q, a, adot/a);

  double x = k*tau;
  double k2 = k*k;

#warning should this be zero??
  double Omega_nu = cosmos->rho2omega(cosmos->tau2rho_nuNR(tau)+cosmos->tau2rho_nu(tau),tau);
  double Omega_gamma = cosmos->rho2omega(cosmos->tau2rho_g(tau),tau);
  double Omega_cdm = cosmos->rho2omega(cosmos->tau2rho_cdm(tau),tau);
  double Omega_b = cosmos->rho2omega(cosmos->tau2rho_b(tau),tau);
  double w=quintcosmos->tau2w_q(tau);
  double w_nuNR=quintcosmos->tau2w_nuNR(tau);
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

    normalize = Psi * (4.0*Omega_nu + 15.0)/10;
    normalize = -1/normalize;
	//cout << "CnInvariant:initialScalar() normalize: " << normalize << endl; 

    break;
  case ControlPanel::isoCDM:
  case ControlPanel::isoNeutrino:
  case ControlPanel::mixed:
  default:
    cout << "CnInvariant::initialScalarPerturbations() - Initial condition not yet implemented" << endl;
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
    for (int l = 3; l <= lmaxnr; l++)
      N[l] =0;
  }

  // fill dark energy. Please note that for DE, we use longitudinal gauge
  double deltaq = Dgq - 3*(1+w)*Phi;
  y[qidx]= deltaq;
  y[qidx+1]=Vq * (1.0 + w);  // we propagate (1+w)*Vq(longitudinal), not Vq

  // field perturbations
  double X = 0;
  double Xprime = 0;
  if (rhoq != 0) {
    double phi = quintcosmos->Tau2Phi->fastY(tau);
    double phiDot = quintcosmos->Tau2PhiDot->fastY(tau);
    double a2 = tau2a(tau)*tau2a(tau);
    // 19.10.10 modified the equation for Xprime
    X = Vq/k*phiDot;
    Xprime = phiDot/(1.+w)*(Dgq-X*Ucommaq/rhoq)+phiDot*(Psi-3.*Phi);
  //Xprime = phiDot*Vqprime/k + (1./k)*(-2*adotoa*phiDot - a2*Ucommaq 
    //     + a2*beta*(rho_nuNR - 3*p_nuNR) )*Vp; 
  //Xprime /= Mp;
  }
  y[xidx] = X;
  y[xidx+1] = Xprime;

  y[nuNRidx] = y[nuNRidx+1] = 0;

  if (cosmos->haveMassiveNeutrinos() && (control.initialCond != ControlPanel::adiabatic)) {
    throw Bad_Error("For massive neutrinos, only adiabatic initial conditions are implemented.");
  } else if (cosmos->haveMassiveNeutrinos()) {
    // massive neutrinos are in longitudinal gauge, too
    double deltan = Dgn - 3.*(1.+cosmos->tau2w_nuNR(tau))*Phi;
    y[nuNRidx  ]  = deltan;
    y[nuNRidx+1]  = Vn * (1.0 + w_nuNR);        // and also here, we're propagating (1+w_nuNR)*Vn_NR

    if (cosmos->nuR()>0) {
     MassiveNeutrinosNew::setFirstIndex((N-y)+lmaxnr+1);
    } else {
     MassiveNeutrinosNew::setFirstIndex((E-y)+lmaxg+1);
    }
   MassiveNeutrinosNew::setIninitalLongitudinalScalarPerturbations(y, cosmos->tau2a(tau), massNu_kbT, deltan, Vn, Pi_nu);
  }

#warning debug
  PerturbationTracker* t = PerturbationTracker::self(this);
  //  t->track("initial-X", tau, X/Mp);
  //t->track("initial-Xprime", tau, Xprime/Mp);
}

void CnInvariant::getReady(const ControlPanel& control)
{
  SpeedyInvariant::getReady(control);  // now call parent objects get ready
  xidx = nvar-5;
  qidx = nvar-3;
  nuNRidx = nvar-1;
}

/*!
  For the sources, to be precise the ISW effect, one needs d/d tau
  of Phi and Psi. In the QuintCnInvariant class, these are needed anyhow
  within the fderivs() function. However, here they are only calculated
  for the sources, which is why we need to do it.

  Please note that fillPhiDot() assumes that it is called after fderivsWrapper() 
  has been called, i.e. Phi and Psi and all other stuff is up to date.
*/
void CnInvariant::fillPhiDot(double tau)
{
  return;
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
 double Pin = 12.0/5.0*N[2];  // neutrinos
 
 // (rho + p)*V of all species and from this and eqn 2.50 d/dtau Phi
 double RhoPV = tau2rho_cdm(tau) * Vc + tau2rho_b(tau)*Vb;
 RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
 RhoPV += quintcosmos->tau2rho_q(tau)* Vq; // simple :-)
 PhiDot  = adotoa*Psi - Gpi4a2*(RhoPV)/k;
 // now d/dtau Psi from d/dtau of eqn 2.52, first d/dtau (p*Pi) p = rho/3
 double Psi2nd = 12.0/5.0*a2*(tau2rhodot_g(tau)*M[2] + tau2rhodot_nu(tau)*N[2] 
			      + tau2rho_g(tau)*Mprime[2] + tau2rho_nu(tau)*Nprime[2])/3.0;
 // add to this d/dtau a^2 term
 Psi2nd += 2.0*a*tau2adot(tau)*(tau2p_g(tau)*Pip + tau2p_nu(tau)*Pin);
 PsiDot = -PhiDot - Gpi8()/k2 * Psi2nd; 
}


void CnInvariant::calcPerturbationNr(const ControlPanel &control) {
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
      nvar +=MassiveNeutrinosNew::pertArraySize();
    }

    nvar += 2; //  2 more for dark energy
    nvar += 2; // and 2 more fore the massive neutrino fluid
    nvar += 2; // for derivatives of field perturbation
  } else nvar = 0;

  if (control.tensor) nvart = 2*lmaxt + 4 + lmaxnr+1; else nvart = 0;

  //cout << "Allocated " << nvar << " entries for k=" << k << endl;

  y = new double[nvar+1];
  yprime = new double[nvar+1];
  yt = new double[nvart+1];
  ytprime = new double[nvart+1];

  std::fill(y, y+nvar+1, 0.);
  std::fill(yprime, yprime+nvar+1, 0.);
  std::fill(yt, yt+nvart+1, 0.);
  std::fill(ytprime, ytprime+nvart+1, 0.);
}


