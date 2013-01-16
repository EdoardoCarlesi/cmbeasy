#include "invariant.h"
#include "cosmos.h"

// uncomment this, if you would like to write to a file during fderivs
//#define OUTPUTFILE
 
Invariant::Invariant(Cosmos *cosmos) : Perturbation(cosmos),  lmaxnu(0), s5(sqrt(5.0)), s6(sqrt(6.0)) {
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
}

/*!
  This are the coefficients kappa with indices s, l and m.  Hence
  the rather strange name.
*/
double Invariant::sklm(double s, double l, double m) {
  return sqrt( (l*l - m*m)*(l*l - s*s) / (l*l));
}

/*!
  tau-derivatives of Gauge-Invariant Variables. y and yprime come
  from odeint(), meaning that they usually will not coincide with our
  objects y and yprime. However, fderivsWrapper() will at the end
  of the integration fill y and yprime of our object with the result. 
  In other words: scalarSources() etc will have a fully valid y and 
  yprime arrary and when scalarSources() calls fderivs(), it does so
  with our objects y and yprime
*/
void Invariant::fderivs(const double tau, const double *y, double *yprime) {
  // initialize E,N,M ... etc makes bookkeeping of Multipoles quite easy 
  // note:
  // y and yprime come from Miscmath::odeint() and are not the same as our
  // object's y and yprime (our object's y will however later contain 
  // the result of the odeint calculation)

  cout << "this is fderivs!\n";

  setPointers(y,yprime); 
  double Dgb = y[1];
  double Vb = y[2];
  double Dgc = y[3];
  double Vc = y[4]; 
 
  // Dg and V for photons
  double Dgp = 4*M[0];
  double Vp = M[1];
  double Dgn = 4*N[0];
  double Vn = N[1];
  if (inform)  cout << "dgp, vp: " << Dgp << "  " << Vp << endl;
  double a= tau2a(tau);   // we take it from cosmos :-), as this is also mother, children will be fast
  double a2 = a * a;
  double adotoa = tau2adot(tau)/a;

  double Gpi4 = Gpi8()*0.5;
  double Gpi4a2 = Gpi4*a2;
 
  // first, we  calculate Phi
  double a2ak = adotoa / k;
  double dgrho = tau2rho_cdm(tau) *(Dgc + 3 *a2ak*Vc);  //we collect all species, first: cdm
  dgrho += tau2rho_b(tau) * (Dgb + 3*a2ak*Vb);              //baryons
  dgrho += tau2rho_g(tau) * (Dgp + 4*a2ak*Vp);              //photons
  dgrho += tau2rho_nu(tau) * (Dgn + 4*a2ak*Vn);           // massless neutrinos
  dgrho *= Gpi4a2;
  double phicontrib = 3*(tau2rho_b(tau) + tau2rho_cdm(tau)) + 4*(tau2rho_g(tau) + tau2rho_nu(tau));
  phicontrib *= Gpi4a2;
  Phi = dgrho / (k*k + phicontrib);  
  

  // then we collect all the shear sources there are
  double Pip = 12.0/5.0*M[2];  // photons
  double Pin = 12.0/5.0*N[2];   // neutrinos

  // then we get Psi, first we get the total shear (times perturbation)
  double dgpres = cosmos->tau2p_g(tau)*Pip + cosmos->tau2p_nu(tau)*Pin;
  dgpres *= 2.0*Gpi4a2;

  Psi = -Phi - dgpres / k2;   
   

  double R = cosmos->tau2R(tau);  // 4/3 * rho_photon / rho_baryon
  double cs2 = soundSpeed(tau);  // baryon soundspeed (almost 0, for tau > 1e-2) 
  // d/dtau for Dgb
  yprime[1] = -k*Vb - 3*cs2*adotoa*Dgb ; 
  // and d/dtau for M[0], which is 1/4 of Dgp 		  
  Mprime[0] = -k/3.0 * Vp;

  //now we calculate d/dtau of Dg and V for cdm
  yprime[3] = -k*Vc;
  yprime[4] =  k*Psi - adotoa * Vc;

  bool coupled = isTightCoupling(tau,photon);

  double Vb_prime; // d/dtau baryon velocity
  double Vp_prime=-1e100; // d/dtau  photon velocity (initialize to silence compiler warning) 
  if (isTightCoupling(tau,photon)) {
    //Vb_prime = ( adotoa *Vb*(3*cs2 - 1) + k*R*(M[0]) + (R+1)*k*Psi - k*Phi*(R + 3*cs2) + k*cs2*Dgb) / (1+R);
    Vb_prime = (  -adotoa *Vb + k*R*(M[0]) + (R+1)*k*Psi - k*Phi*(R + 3*cs2) + k*cs2*Dgb) / (1+R);
    // the slip is d/dtau ( V_b - V_photon)  the gauge-invariant formulation of eqn 74 in Ma & Bertschinger
    // in principle we need two additional things, (d^2/dtau^2 a)/a  **and** d/dtau phi
     // (rho + p)*V of all species and from this  d/dtau Phi
    double RhoPV = tau2rho_cdm(tau) * Vc + tau2rho_b(tau)*Vb;
    RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
    double PhiDot  = adotoa*Psi - Gpi4a2*(RhoPV)/k;
    double adotdota = (adotoa * adotoa -  Gpi8()*a2*tau2p(tau)) * .5;
    double slip = adotoa*2*R*(Vb - Vp); 
    slip += ( -adotdota*Vb - adotoa*k*(0.5*Dgp + Psi - 2*Phi) + k*(cs2*yprime[1] - Mprime[0] + PhiDot*(1-3*cs2) ) )/opac(tau);
    slip /= (1+R);
    Vb_prime += R/(1+R)*slip; 
    
#ifdef OUTPUTFILE
    if (tau < 1) 
      (*ofs) << tau <<  "  " << Vb - Vp << " " << k*1/( opac(tau)*(1+R) *tau*tau*tau) << "  " << 1/( opac(tau)*(1+R))*adotoa*k<<  " " <<  adotoa*k*(0.5*Dgp + Psi - 2*Phi)<< "  " << k*(cs2*yprime[1] - Mprime[0] + PhiDot*(1-3*cs2) ) <<  " " << adotoa << " " << Dgp << endl;
#endif
    
    
  } else {
    Vb_prime  = k*(Psi -3*cs2*Phi) + Dgb*k*cs2 - (adotoa - 3*cs2)*Vb + opac(tau)*R*(Vp - Vb);
  }
  

  yprime[2] = Vb_prime;
  
  if (coupled) {
    // for the photons, we solve the equivalent of eqn 70 Ma + Bertschinger for Vp_prime
    Vp_prime = (- Vb_prime - adotoa*Vb + k*(cs2*Dgb - 3*cs2*Phi + (R+1)*Psi) )/R + k*(M[0]-2.0/5.0*M[2] - Phi); 
    for (int l = 2; l <= lmaxg; l++) {
      Mprime[l] = 0;
      Eprime[l]  = 0;
    }
  } else {
    //Vp_prime = -(Vb_prime + adotoa*Vb - k*cs2*Dgb + 3*cs2*k*Phi)/R + k*(M[0]-2.0/5.0*M[2]) - k*Phi + (R+1)/R * k*Psi; 
    Vp_prime = opac(tau)*(Vb - Vp) + k*(Psi - Phi) + k *(M[0]-2.0/5.0*M[2]);
    Mprime[2] = -opac(tau)*(9.0/10.0 * M[2] + s6/10.0*E[2]) + k*(2.0/3.0 * Vp - 3.0/7.0*M[3]);
    for (int l = 3; l < lmaxg; l++) 
      Mprime[l]  = - opac(tau) * M[l] + k*( denom1[l] * M[l-1] - denom2[l]*M[l+1]);
    // truncation of the scheme similar to  Ma & Bertschinger 
    Mprime[lmaxg] =  ((double) 2*lmaxg +1)/((double) 2*lmaxg-1)*k*M[lmaxg - 1] - M[lmaxg]*( (lmaxg+1)/tau + opac(tau));
 
    // propagation of E - Polarization for photons
    Eprime[2] = -k*( s5 / 7.0 *E[3] ) - opac(tau) * (E[2] + s6/10.0 * (M[2] - s6*E[2]));
    for (int l = 3; l < lmaxg; l++) 
      Eprime[l]  = k*( denomE1[l] * E[l-1] - denomE2[l]*E[l+1] ) - opac(tau)*E[l];
    // truncation of the scheme exactly like in Ma & Bertschinger 
    Eprime[lmaxg] =  ((double) 2*lmaxg +1)/((double) 2*lmaxg-1)*k*E[lmaxg - 1] - E[lmaxg]*( (lmaxg+1)/tau + opac(tau));
  }    
  Mprime[1] = Vp_prime;

  // massless neutrinos
  Nprime[0] = -k/3.0 * N[1];
  Nprime[1] = k*(N[0] - 2.0/5.0*N[2] + Psi - Phi);
  for (int l = 2; l < lmaxnr; l++)
    Nprime[l] = k*(denom1[l]*N[l-1] - denom2[l]*N[l+1]);
  // truncation just like photons
  Nprime[lmaxnr] = ((double) 2*lmaxnr +1)/((double) 2*lmaxnr-1) *k*N[lmaxnr - 1] - N[lmaxnr]* (lmaxnr+1)/tau;
  
} 

void Invariant::initialScalarPerturbations(const ControlPanel &control, const double tau)
{

  cout << "initialScalarPerturbations\n";

#ifdef OUTPUTFILE  
  ofs = new ofstream("colddark.dat");
#endif

  // unfortunatly, we need to define our M,N,E ourselves, cause
  // they are non -constant (we have to initialize them)
  // see setPointers() for more
  double *M = &y[5];                    // Temperature
  double *E = &M[lmaxg-1];  // E[2] is first after M[lmaxg]
  double *N = &E[lmaxg+1];  //  massles neutrinos

  double Dgp, Dgc, Dgb, Dgn, Vb,Vc,Vp,Vn,Psi,Pi_nu;
  k2 = k*k;

  double x = k*tau;

  double Omega_nu = cosmos->rho2omega(cosmos->tau2rho_nu(tau),tau);
  double Omega_gamma = cosmos->rho2omega(cosmos->tau2rho_g(tau),tau);
  double Omega_cdm = cosmos->rho2omega(cosmos->tau2rho_cdm(tau),tau);
  double Omega_b = cosmos->rho2omega(cosmos->tau2rho_b(tau),tau);

   double Q = 1.0/(3.0*Omega_cdm+8.0*(Omega_nu + Omega_gamma));

  double P = 1.0/(4.0*Omega_nu + 15.0);
  double T=1/4;
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
    Psi= -Phi-Omega_nu*Gpi8()*Pi_nu/(x*x);
    normalize = -1.0/Psi;  // this is for tensor / scalar ratio , only for adiabatic see below
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
    throw Bad_Error("Invariant::initialScalarPerturbations():  This initial condition is notsupported.");
  }

  
  // Mmmhh. I know, everyone will hate me for this, but:
  // The initial conditions may be nice, yet I didn't derive them
  // under the condition that Psi = -1
  // However, in order to ensure that the Tensor modes (which
  // are ripped from the synchronous original cmbfast) have the
  // right amplitude w.r.t. the scalars, I scale here such that
  // Psi (which is actually not needed, as it is computed within
  // fderivs() ) equals -1 

  // this is only in effect for adiabatic conditions, since else normalize == 1

  Dgp *= normalize;
  Dgn *= normalize;
  Dgb *= normalize;
  Dgc *= normalize;
  Vp *= normalize;
  Vb *= normalize;
  Vn *= normalize;
  Vc *= normalize;
  Pi_nu *= normalize;
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
} 

void Invariant::initialTensorPerturbations() {
  yt[1] = 1;   
  yt[2] = 0.;
  int ind1 = 3;
  int ind2 = ind1 + lmaxt +1;
  for (int l = 0; l <= lmaxt; ++l) {
    yt[ind1 + l] = 0.;
    yt[ind2 + l] = 0.;
  }
} 



#define LM2 7 // 7
#define LM3 4
#define NQ1 15 

void Invariant::getReady(const ControlPanel& control) {
  if (control.scalar) {
    lmaxg = LMAX0;
    if (control.highPrecisionTransfer) lmaxnr = LMAXNR0; else lmaxnr = LM2;
    
  }
  
  if (control.tensor) lmaxt = LMAXT0; else lmaxt = 0;
  
  /*     Initialize neutrino mass and storage. */
  if (nuNR() == 0. || cosmos->omega_nu() == 0.) {
    nqmax = 0;
    lmaxnu = 0;
  } else {
    if (control.highPrecisionTransfer) {
      nqmax = NQMAX0;
      lmaxnu = LMAXNU0;
    } else {
      nqmax = NQ1;
      lmaxnu = LM3;
    }
  }

  // Calculate qdn 
  double dq = 1.,q;
  for (int i = 0; i <  NQMAX0; ++i) {
    q = i + .5;
    // 	  aq=a*amnu/q 
    // 	  v=1.0d0/sqrt(1.0d0+aq*aq) 
    qdn[i] = dq * q * q * q / (exp(q) + 1.);
  }
  epsw = 100. / cosmos->tau_0();
  Perturbation::getReady(control);  // now call parent objects get ready
}

/*!
  For the sources, to be precise the ISW effect, one needs d/d tau
  of Phi and Psi. In the QuintInvariant class, these are needed anyhow
  within the fderivs() function. However, here they are only calculated
  for the sources, which is why we need to do it.

  Please note that fillPhiDot() assumes that it is called after fderivsWrapper() 
  has been called, i.e. Phi and Psi and all other stuff is up to date.
*/
void Invariant::fillPhiDot(double tau) {
  double Vb = y[2];
  double Vc = y[4]; 
  double Vp = M[1];
  double Vn = N[1];
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
  PhiDot  = adotoa*Psi - Gpi4a2*(RhoPV)/k;
  // now d/dtau Psi from d/dtau of eqn 2.52, first d/dtau (p*Pi) p = rho/3
  double Psi2nd = 12.0/5.0*a2*(tau2rhodot_g(tau)*M[2] + tau2rhodot_nu(tau)*N[2] 
			       + tau2rho_g(tau)*Mprime[2] + tau2rho_nu(tau)*Nprime[2])/3.0;
  // add to this d/dtau a^2 term
  Psi2nd += 2.0*a*tau2adot(tau)*(tau2p_g(tau)*Pip + tau2p_nu(tau)*Pin);
  PsiDot = -PhiDot - Gpi8()/(k*k) * Psi2nd;
}


void Invariant::scalarSources(double tau, double *d, double *dp, double *dk) {  

  cout << "das ist scalarSources!\n";

  fderivsWrapper(tau,y,yprime);  // get Phi and all derivatives of M[l] etc
  fillPhiDot(tau);   // get PhiDot and PsiDot right
 
  double Dgp = 4*M[0];
  double Vb = y[2];
  
 
  double C = (M[2] - s6*E[2])/10.0; 
  double Cdot = (Mprime[2] - s6*Eprime[2])/10.0;

  double M2ddot = -dopac(tau)*(9*M[2] + s6*E[2])/10.0 - opac(tau)*(9*Mprime[2] + s6*Eprime[2])/10.0;
  M2ddot += k*(2.0/3.0*Mprime[1] - 3.0/7.0*Mprime[3]);
  double E2ddot = -k*s5/7.0 *Eprime[3] - dopac(tau)*(E[2] + s6*C) - opac(tau)*(Eprime[2] + s6*Cdot);
  double Cddot =  (M2ddot - s6*E2ddot)/10.0;

 
   double Vbdot = yprime[2];

   *d = 0; // Mona  

   // double delta_gamma = Dgp - 4*Phi; //  Michael: das ist in synchroner Eichung -- hoffentlich ;-) -> should use that to see the different contributions to cmb (because other people used that, too). The contributions get redistributed when changing the gauge...


   *d += -expmmu(tau)*(PhiDot - PsiDot);  // ISW
   *d += visibility(tau)*( -Phi + Psi  + 0.25*Dgp + Vbdot/k + C/2.0 + 3.0/(2.0*k*k)*Cddot );
   *d += dvisibility(tau)*(Vb / k + 3.0*Cdot/(k*k));
   *d += ddvisibility(tau)*3.0/(2.0*k*k)*C;   

   // *d =  visibility(tau)*0.25*delta_gamma; // should be the acoustic contribution  

   // Mona start: give out sources:
   double isw_mona = -expmmu(tau)*(PhiDot - PsiDot);
   double acoustic_mona =  visibility(tau)*( -Phi + Psi  + 0.25*Dgp + Vbdot/k + C/2.0 + 3.0/(2.0*k*k)*Cddot );
   double doppler_mona = dvisibility(tau)*(Vb / k + 3.0*Cdot/(k*k));
   double fourth_mona = ddvisibility(tau)*3.0/(2.0*k*k)*C;   

   //  this does not work -- seems to be working only for fderifs!!!!
   /*
#ifdef OUTPUTFILE
   (*ofs) << tau << "  " << isw_mona << "  " << acoustic_mona << "  " << doppler_mona << "  " << fourth_mona << "  " << acoustic_mona + doppler_mona << "  " << isw_mona + acoustic_mona + doppler_mona << "  " << isw_mona + acoustic_mona + doppler_mona + fourth_mona << endl;
#endif
   */
   

  
   // Mona end

   double x = k * (tau_0() - tau);
   if (x > 0.) *dp = visibility(tau) * 3. / 2. * C / (x * x); else *dp = 0.;
  
} 

void Invariant::fderivsTensor(const double tau,   double *y, double *yprime) {
    static double  psie, htpr, a;
    static int l;
    static double a2, htdpr, rhonu, ep, ht, adotoa;
  
    static double cs2;
    static double deltap0, deltat0, tcp, pnu;
    static int ind1, ind2;
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
    if (cosmos->nuNR() == 0.) {
      rhonu = 1.;
      pnu = .33333333333333331;
    } else {
      throw Bad_Error("fix massive nu");
//X       cosmos->nu1(a, &rhonu, &pnu);
    }

    //  8*pi*G*rho*a**2 and 8*pi*G*P*a**2. 
    // Tensors 

    adotoa = tau2adot(tau)/a;
    ht = y[1];
    htpr = y[2];
    yprime[1] = htpr;
    htdpr = adotoa * -2 * htpr - k2 * ht;
    yprime[2] = htdpr;

    // Photon perturbations 
    ind1 = 3;
    ind2 = ind1 + lmaxt + 1;
    psie = y[ind1] / 10. + y[ind1 + 2] / 7. + y[ind1 + 4] * 3. / 70. - y[ind2]
	     * 3. / 5. + y[ind2 + 2] * 6. / 7. - y[ind2 + 4] * 3. / 70.;
 
    if (tcp == 1.) {
      // no tight coupling approx 
	yprime[ind1] = -k * y[ind1 + 1] - opac(tau) * y[ind1] + opac(tau) * psie - htpr;
 	yprime[ind2] = -k * y[ind2 + 1] - opac(tau) * y[ind2] - opac(tau) * psie;
	// l=1...lmaxt 
	for (l = 1; l < lmaxt ; ++l) {       //! Comment above no true: LMAXT IS NOT REACHED!!!
	    yprime[ind1 + l] = k * denl[l] * (l * y[
		    ind1 - 1 + l] - (l + 1) * y[ind1 + 1 + l]) - opac(tau) * y[ind1 + l];
	    yprime[ind2 + l] = k * denl[l] * (l * y[
		    ind2 - 1 + l] - (l + 1) * y[ind2 + 1 + l]) - opac(tau) * y[ind2 + l];
	}

	// Truncate moment expansion  , now LMAXT is REACHED
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
}  

void Invariant::tensorSources(double tau, double *dt, double *dte, double *dtb) {
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
    *dt =  (- expmmu(tau) * htpr + visibility(tau)* psie) / x2;  
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




// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
void Invariant::nu2(const double a, double *drhonu, double *fnu, double *dpnu, double *shearnu, const double *psi0, const double * psi1, const double *psi2) {
  throw Bad_Error("fix massive nu");
//X   //cout << "nu2" << endl;
//X     static double q, v, g0[4], g1[NQMAX0+2], g2[NQMAX0+2], g3[NQMAX0+2], g4[NQMAX0+2], aq;
//X     static int iq;
//X     static double gf1, gf2, gf3, gf4;
//X    
//X     //  Compute the perturbations of density, energy flux, pressure, and 
//X     //  shear stress of one flavor of massive neutrinos, in units of the mean 
//X     //
//X     //  density of one flavor of massless neutrinos, by integrating over 
//X     //  momentum. 
//X 
//X     //  const=7*pi**4/120.  = 5.68219698
//X 
//X     if (nqmax == 0) {
//X       *drhonu = 0.;
//X       *fnu = 0.;
//X       *dpnu = 0.;
//X       *shearnu = 0.;
//X       return;
//X     }
//X 
//X     //  q is the comoving momentum in units of k_B*T_nu0/c. 
//X     g1[0] = 0.;
//X     g2[0] = 0.;
//X     g3[0] = 0.;
//X     g4[0] = 0.;
//X     for (iq = 2; iq <= NQMAX0+1; ++iq) {
//X       q = iq -  1.5;
//X       aq = a * cosmos->amnu / q;
//X       v = 1. / sqrt(aq * aq + 1.);
//X       g1[iq-1] = qdn[iq - 2] * psi0[iq-2] / v;
//X       g2[iq-1] = qdn[iq - 2] * psi0[iq-2] * v;
//X       g3[iq-1] = qdn[iq - 2] * psi1[iq-2];
//X       g4[iq-1] = qdn[iq - 2] * psi2[iq-2] * v; 
//X     }
//X     Miscmath::splint(g1, g0, NQMAX0+1);
//X     Miscmath::splint(g2, &g0[1], NQMAX0+1);
//X     Miscmath::splint(g3, &g0[2], NQMAX0+1);
//X     Miscmath::splint(g4, &g0[3], NQMAX0+1);
//X     gf1 = g1[NQMAX0];
//X     gf2 = g2[NQMAX0];
//X     gf3 = g3[NQMAX0];
//X     gf4 = g4[NQMAX0];
//X     *drhonu = (g0[0] + gf1 * 2. / 14.5) / 5.68219698;
//X     *fnu = (g0[2] + gf3 * 2. / 14.5) / 5.68219698;
//X     *dpnu = (g0[1] + gf2 * 2. / 14.5) / 5.68219698 / 3.;
//X     *shearnu = (g0[3] + gf4 * 2. / 14.5) / 5.68219698 * 2. / 3.;
//X     //cout << a << " pnu: " << *dpnu << "    " << *drhonu << "   "  << *shearnu << endl;
} // nu2_ 

double Invariant::delta_nr(double tau) {
  throw Bad_Error("fix massive nu");
//X   double a = tau2a(tau);
//X   double rhonu,pnu,shearnu,fnu,drhonu,dpnu;
//X   cosmos->nu1(a,&rhonu,&pnu);
//X   nu2(a,&drhonu,&fnu,&dpnu,&shearnu,&y[iq0],&y[iq1],&y[iq2]);
//X   return drhonu/rhonu/(k*k);
}


void Invariant::calcPerturbationNr(const ControlPanel &control) {
  //   Calculate number of equations 
  if (control.scalar) {
    iq0 = (lmaxg << 1) + 9 + lmaxnr;
    iq1 = iq0 + nqmax;
    iq2 = iq1 + nqmax;

    nvar = 4;  // CDM and Baryon Delta and Velocity
    nvar += lmaxg +1; // Photons  0.. lmaxg
    nvar += lmaxg -1; // Polarization 2 .. lmaxg
    nvar += lmaxnr +1; // massless neutrinos 
  } else nvar = 0;

  if (control.tensor) nvart = (lmaxt << 1) + 4; else nvart = 0;

  y = new double[nvar+1];
  yprime = new double[nvar+1];
  yt = new double[nvart+1];
  ytprime = new double[nvart+1];
  
}

/*!
  Set M, E, N, Mprime, Eprime, Nprime such that 
  they have the right position within y and yprime. 
  This just centralises a bit the bookkeeping for the
  fderivs() and initialScalarPerturbation()
*/
void Invariant::setPointers(const double *y, double *yprime) {
  // set pointers M, E, N for Temperature, Polarization and 
  // massless neutrinos. Just convenience, but diff eqn. much more
  // readable with this
  M = &y[5];                    // Temperature
  // Polarization, E[2] is the first one we need and we 
  // therefore readjust the beginning such that E[2] points to the first
  // free variable after the M's
  E = &M[lmaxg-1];  // E[2] is first after M[lmaxg]
  N = &E[lmaxg+1]; // N[0] is firstafter E[lmaxg]
  Mprime = &yprime[5];
  Eprime = &Mprime[lmaxg-1];
  Nprime = &Eprime[lmaxg+1];
}

void Invariant::nuder(const double a, const double adotoa, const double rhonu, double *rhonudot, double *shearnudot, const double *psi2, const double *psi2dot) {
  throw Bad_Error("fix massive nu");
//X 
//X   // Local variables 
//X     static double vdot, d;
//X     //static int i;
//X     double q, v, aqdot, g0, g1[NQMAX0+1], aq;
//X     //static int iq;
//X     static double gf1;
//X     Nu1d_ &nu1d = cosmos->nu1d;
//X    
//X     //  Compute the time derivative of the mean density in massive neutrinos 
//X     //
//X     //  and the shear perturbation. 
//X 
//X 
//X     // Parameter adjustments 
//X     --psi2dot;
//X     --psi2;
//X 
//X     if (nqmax == 0) {
//X       *rhonudot = 0.;
//X       *shearnudot = 0.;
//X       return;
//X     }
//X 
//X     //  q is the comoving momentum in units of k_B*T_nu0/c. 
//X     g1[0] = 0.;
//X     for (int iq = 1; iq <= NQMAX0; ++iq) {
//X       q = iq - 1.5;
//X       aq = a * cosmos->amnu / q;
//X       aqdot = aq * adotoa;
//X       v = 1. / sqrt(aq * aq + 1.);
//X       vdot = -aq * aqdot / pow(aq * aq + 1. , 1.5);
//X       g1[iq] = qdn[iq - 1] * (psi2dot[iq] * v + psi2[iq] * vdot);
//X     }
//X     Miscmath::splint(g1, &g0, NQMAX0+1);
//X     gf1 = g1[NQMAX0];
//X     *shearnudot = (g0 + gf1 * 2. / 14.5) / 5.68219698 * 2. / 3.;
//X 
//X     d = log(a / nu1d.amin) / nu1d.dlna + 1.;
//X     int i = (int) d;
//X     d -= i;
//X     if (i < 1) {
//X       //  Use linear interpolation 
//X       *rhonudot = nu1d.dr1[0] + (d - 1) * nu1d.ddr1[0];
//X     } else if (i > NRHOPN) {
//X       //  This should not happen, unless the user evolves to z<0! 
//X       *rhonudot = nu1d.dr1[NRHOPN-1] + (d + i - NRHOPN) * nu1d.ddr1[NRHOPN-1];
//X     } else {
//X       //  Cubic spline interpolation for rhonudot. 
//X       *rhonudot = nu1d.dr1[i - 1] + d * (nu1d.ddr1[i - 1] + d * ((
//X 		nu1d.dr1[i] - nu1d.dr1[i - 1]) * 3. - nu1d.ddr1[i - 1] *
//X 		 2. - nu1d.ddr1[i] + d * (nu1d.ddr1[i - 1] + nu1d.ddr1[
//X 		i] + (nu1d.dr1[i - 1] - nu1d.dr1[i]) * 2.)));
//X     }
//X     *rhonudot = rhonu * adotoa * *rhonudot / nu1d.dlna;
} // nuder_  

