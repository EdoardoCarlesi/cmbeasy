//class CoupledInvariant implementation
// version 0.9 (initial conditions are not yet right)
#include "coupledinvariant.h"
#include "coupledquintcosmos.h"
#include "originalcoupledquintcosmos.h"
#include <exception>
#include <iomanip>

// uncomment this, if you would like to write to a file during fderivs
//#define OUTPUTFILE

//Georg: Christian and Mona forgot to add 2 for the Nr of Perturbations in calcPerturbationNr ->diff of about 1.5% at l=1400
bool CoupledInvariant::UseWrongNumberOfPertEquationLikeOriginal = true;

CoupledInvariant::CoupledInvariant(QuintCosmos *cosmos) : Invariant(cosmos) {
    quintcosmos = dynamic_cast<OriginalCoupledQuintCosmos*>(cosmos);
    if (!quintcosmos)
        throw Bad_Error("CoupledInvariant must be constructed with an \
                         OriginalCoupledQuintCosmos as argument given to the constructor.");
}

void CoupledInvariant::fderivs(const double tau, const double *y, double *yprime) {
  //inform = true;
  if (inform) cout << endl << "invarinat::fderivs()   " << tau << endl;

  if(quintcosmos->bprimespline==0) cout << "Error in CoupledInvariant::fderivs: bprimespline not set!" << endl;
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

  // QUINTESSENCE
  // Update quintessence - field so it can later tell us about energy fluctuations etc...
  quint->touch(tau);

  // Dg and V for photons
  double Dgp = 4*M[0];
  double Vp = M[1];

  // Dg and V for neutrinos
  double Dgn = 4*N[0];
  double Vn = N[1];
  if (inform)  cout << "dgp, vp: " << Dgp << "  " << Vp << endl;

  double a= tau2a(tau);   // we take it from cosmos :-), as this is also mother, children will be fast
  double adot= tau2adot(tau);
  double a2 = a * a;
  double adotoa = tau2adot(tau)/a;

  double Gpi4 = Gpi8()*0.5;
  double Gpi4a2 = Gpi4*a2;
 
  // Dg and V for quintessence
  double Dgq = y[qidx];
  double Vq = y[qidx+1];
  
  double q = quint->q(a);
  double qdot = quint->qDot(a);  // d/dtau of background field value
 
  double a2ak = adotoa / k;
  // shear times pressure needed
   // then we collect all the shear sources there are
  double Pip = 12.0/5.0*M[2];  // photons
  double Pin = 12.0/5.0*N[2];   // neutrinos

  double gpres = cosmos->tau2p_g(tau)*Pip + cosmos->tau2p_nu(tau)*Pin;
  gpres *= 2.0*Gpi4a2;

  // first, we  calculate Phi eqn. 2.50
  // as mentioned in Appendix of PhD thesis,
  // Phi  =  A / ( k^2 + B) 
  // and dgrho corresponds to A, phicontrib corresponds to B
  double dgrho = 0;
  dgrho += tau2rho_cdm(tau) *(Dgc);  //we collect all species, first: cdm // Mona: 1a
  dgrho += tau2rho_cdm(tau) *(3 *a2ak*Vc); // Mona: 1b
  dgrho += tau2rho_b(tau) * (Dgb + 3*a2ak*Vb);    //baryons // Mona: 2
  dgrho += tau2rho_g(tau) * (Dgp + 4*a2ak*Vp);    //photons // Mona: 3
  dgrho += tau2rho_nu(tau) * (Dgn + 4*a2ak*Vn);   // massless neutrinos // mona:4 
  dgrho += quintcosmos->tau2rho_q(tau) * (Dgq + 3.0*(1.0+quintcosmos->tau2w_q(tau))*a2ak*Vq); //quintessence // Mona: 5
 
  double phicontrib = 0;
  phicontrib += 3.0*tau2rho_b(tau) + 3.0*tau2rho_cdm(tau)+4.0*tau2rho_g(tau)+4.0*tau2rho_nu(tau); // Mona: a
  phicontrib += 3.0*(1.0+quintcosmos->tau2w_q(tau))*quintcosmos->tau2rho_q(tau); // Mona: b

  //  double dgrhobary =  tau2rho_b(tau) * (Dgb + 3*a2ak*Vb);   // VM 
  //  double phicontribbary =  3.0*tau2rho_b(tau);                // VM
  //  double Phibaryon =  dgrhobary / (k*k/Gpi4a2 + phicontribbary);     // VM

  Phi = dgrho / (k*k/Gpi4a2 + phicontrib);  // note: Gpi4a2 = a^2/(2 Mpl^2)
  Psi = -Phi - gpres / k2;   // eqn. 2.52                               
//X   cout << setprecision(17) << " " << Phi << " " << Psi << " " << gpres << " " << " " << Pin << " " << Pip<< " " << M2 << M[2] << endl;
//X   exit(-1);
  
  // VM start  
  
  double vmVb = y[2];
  double vmVc = y[4]; 
  double vmVp = M[1]; 
  double vmVn = N[1];
  double vmVq = y[qidx+1];
  double vma= tau2a(tau);   // we take it from cosmos :-), as this is also mother, children will be fast
  double vma2 = vma * vma;
  double vmadotoa = tau2adot(tau)/vma;

  double vmGpi4 = Gpi8()*0.5;
  double vmGpi4a2 = vmGpi4*vma2;

  // double vmPip = 12.0/5.0*M[2];  // photons
  // double vmPin = 12.0/5.0*N[2];   // neutrinos

   // (rho + p)*V of all species and from this and eqn 2.50 d/dtau Phi
  double vmRhoPV = tau2rho_cdm(tau) * vmVc + tau2rho_b(tau)*vmVb; 
  vmRhoPV += (1+quintcosmos->tau2w_q(tau))*quintcosmos->tau2rho_q(tau)*vmVq;
  vmRhoPV += 4.0/3.0*(tau2rho_g(tau)*vmVp + tau2rho_nu(tau)*vmVn);
  double vmPhiDot  = vmadotoa*Psi - vmGpi4a2*(vmRhoPV)/k;

  // Mona: get Psidot 
  // now d/dtau Psi from d/dtau of eqn 2.52, first d/dtau (p*Pi) p = rho/3
  double vmPsi2nd = 12.0/5.0*a2*(tau2rhodot_g(tau)*M[2] + tau2rhodot_nu(tau)*N[2] 
			       + tau2rho_g(tau)*Mprime[2] + tau2rho_nu(tau)*Nprime[2])/3.0;
  // add to this d/dtau a^2 term
  vmPsi2nd += 2.0*a*tau2adot(tau)*(tau2p_g(tau)*Pip + tau2p_nu(tau)*Pin);
  double vmPsiDot = -vmPhiDot - Gpi8()/(k*k) * vmPsi2nd;


  // VM end

  //#ifdef OUTPUTFILE   // VM
  //  (*ofs) << tau << " " << a << " " << Dgq << "  " << Vq << " " << Dgq - 3*(1+quintcosmos->tau2w_q(tau))*Phi << endl;   // VM
  //#endif      // VM

  //#ifdef OUTPUTFILE    // VM
  //  (*ofs) << tau << " " << k << " " << Dgq << "  " << Vq << " " << Dgc << "  " << Vc << " " << Phi << "  " << Phibaryon << endl;   // VM
  //#endif   // VM

  


  // introduce some quantities of use for nonzero coupling (from Kodama & Sasaki modified to include field-
  // dependent couplings)
  double Mpl = cosmos->M_p();
  double wq= quintcosmos->tau2w_q(tau);
  double beta = quintcosmos->getB(quint->q(a));  //coupling
  double betaphi = quintcosmos->getBprime(quint->q(a)); // first derivative of the coupling w.r.t. the field
  double betaphiphi = quintcosmos->getBprimeprime(quint->q(a)); // second derivative
  double S = -beta/Mpl*tau2rho_cdm(tau);
  double Q_M = S/a*qdot;
  double Q_q = - Q_M;
  double Sdot = - betaphi*qdot/Mpl*tau2rho_cdm(tau) - beta/Mpl*(-3.0*adotoa*tau2rho_cdm(tau)-beta/Mpl*tau2rho_cdm(tau)*qdot);
  double Qdot_M = Sdot/a*qdot-3.0*adotoa/a*S*qdot-S*a*quint->Vprime(quint->q(a),a,adotoa)-S*S*a;
  double q_M = S*a*qdot/(3*adot*tau2rho_cdm(tau));
  double h_q = 1/(a*a)*qdot*qdot;
  double q_q = - S*a*qdot/(3*adot*h_q);
  // cout << beta <<  " " << betaphi*q << " " << betaphiphi*q << endl; 
  //  double DS = -tau2rho_cdm(tau)/Mpl*(beta+betaphi*q)*(Dgc-3*(1-q_M)*Phi)-1/Mpl*tau2rho_cdm(tau)*
  //  (2*betaphi+betaphiphi*q)*qdot*Vq/k;
  double DS = -tau2rho_cdm(tau)/Mpl*beta*(Dgc-3*(1-q_M)*Phi)-1/Mpl*tau2rho_cdm(tau)*betaphi*qdot*Vq/k;
  double Xdot = qdot*(-3*(1-q_q)*Phi+Psi+Dgq/(1+wq)-quint->Vprime(quint->q(a),a,adotoa)*
		      qdot/(k*quintcosmos->tau2rho_q(tau)*(1+wq))*Vq); 
  double E_M = (1/a*S*Xdot+1/a*qdot*DS+Q_q*Psi-1/k*Qdot_M*Vc); 
  // This has been multiplied with Q_M from the original reference
  double E_SM= E_M + Qdot_M/(k)*Vc;  // to cope with divergences if Q_M =0
  double A = a*a/(2*k*adotoa*Mpl*Mpl)*(tau2rho_b(tau)*Vb+tau2rho_cdm(tau)*Vc+4/3*tau2rho_g(tau)*Vp+
	              4/3*tau2rho_nu(tau)*Vn+(1+wq)*quintcosmos->tau2rho_q(tau)*Vq)-a*a/(2*adotoa*adotoa*Mpl*Mpl)*Phi*phicontrib/3.0; // Mona: Vn forgotten in  4/3*tau2rho_nu(tau)*Vn, I've put it there
//X   cout << setprecision(17) << " " << DS << " " << Xdot << " "<< E_M << " " << E_SM << " " << A << endl;
//X   exit(-1);
  double h_M = tau2rho_cdm(tau);
  double F_M = 1/h_M*S*qdot/adotoa*(Vq-Vc);         

  //now we calculate d/dtau of Dg and V for cdm

  yprime[3] = -3*adotoa*Dgc*q_M-k*Vc + a/tau2rho_cdm(tau)*(A*Q_M-Qdot_M*Phi/adotoa+E_SM);
  yprime[4] =  k*Psi - adotoa * Vc + adotoa*F_M;



  double coupled = isTightCoupling(tau,photon);

  double R = cosmos->tau2R(tau);  // 4/3 * rho_photon / rho_baryon
  double cs2 = soundSpeed(tau);  // baryon soundspeed (almost 0, for tau > 1e-2) 
  double Vb_prime; // d/dtau baryon velocity
  double Vp_prime; // d/dtau  photon velocity 

  // the equations are simpler if still in tight coupling
  if (coupled) {
    Vb_prime = ( adotoa *Vb*(3*cs2 - 1) + k*R*(M[0]) + (R+1)*k*Psi - k*Phi*(R + 3*cs2) + k*cs2*Dgb) / (1+R);
    Vp_prime = Vb_prime;
    for (int l = 2; l <= lmaxg; l++) {
      Mprime[l] = 0;
      Eprime[l]  = 0;
    }
  }   else {
    // den cs2 kram kann man sich eigentlich schenken fuer die grossen taus hier.
    Vb_prime  = k*(Psi -3*cs2*Phi) + Dgb*k*cs2 - (adotoa - 3*cs2)*Vb + opac(tau)*R*(Vp - Vb);
  }

  yprime[1] = -k*Vb ; 
  yprime[2] = Vb_prime;

 
  // now we propagate the Photon Multipoles  
  Mprime[0] = -k/3.0 * Vp;
  
  // if not tightly coupled, we need to compute the Boltzmann hierarchy
  if (! coupled) {
    Vp_prime = -(Vb_prime + adotoa*Vb*(1 - 3*cs2) - k*cs2*Dgb + 3*cs2*k*Phi)/R +
      k*(M[0]-2.0/5.0*M[2]) - k*Phi + (R+1)/R * k*Psi; 
    Mprime[2] = -opac(tau)*(9.0/10.0 * M[2] + s6/10.0*E[2]) + k*(2.0/3.0 * Vp - 3.0/7.0*M[3]);
    for (int l = 3; l < lmaxg; l++) 
      Mprime[l]  = - opac(tau) * M[l] + k*( denom1[l] * M[l-1] - denom2[l]*M[l+1]);
    
    // truncation of the scheme exactly like in Ma & Bertschinger 
    Mprime[lmaxg] =  ((double) 2*lmaxg +1)/((double) 2*lmaxg-1)*k*M[lmaxg - 1] - M[lmaxg]*( (lmaxg+1)/tau + opac(tau));
 
    // propagation of E - Polarization for photons
    Eprime[2] = -k*( s5 / 7.0 *E[3] ) - opac(tau) * (E[2] + s6/10.0 * (M[2] - s6*E[2]));
    for (int l = 3; l < lmaxg; l++) 
      Eprime[l]  = k*( denomE1[l] * E[l-1] - denomE2[l]*E[l+1] ) - opac(tau)*E[l];
    // truncation of the scheme exactly like in Ma & Bertschinger 
    // OLD Eprime[lmaxg] = k*E[lmaxg - 1] - E[lmaxg]*( (lmaxg+1)/tau + opac(tau));
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

  // now, we can propagate Dgq and Vq

  double wdot = quintcosmos->tau2wdot_q(tau);
  double cq = wq - wdot/(3.0*adotoa*(1+wq)*(1-q_q));  
  double wGamma = (1-cq)*(Dgq-3.0*(1+wq)*(1-q_q)*(Phi-adotoa/k*Vq));
  double Qdot_q= - Qdot_M;
  double E_SQ  = Qdot_q/(k)*Vq+(-1/a*qdot*DS-1/a*S*Xdot-Qdot_q/k*Vq-Q_q*Psi); // this has
  // been multiplied with Q_q, otherwise we get into trouble if there is no coupling

  
  yprime[qidx] = -3.0*adotoa*Dgq*(cq-wq+q_q*(1+wq))-k*(1+wq)*Vq-
    3.0*adotoa*wGamma+a/quintcosmos->tau2rho_q(tau)*(Q_q*A-1/adotoa*Qdot_q*Phi+E_SQ); 
  yprime[qidx+1] = k*(Psi-3.0*(1-q_q)*Phi)+k/(1+wq)*Dgq+adotoa*(3.*(1-q_q)-1)*Vq; 


//X #warning added as test
//X   Nprime[2] = 0;
  static bool say = true;
  if (say) cout <<  yprime[1] << "##c##" << Vc << " is Vgc " << endl;
  say = false;
#ifdef OUTPUTFILE
  (*ofs) << setprecision(18) << tau << " " << k << " " << N[2] << " " << Psi << " " << Pin << " " << Pip << endl;
//X   for (int i = 0; i <= nvar; ++i)
//X      (*ofs)  << y[i] << " " << yprime[i] << " ";
//X   (*ofs) << endl;
#endif //OUTPUTFILE
//X  (*ofs) << tau << " " << k << " " << Pin << " " << Pip << endl;

  //  cout <<  1/a << " " << -1/a*qdot*DS <<" " << qdot << " "<< DS << " " << yprime[qidx] <<" " << Q_q*Psi <<endl; 
  //  cout << yprime[qidx] <<" " << Q_q << " " << E_SQ << "  " << Qdot_q  << endl;
  //  cout << q << " " << DS << " " << betaphiphi*qdot*Vq/k*q*tau2rho_cdm(tau)/Mpl <<  " " << Mpl << endl;
//X #ifdef OUTPUTFILE
//X //X //X   (*ofs) << tau << " " << cq << " " << wGamma << " " << Qdot_q << " " << E_SQ << endl;
//X   (*ofs) << tau
//X   << " " << k
//X   << " " << Vq
//X   << " " << Dgq
//X   << " " << Vc
//X   << " " << Dgc
//X   << " " << Dgb
//X   << " " << Psi
//X   << " " << Phi
//X   << " "  << qdot
//X  << " "  <<adot
//X  << " "  <<Mpl
//X  << " "  <<wq
//X  << " "  <<beta
//X  << " "  <<betaphi
//X  << " "  <<betaphiphi
//X  << " "  <<S
//X  << " "  <<Q_M
//X  << " "  <<Q_q
//X  << " "  <<Sdot
//X  << " "  <<Qdot_M
//X  << " "  <<q_M
//X  << " "  <<h_q
//X  << " "  <<q_q
//X  << " "  <<DS
//X  << " "  <<Xdot
//X  << " "  <<E_M
//X  << " "  <<E_SM
//X  << " "  <<A
//X  << " "  <<h_M
//X  << " "  <<F_M
//X  << endl;
//X #endif
}

/*!
  For the sources, to be precise the ISW effect, one needs d/d tau
  of Phi and Psi. In the CoupledInvariant class, these are needed anyhow
  within the fderivs() function. However, here they are only calculated
  for the sources, which is why we need to do it.

  Please note that fillPhiDot() assumes that it is called after fderivsWrapper() 
  has been called, i.e. Phi and Psi and all other stuff is up to date.
*/
void CoupledInvariant::fillPhiDot(double tau) {
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
  RhoPV += (1+quintcosmos->tau2w_q(tau))*quintcosmos->tau2rho_q(tau)*Vq;
  RhoPV += 4.0/3.0*(tau2rho_g(tau)*Vp + tau2rho_nu(tau)*Vn);
  PhiDot  = adotoa*Psi - Gpi4a2*(RhoPV)/k;

  // now d/dtau Psi from d/dtau of eqn 2.52, first d/dtau (p*Pi) p = rho/3
  double Psi2nd = 12.0/5.0*a2*(tau2rhodot_g(tau)*M[2] + tau2rhodot_nu(tau)*N[2] 
			       + tau2rho_g(tau)*Mprime[2] + tau2rho_nu(tau)*Nprime[2])/3.0;
  // add to this d/dtau a^2 term
  Psi2nd += 2.0*a*tau2adot(tau)*(tau2p_g(tau)*Pip + tau2p_nu(tau)*Pin);
  PsiDot = -PhiDot - Gpi8()/(k*k) * Psi2nd;
}

void CoupledInvariant::calcPerturbationNr(const ControlPanel &control) {
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

  if (!UseWrongNumberOfPertEquationLikeOriginal)
    nvar += 2; //  2 more for dark energy

  y = new double[nvar+1];
  yprime = new double[nvar+1];
  yt = new double[nvart+1];
  ytprime = new double[nvart+1];

}

void CoupledInvariant::getReady(const ControlPanel& control) {
  quint = quintcosmos->quintessence(); 
  Invariant::getReady(control);  // now call parent objects get ready
  qidx = nvar-1;
}

/*! Set the initial conditions; only adiabatic ones will be accepted, others are not yet implemented in
CoupledInvariant. */
void CoupledInvariant::initialScalarPerturbations(const ControlPanel &control, const double tau)
{

  //  cout << "initialScalarPerturbations from coupledinvariant\n";

#ifdef OUTPUTFILE
  ofs = new ofstream("coupledinvcolddark.dat");
#endif


  mona_a = tau2a(tau); // Mona

  
  // unfortunatly, we need to define our M,N,E ourselves, cause
  // they are non -constant (we have to initialize them)
  // se setPointers() for more  
  double *M = &y[5];                    // Temperature
  double *E = &M[lmaxg-1];  // E[2] is first after M[lmaxg]
  double *N = &E[lmaxg+1];  //  massles neutrinos

  double Dgp, Dgc, Dgb, Dgn, Dgq,Vb,Vc,Vp,Vn,Psi,Pi_nu,Vq;
  k2 = k*k;
 
  double x = k*tau;

  double Omega_nu = cosmos->rho2omega(cosmos->tau2rho_nu(tau),tau);

  double wq=quintcosmos->tau2w_q(tau);
  double P  = 1.0/(4.0*Omega_nu + 15.0);
 
  double normalize = 1.0; // see comment below on rescaling

  switch (control.initialCond) {
  case ControlPanel::adiabatic: 
    Dgp = 1;
    Dgn = Dgp;
    Dgb = 3.0/4.0*Dgp;
    Dgc = Dgb;
    Dgq = 3.0/4.0*(1.0+wq)*Dgp;
    Vp = -5.0/4.0*P*x*Dgp;
    Vb = Vp;
    Vn = Vp;
    Vc = -5.0/4.0*P*x*Dgp;
    Vq = -5.0/4.0*P*x*Dgp;
    Pi_nu = -Dgp*P*x*x;
    Phi = (5.0+2.0*Omega_nu)/(30.0+8.0*Omega_nu)*Dgp; 
    Psi= -Phi-Omega_nu*Gpi8()*Pi_nu/(x*x);
    normalize = -1.0/Psi;  // this is for tensor / scalar ratio , only for adiabatic see below
    break;
  default:
    cout << "Unknown initial conditions - CoupledInvariant takes only adiabatic initial conditions" << endl;
    throw exception();
    break;
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
  Dgq *= normalize;
  Vp *= normalize;
  Vb *= normalize;
  Vn *= normalize;
  Vc *= normalize;
  Vq *= normalize;
  Pi_nu *= normalize;
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
  N[0] = Dgn / 4.0;
  N[1] = Vn;
  N[2] = 5.0/12.0*Pi_nu;
  
  // set all higher moments to 0, obs: E[2] is the first one needed
  for (int l = 2; l <= lmaxg; l++) {M[l] =0; E[l] = 0;}
  for (int l = 3; l <= lmaxnr; l++) N[l] =0; 

  // fill quintessence
  y[qidx]=Dgq;
  y[qidx+1]=Vq;
} 

void CoupledInvariant::scalarSources(double tau, double *d, double *dp, double *dk) { // added by Mona!!!!!!!! 

  //  cout << "das ist scalarSources! von coupledinvariant!\n";

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
   
#ifdef OUTPUTFILE
   (*ofs) << tau << "  " << isw_mona << "  " << acoustic_mona << "  " << doppler_mona << "  " << fourth_mona << "  " << acoustic_mona + doppler_mona << "  " << isw_mona + acoustic_mona + doppler_mona << "  " << isw_mona + acoustic_mona + doppler_mona + fourth_mona << endl;
#endif
   
   

  
   // Mona end

   double x = k * (tau_0() - tau);
   if (x > 0.) *dp = visibility(tau) * 3. / 2. * C / (x * x); else *dp = 0.;
  
} 
