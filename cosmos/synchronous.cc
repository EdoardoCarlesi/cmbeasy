#include "synchronous.h"
#include "cosmos.h"

//#define OUTPUTFILE
Synchronous::Synchronous(Cosmos *cosmos) : Perturbation(cosmos), Psi(0) {
  // to speed up, precompute  1 / (2*l + 1) 
  for (int j = 1; j <= LMX0; ++j)  denl[j] = 1. / (double) ((j << 1) + 1);
}

/*!
  This is (more or less) the cmbfast fderivs() function. The difference is (at least
  should be), that Cosmos provides background quantities. 
  Please note: Due to some mystery, it makes a tiny difference, approx 0.05 %
  whether one chooses the variables here as static or not. 
  I have not the slightest idea, why this should be so, so if someone knows
  the answer, please let me know. Most propable it is a compiler implementation
  issue.
  The bad thing about static variables is that it may prevent parallalization of the
  code (I am not entirely sure on this issue, but I believe so). Hence, when 
  thinking about a threaded cmbeasy, one should most propably change to 
  none - static variables. Static, btw is about 1/2 % faster
*/
void Synchronous::fderivs(const double tau, const double *y, double *yprime) {

  called++;

   double drag, dpnu, slip, adotdota, a;
   int l;
    double dgrho, gpres, a2, deltabdot, deltagdot,
    thetabdot, deltardot, thetagdot, thetardot,
    deltac, deltab, deltag, thetac, thetab, shearg, thetag,
    deltar, etadot, shearr, thetar;

   double drhonu, polter, cs2;

   double pb43, eta;
   double fnu, dgtheta, pnu;
   bool coupled;

  //  Evaluate the time derivatives of the perturbations. 
  // ep is used to stop the tight coupling approximation. 

  int lmaxg2 = lmaxg << 1;        // 2*lmaxg is frequently needed

  a= tau2a(tau);   // we take it from cosmos, as this is also mother, children will be fast
  a2 = a * a;
  adotoa = tau2adot(tau)/a; 

  eta = y[1];

  // CDM
  deltac = y[2];   // delta of cold dark matter
  thetac = y[3];

  if (thetac != 0) throw Bad_Error("none");

  // Baryons
  deltab = y[4];  // delta of baryons
  thetab = y[5];


  // Photons. 
  deltag = y[6];
  thetag = y[7]; 
  shearg = 0.5*y[8];

  //  Polarization term.
  polter = y[8] + y[ lmaxg + 7] + y[ lmaxg + 9];

  //  Massless neutrinos.
  deltar = y[lmaxg2 + 8];
  thetar = y[lmaxg2 + 9];
  shearr = 0.5*y[lmaxg2 + 10];

  if (cosmos->nuR()==0) {
    deltar=thetar=shearr=0;
  }

  // Get soundspeed squared from interpolation:  cs2

  cs2 = soundSpeed(tau);

  coupled = isTightCoupling(tau,photon);

  //  Photon mass density over baryon mass density.
  pb43 = cosmos->tau2R(tau);

  double neutrinoMass_kbT = 0;
  if (cosmos->haveMassiveNeutrinos()) {
    neutrinoMass_kbT = cosmos->mass_nu_kbT();
    MassiveNeutrinos::nu1(neutrinoMass_kbT*a, &rhonu, &pnu);
    MassiveNeutrinos::nu2(neutrinoMass_kbT*a, 0. /*no coupling to DE*/,
                          rhonu, pnu, 0. /* field perturbation, only needed when coupled to dark energy*/,
                          &drhonu, &fnu, &dpnu, &shearnu, &y[iq0], &y[iq1], &y[iq2]);
  } else {
    rhonu = 1.;
    pnu = 1.0/3.0;
    drhonu = 0.;
    fnu = 0.;
    dpnu = 0.;
    shearnu = 0.;
  }

  gpres = Gpi8()*a2*tau2p(tau);

  //  Evaluate metric and massive neutrino perturbations.
  //deltan=drhonu/rhonu
  //thetan=ak*fnu/(rhonu+pnu)

  //  8*pi*G*delta_rho*a**2 and 8*pi*G*delta_P*a**2.

  // [rho_xx]  Mpc^-4
  // [Gpi8   ]  Mpc^ 2
  // -- >                        [dgrho ]   Mpc^-2

  dgrho = Gpi8()*a2* (tau2rho_cdm(tau) * deltac + tau2rho_b(tau)*deltab
                      + tau2rho_g(tau)*deltag + tau2rho_nu(tau)*deltar  + drhonu/rhonu*tau2rho_nuNR(tau));

  //  Add a seed if desired.
  if (seedFlag() )  dgrho += Gpi8()*rho_0() / a; // COM dgrho += grhom / a;

  //  Force energy conservation.
  hdot = (k2 * 2. * eta + dgrho) / adotoa;

  // dgtheta is 8pi G a^2 times  (rho + p) theta_species

  dgtheta = Gpi8()*a2*( tau2rho_cdm(tau) * thetac + tau2rho_b(tau) * thetab
                        + 4.0/3.0 * (tau2rho_g(tau) * thetag + tau2rho_nu(tau)*thetar)
                        + k*fnu/rhonu*tau2rho_nuNR(tau));
  etadot = dgtheta * .5 / k2;
  yprime[1] = etadot;

  //  8*pi*G*(rho+P)*sigma*a**2.
  dgshear = Gpi8()*a2*( 4.0/3.0*  (tau2rho_g(tau) * shearg + tau2rho_nu(tau)*shearr)
                        + shearnu/rhonu*tau2rho_nuNR(tau));

  //  CDM equations of motion.

  yprime[2] =  -thetac - hdot * .5;  //  deltacdot;
  yprime[3] = -adotoa * thetac;      //  thetacdot;

  if (yprime[3] !=0) throw Bad_Error("nein!!!");

  //  Baryon equations of motion.
  deltabdot = -thetab - hdot * .5;
  yprime[4] = deltabdot;

  //  Need photon perturbation for first-order correction to tightly-coupled
  //  baryon-photon approximation.
  deltagdot = (-thetag - hdot * .5) * 4.0/3.0;
  drag = opac(tau) * (thetag - thetab);
  if (! isTightCoupling(tau,baryon)) {   //  Treat baryons as uncoupled from photons
    thetabdot = -adotoa * thetab + k2 * cs2 * deltab + pb43 * drag;
  } else {
    //  Treat baryons and photons as tightly coupled. 
    //  Zeroth-order approximation to baryon velocity. 
    thetabdot = (-adotoa * thetab + k2 * cs2 * deltab 
		 + k2 * pb43 * (deltag * .25 - shearg)) / (pb43 + 1.);
    //  (\ddot a)/a. 
    adotdota = (adotoa * adotoa - gpres) * .5;

    //  First-order approximation to baryon-photon slip, thetabdot-thetagdot.
    // this is the r.h.s of eqn (74), Ma & Bertschinger
    slip = pb43 * 2. / (pb43 + 1.) * adotoa * (thetab - thetag) + 
      1. / opac(tau) * (-adotdota * thetab - adotoa * k2 * .5 * deltag 
			+ k2 * (cs2 * deltabdot - deltagdot * .25)) / (pb43 + 1.);
    //  First-order approximation to baryon velocity. 
    thetabdot += pb43 / (pb43 + 1.) * slip;
  }
  yprime[5] = thetabdot;

  //  Photon total intensity and polarization equations of motion.
  yprime[6] = deltagdot;
  thetagdot = (-thetabdot - adotoa * thetab + k2 * cs2 *
               deltab) / pb43 + k2 * (deltag * .25 - shearg);
  yprime[7] = thetagdot;

  if (!coupled) {           //  Treat photons as uncoupled from baryons
    // eqn 63, Ma & Bertschinger, F^g_2 (shearg)
    yprime[8] = 8.0/15.0 * thetag - 0.6* k * y[9]
                - opac(tau) * y[8] + 4.0/15.0 * hdot + 8.0/5.0*etadot + opac(tau) * .1 * polter;
    //  Polarization equations for l = 0, 1, 2.
    yprime[lmaxg + 7] = -k * y[ lmaxg + 8] - opac(tau) * y[ lmaxg + 7] + opac(tau) * .5 * polter;
    yprime[lmaxg + 8] = k / 3. * (y[lmaxg + 7] - y[ lmaxg + 9] * 2.) - opac(tau) * y[ lmaxg + 8];
    yprime[lmaxg + 9] = k * (y[ lmaxg + 8] * .4 - y[ lmaxg + 10] * .6) - opac(tau) * y[ lmaxg + 9] + opac(tau) * .1 * polter;
    for (l = 3; l < lmaxg; ++l) {
      yprime[l + 6] = k * denl[l] * (l * y[ l + 5] - (l + 1) * y[ l + 7]) - opac(tau) * y[ l + 6];  // d/dtau F^g_l = ...
      yprime[lmaxg + 7 + l] = k * denl[l] * (l * y[ lmaxg + 6 + l] - (l + 1) * y[ lmaxg + 8 + l]) - opac(tau) * y[ lmaxg + 7 + l];
    }
  } else {
    //  Treat baryons and photons as tightly coupled (with no polarization). 
    yprime[8] = 0.;
    yprime[lmaxg + 7] = 0.;
    yprime[lmaxg + 8] = 0.;
    yprime[lmaxg + 9] = 0.;
    for (l = 3; l < lmaxg; ++l) {
      yprime[l + 6] = 0.;
      yprime[lmaxg + 7 + l] = 0.;
    }
  }

  //  Truncate moment expansion.
  //    yprime(8+lmax)=ak*lmax*y(7+lmax)/(2*lmax+1)-opac*y(8+lmax)
  //    yprime(9+2*lmax)=ak*lmax*y(8+2*lmax)/(2*lmax+1)-opac*y(9+2*lmax)
  yprime[lmaxg + 6] = k * y[lmaxg + 5] - (lmaxg + 1) / tau * y[lmaxg + 6] - opac(tau) * y[lmaxg + 6];
  yprime[lmaxg2 + 7] = k * y[lmaxg2 + 6] - (lmaxg + 1) / tau * y[lmaxg2 + 7] - opac(tau) * y[lmaxg2 + 7];

  if (cosmos->nuR()>0) {
    //  Massless neutrino equations of motion.
    deltardot = 4.0/3.0 * (-thetar - hdot * .5);
    yprime[lmaxg2 + 8] = deltardot;
    thetardot = k2 * (deltar * .25 - shearr);
    yprime[lmaxg2 + 9] = thetardot;

    // F^nu_2   :   y[lmaxg2 + 10]
    // F^nu_3   :   y[lmaxg2 + 11]
    yprime[lmaxg2 + 10] = 8.0/15.0*thetar - 0.6*k * y[lmaxg2 + 11] + 4.0/15.0*hdot + 8.0/5.0*etadot;

    // F^nu_l    :   y[lmaxg2 + 8 +  l]
    for (l = 3; l < lmaxnr; ++l) {
      yprime[lmaxg2 + 8 + l] = k * denl[l] * (l * y[lmaxg2 + 7 + l] - (l + 1) * y[lmaxg2 + 9 + l]);
    }
    //  Truncate moment expansion.
    yprime[lmaxg2 + 8 + lmaxnr] = k * y[lmaxg2 + 7 +
      lmaxnr] - (lmaxnr + 1) / tau * y[lmaxg2 + 8 + lmaxnr];
  } else {
    deltardot = yprime[lmaxg2 + 8] = thetardot
      = yprime[lmaxg2 + 9] = yprime[lmaxg2 + 10]=0;
    for (l = 3; l < lmaxnr; ++l) {
      yprime[lmaxg2 + 8 + l] = 0;
    }
    yprime[lmaxg2 + 8 + lmaxnr] = 0;
  }
#ifdef OUTPUTFILE
  double xx = k*tau;
  (*ofs) << xx << " " <<   deltag  << endl;
#endif
  //double Phi = eta - 0.5 / k2 * adotoa * ( hdot + 6 * etadot );
  Psi = -1.5 * dgshear / k2 - adotoa + eta - adotoa / ( 2 * k2 ) * ( hdot + 6 * etadot );


  if (cosmos->haveMassiveNeutrinos()) {
    MassiveNeutrinos::propagateMassiveNeutrinoMomentsSynchronous(y, yprime, a, tau,
                                   neutrinoMass_kbT, k, hdot, etadot);
  }
}

void Synchronous::initialScalarPerturbations(const ControlPanel &control, const double tau)
{
#ifdef OUTPUTFILE
  ofs = new ofstream("colddark.dat");
#endif
   double  a, C, h;
   int l;
   double a2, f1, delta0,
          deltac, deltab, deltag, thetab, thetac, deltan,
          deltar, thetag, thetan, shearr, thetar, ahdot;
   double eta;
   double psi, kt2;

  //  Initial conditions.

  k2 = k*k;
  a = tau2a(tau);
  a2 = a * a;

  double sum = tau2rho(tau) + tau2p(tau); // COM s = grho + gpres
  double  fracnu = (tau2rho_nu(tau) + tau2rho_nuNR(tau) + tau2p_nu(tau) + tau2p_nuNR(tau) )  / sum;

  double yrad = tau2rho_m(tau) /  ( tau2rho_g(tau) + tau2rho_nu(tau) + tau2rho_nuNR(tau));

  double Omega_nu = cosmos->rho2omega(cosmos->tau2rho_nu(tau),tau) + cosmos->rho2omega(cosmos->tau2rho_nuNR(tau), tau);

  switch (control.initialCond) {
  case ControlPanel::adiabatic:

    //-----------------------------------------------------------
    //  First case.
    //  Isentropic ("adiabatic") initial conditions.
    psi = -1 / (4*Omega_nu + 15.0)*10;
    C = (fracnu * 4. + 15.) / 20. * psi;
    kt2 = pow(k*tau,2);
    h = C * kt2 * (1. - yrad * .2);
    eta = C * 2. - (fracnu * 4. + 5.) / 6. / (fracnu * 4. + 15.) * C * kt2 * (1. - yrad / 3.);
    f1 = (fracnu * 4. + 23.) / (fracnu * 4. + 15.);
    deltac = h * -.5;
    deltag = h * -2.0/3.0 * (1. - kt2 / 36.);
    deltab = deltac;
    deltar = h * -2.0/3.0 * (1. - kt2 / 36. * f1);
    thetac = 0.;
    thetag = -C / 18. * kt2 * kt2 / tau;
    thetab = thetag;
    thetar = f1 * thetag;
    shearr = k2 * 4.0/15.0 / (Gpi8()*a2*sum) * psi * (yrad * 7.0/36.0 + 1.);
    ahdot = C * 2. * k2 * tau * a * (1. - yrad * .3);
    break;
  case ControlPanel::isoCDM:
    //  Isocurvature CDM initial conditions: perturb only CDM as a --> 0.
    delta0 = 1.;
    h = delta0 * yrad * (1. / (omega_b() / omega_c() + 1.) - yrad * .5);
    deltac = delta0 - h * .5;

    //  Compensate perturbation with everything else.
    deltag = h * -.66666666666666663;
    deltab = deltag * .75;
    deltar = deltag;
    thetac = 0.;
    thetag = -h / 12. * k2 * tau;
    thetab = thetag;
    thetar = thetag;
    shearr = 0.;
    //ahdot = adotrad * h * (1. - yrad * .5);
    eta = -h / 6.;
    break;
  case ControlPanel::isoBaryon:
    //  Isocurvature baryon initial conditions: perturb only baryons as a ->0.
    delta0 = 1.;
    h = delta0 * yrad * (1. / (omega_c() / omega_b() + 1.) - yrad * .5);
    deltab = delta0 - h * .5;

    //  Compensate perturbation with everything else.
    deltac = h * -.5;
    deltag = h * -.66666666666666663;
    deltar = deltag;
    thetac = 0.;
    thetag = -h / 12. * k2 * tau;
    thetab = thetag;
    thetar = thetag;
    shearr = 0.;
    //ahdot = adotrad * h * (1. - yrad * .5);
    eta = -h / 6.;
    break;
  case ControlPanel::isoSeed:
    //  Fourth case. 
    //  Isocurvature seed initial conditions:everything is unperturned as a->0
    delta0 = 1.;
    h = delta0 * yrad * (1. / (omega_c() / omega_b() + 1.) - yrad * .5);

    //  Compensate perturbation with everything else.
    deltab = h * -.5;
    deltac = h * -.5;
    deltag = h * -.66666666666666663;
    deltar = deltag;
    thetac = 0.;
    thetag = -h / 12. * k2 * tau;
    thetab = thetag;
    thetar = thetag;
    shearr = 0.;
    //ahdot = adotrad * h * (1. - yrad * .5);
    eta = -h / 6.;
    break;
  default:
    throw Bad_Error("fintial, strange intital condition : unknown");
  }


  deltan = deltar;
  thetan = thetar;

  y[1] = eta;

  // CDM.
  y[2] = deltac;
  y[3] = thetac;

  // Baryons.
  y[4] = deltab;
  y[5] = thetab;

  // Photons (total intensity and polarization).
  y[6] = deltag;
  y[7] = thetag;
  // shearg=0.0d0
  y[lmaxg + 7] = 0.;
  y[lmaxg + 8] = 0.;

  for (l = 2; l <=lmaxg ; ++l) {
    y[l + 6] = 0.;
    y[lmaxg + 7 + l] = 0.;
  }

  //  Massless neutrinos.
  if (cosmos->nuR()>0) {
    y[(lmaxg << 1) + 8] = deltar;
    y[(lmaxg << 1) + 9] = thetar;
    y[(lmaxg << 1) + 10] = shearr * 2.;

    for (l = 3; l <=lmaxnr ; ++l)
      y[(lmaxg << 1) + 8 + l] = 0.;
  }

  //  Massive neutrinos.
  if (cosmos->haveMassiveNeutrinos() && (control.initialCond != ControlPanel::adiabatic)) {
    throw Bad_Error("For massive neutrinos, only adiabatic initial conditions are implemented.");
  } else if (cosmos->haveMassiveNeutrinos()) {
    MassiveNeutrinos::setFirstIndex(iq0);
    const double mass_nu = cosmos->mass_nu_kbT();
    const double Vn=thetan/k;
    const double Pi_nu=6.*shearr;
    MassiveNeutrinos::setIninitalLongitudinalScalarPerturbations(y, cosmos->tau2a(tau), mass_nu, deltan, Vn, Pi_nu);
  }

}

#define LM2 7
#define LM3 4
#define NQ1 15 

void Synchronous::getReady(const ControlPanel& control) {
  called=0;
  if (control.scalar || control.power_cdm) {
    lmaxg = LMAX0;
    if (control.highPrecisionTransfer) {
      lmaxnr = LMAXNR0;
    } else {
      lmaxnr = LM2;
    }
  }

  if (control.tensor) lmaxt = LMAXT0; else lmaxt = 0;

  Perturbation::getReady(control);  // now call parent objects get ready
}

void Synchronous::scalarSources(double tau, double *d, double *dp, double *dk) {
  double chir, alphadot, coupldot, rhonudot, a, alpha, x, glens;
  double alphaddot, a2, thetabdot, adotoadot, sheargdot, thetagdot;
  double s1, s2, shearrdot, polterdot, deltag, thetab, etadot;
  double dgsheardot, polter, shearnudot, polterddot, chi, eta, phi;

  fderivs(tau,y,yprime);   // first, we do the fderivs 

  x = k * (tau_0() - tau);
  a = tau2a(tau);
  a2 = a * a;

  // ahdot=y(2) 
  eta = y[1];
  etadot = yprime[1];
  alpha = (hdot + etadot * 6) / (k2 * 2.);
  alphadot = -3.0* dgshear / (k2 * 2.) + eta - adotoa * 2. * alpha;

  // Baryons. 
  // deltab=y(6) 
  thetab = y[5];
  thetabdot = yprime[5];

  //  Photons. 
  deltag = y[6];

  // thetag=y(9) 
  // shearg=y(10)/2.0d0 

  thetagdot = yprime[7];
  sheargdot = 0.5*yprime[8];

  //  Polarization term. 
  polter = y[8] + y[lmaxg + 7] + y[lmaxg + 9];
  
  // coupl=8.0d0*(thetag+ak2*alpha)/15.0d0-ak*0.6d0 
  // 2 *(y(11)+y(10+lmaxg)+y(12+lmaxg)) 
  coupldot = (thetagdot + k2 * alphadot) * 8. / 15. - k * .6 * (yprime[9] + yprime[lmaxg + 8] + yprime[lmaxg + 10]);
  polterdot = yprime[8] + yprime[lmaxg + 7] + yprime[lmaxg + 9];
  polterddot = coupldot - (dopac(tau) * polter + opac(tau)  * polterdot) * .3;

  // Massless neutrinos. 
  //   deltar=y(10+2*lmaxg) 
  //   thetar=y(11+2*lmaxg) 
  //   shearr=y(12+2*lmaxg)/2.0d0 
  shearrdot = 0;
  if (cosmos->nuR() != 0) {
    shearrdot = 0.5*yprime[(lmaxg << 1) + 10];
  }

  rhonudot = 0.;
  if (cosmos->nuNR() == 0.0) {
    shearnudot = 0.;
  } else {
    shearnudot = MassiveNeutrinos::nuder(a*cosmos->mass_nu_kbT(), adotoa,
                                                0 /*no coupling*/, tau2rho_nuNR(tau),
                                                0 /*field pert, only needed for coupling*/,
                                                y+iq2, yprime+iq2);
  }

  adotoadot = tau2grhodot(tau) / (adotoa * 6.0);  // d/dtau adotoa

  //  Derivative of the shear
  dgsheardot = Gpi8()*a2* (4.0/3.0*(tau2rho_g(tau) * sheargdot + tau2rho_nu(tau)*shearrdot)
                            + shearnudot/rhonu*tau2rho_nuNR(tau)) - 2.0 * adotoa * dgshear;
  //  Calculation of the sources 
  alphaddot = -3.0 * dgsheardot / (k2 * 2.) + etadot - 2.0*( adotoadot * alpha +  adotoa * alphadot );

  s1 = etadot + alphaddot; // ISW
  s2 = alphadot * 2;

 *d = expmmu(tau) *s1 
    + visibility(tau)* (deltag * .25 + s2 + polter / 16 + thetabdot / k2 + polterddot * .1875 / k2)
    + dvisibility(tau) * (alpha + thetab / k2 + polterdot * .375 / k2) 
    + ddvisibility(tau) * 3. / 16. * polter / k2;


  if (x > 0.) *dp = visibility(tau) * 3. / 16. * polter / (x * x); else *dp = 0.;

  // CMBFAST stuff:
  // lensing visibility function; approximate epoch of recombination 
  chi = tau_0() - tau;
  chir = tau_0() - tau_ls(); // maxreion_.taurmax;
  if (chi < chir) {
    // lensing convergence, expmmu is an approximation 
    // one should integrate visibility function over r(chi)/r(tau) 
    // but the error is harmless 
    phi = eta -  adotoa * alpha;
    glens = (chir - chi) * chi / chir;
    *dk = glens * k2 * phi * expmmu(tau);
    // alternative form ala Stebbins; add l(l+1) in the final Cl
    // it is numerically inaccurate at high l, equivalent to above at low l
    //            glens=(chir-chi)/chi/chir
    //            dk=glens*phi*expmmu(j)
    // using Poisson's eq. below gives identical results
    //            deltac=y(4)
    //            dk=glens*grhom*(omegac*deltac+omegab*deltab)/2.0/a*expmmu(j)
  } else {
    *dk = 0.;
  }
}

/*
void Synchronous::fderivsTensor(const double tau,   double *y, double *yprime) {
    // Local variables 
     double  psie, htpr, a;
     int l;
     double a2, htdpr, rhonu, ep, ht, adotoa;
  
     double cs2;
     double deltap0, deltat0, pnu;
     int ind1, ind2;
    
    //  Evaluate the time derivatives of the perturbations. 
    // ep is used to stop the tight coupling approximation. 
    if (k > epsw * .06) ep = .01; else ep = .0117;
        

    a = tau2a(tau);
    a2 = a * a;
    cs2 =soundSpeed(tau);

    // Tight Coupling parameters 
   
    double tcp = 0.;
    double tcp1 = k / opac(tau);
    double tcp2 = 1. / (opac(tau) * tau);
    if (tcp1 > ep || tcp2 > ep) tcp = 1.;
    
    bool uncoupled = (tcp == 1.0);

    //  Compute expansion rate. 
    if (cosmos->amnu == 0.) {
	rhonu = 1.;
	pnu = .33333333333333331;
    } else cosmos->nu1(a, &rhonu, &pnu);
    
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
 
    if (uncoupled) {
      // no tight coupling approx 
	yprime[ind1] = -k * y[ind1 + 1] - opac(tau) * y[ind1] + opac(tau) * psie - htpr;
 	yprime[ind2] = -k * y[ind2 + 1] - opac(tau) * y[ind2] - opac(tau) * psie;
	for (l = 1; l < lmaxt ; ++l) {       // up to lmaxt-1
	    yprime[ind1 + l] = k * denl[l] * (l * y[ind1 - 1 + l] - (l + 1) * y[ind1 + 1 + l]) - opac(tau) * y[ind1 + l];
	    yprime[ind2 + l] = k * denl[l] * (l * y[ind2 - 1 + l] - (l + 1) * y[ind2 + 1 + l]) - opac(tau) * y[ind2 + l];
	}
	// Truncate moment expansion  , now LMAXT is reached
	yprime[ind1 + lmaxt] = k * y[ind1 - 1 + lmaxt] - (lmaxt + 1) / tau * y[ind1 + lmaxt] - opac(tau) * y[ind1 + lmaxt];
	yprime[ind2 + lmaxt] = k * y[ind2 - 1 + lmaxt] - (lmaxt + 1) / tau * y[ind2 + lmaxt] - opac(tau) * y[ind2 + lmaxt];
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

void Synchronous::tensorSources(double tau, double *dt, double *dte, double *dtb) {
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


  Set the initial conditions for the tensor pertubrations at time tau.
  Slightly modified version, for we do not need to integrate a(tau), cause
  we have it already

void Synchronous::initialTensorPerturbations() {
  yt[1] = 1;   
  yt[2] = 0.;
  int ind1 = 3;
  int ind2 = ind1 + lmaxt +1;
  for (int l = 0; l <= lmaxt; ++l) {
    yt[ind1 + l] = 0.;
    yt[ind2 + l] = 0.;
  }
} 
*/


double Synchronous::delta_nr() {
  return 0;
}


void Synchronous::calcPerturbationNr(const ControlPanel &control) {
  nqmax=MassiveNeutrinos::qGridSize();
  //   Calculate number of equations
  if (control.scalar || control.power_cdm) {
    nvar = ( (lmaxg + 1) << 1) + 5;
    iq0 = (lmaxg << 1) + 8;
    if (cosmos->nuR() > 0) {
      nvar += lmaxnr +1; // massless neutrinos
      iq0 += lmaxnr+1;
    }
    iq1 = iq0 + nqmax;
    iq2 = iq1 + nqmax;
    if (cosmos->haveMassiveNeutrinos()) {
      nvar += MassiveNeutrinos::pertArraySize();
    }
  } else {
    nvar = 0;
  }

  if (control.tensor) {
    nvart = 2*lmaxt + 4 + lmaxnr+1;
  } else {
    nvart = 0;
  }

  y = new double[nvar+1];
  yprime = new double[nvar+1];
  yt = new double[nvart+1];
  ytprime = new double[nvart+1];
}




