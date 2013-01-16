#include "quintsynchronous.h"
#include "quintcosmos.h"

QuintSynchronous::QuintSynchronous(QuintCosmos *cosmos) : Synchronous(cosmos) , quintcosmos(cosmos) {}

#define OUTPUTFILE
void QuintSynchronous::fderivs(const double tau, const double *y, double *yprime) {
  //cout << "QuintSynchronous::fderivs" << endl;
  // cout << "fde: " << tau << endl;
    static double drag, dpnu, slip, adotdota, a;
    static int l;
    static  double dgrho, gpres, a2, deltabdot, deltagdot,
                   thetabdot, deltardot, thetagdot, thetardot,
                   deltac, deltab, deltag, thetac, thetab, shearg, thetag,
                   deltar, etadot, shearr, thetar;

     static double drhonu, polter, cs2;

     static double pb43, eta;
     static double fnu, dgtheta, pnu;
     static bool coupled;


     //  Evaluate the time derivatives of the perturbations. 
 
    static int lmaxg2 = lmaxg << 1;        // 2*lmaxg is frequently needed
    
    a= tau2a(tau);   // we take it from quintcosmos :-), as this is also mother, children will be fast
    a2 = a * a;
    adotoa = tau2adot(tau)/a; 
    
    // QUINTESSENCE
    // Update quintessence - field 
    quint->touch(tau);

    double DelQ = y[qidx];  // quintessence fluctuation 
    double DelQDot = y[qidx+1]; // quintessence velocity fluct.

    double sign = 1.0;
    if (quintcosmos->tau2w_q(tau) < -1.0) sign = -1.0; 
 
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
    shearg = y[8] / 2.;

    //  Polarization term. 
    polter = y[8] + y[ lmaxg + 7] + y[ lmaxg + 9];

    //  Massless neutrinos. 
    deltar = y[lmaxg2 + 8];
    thetar = y[lmaxg2 + 9];
    shearr = y[lmaxg2 + 10] / 2.;
   
    // Get soundspeed squared from interpolation:  cs2
   
    // opac = opac(); // (*splineDotmu)(tau);
    cs2 = soundSpeed(tau);   
    coupled = isTightCoupling(tau,photon);

    //  Photon mass density over baryon mass density. 
    pb43 =cosmos->tau2R(tau); // COM slight effect: photbar = grhog / (grhom * omega_b() * a); pb43 = 4.0/3.0 *photbar;

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
    // 	deltan=drhonu/rhonu 
    // 	thetan=ak*fnu/(rhonu+pnu) 
    //  8*pi*G*delta_rho*a**2 and 8*pi*G*delta_P*a**2. 

    // [rho_xx]  Mpc^-4
    // [Gpi8   ]  Mpc^ 2
    // -- >                        [dgrho ]   Mpc^-2



    double delRhoQ = sign*DelQDot*quint->qDot(a)/a2 + DelQ*quint->Vprime(quint->q(a), a, adotoa);
   
    dgrho = Gpi8()*a2* (tau2rho_cdm(tau) * deltac + tau2rho_b(tau)*deltab 
			+ tau2rho_g(tau)*deltag + tau2rho_nu(tau)*deltar  + rho_nu0()*nuNR()*drhonu/(a2*a2)
			+ delRhoQ);
   
    //cout <<ak << "  " << tau << "   " << dgrho << endl;

    //  Add a seed if desired. 
    if (seedFlag() )  dgrho += Gpi8()*rho_0() / a; // COM dgrho += grhom / a; 
   
 
    //  Force energy conservation. 
    hdot = (k2 * 2. * eta + dgrho) / adotoa;

    // ????
      
    // QUINTESSENCE
    //yprime[qidx] = quint->delQDot();
    //yprime[qidx +1] = quint->delQDot_dtau(adotoa, a, k2,hdot);


    //    cout << "sign: " << sign << endl;

    yprime[qidx] = DelQDot;
    yprime[qidx +1] = -2.0*adotoa*DelQDot -  (sign*a2*quint->Vprime2(quint->q(a),a,adotoa,tau) + k2)*DelQ -  0.5*hdot*quint->qDot(a); 
 
    

    //    yprime[qidx+1] *= sign;

    // dgtheta is 8pi G a^2 times  (rho + p) theta_species 
   
    dgtheta = Gpi8()*a2*( tau2rho_cdm(tau) * thetac + tau2rho_b(tau) * thetab
			  + 4.0/3.0 *(tau2rho_g(tau) * thetag + tau2rho_nu(tau)*thetar) 
			  + rho_nu0()*nuNR() *k * fnu/(a2*a2) );   
    double dgthetaQ =sign*Gpi8()*(  k2 * quintcosmos->tau2phidot(tau) * DelQ  );   // Quintessential correction     
        
    dgtheta += dgthetaQ;
   
    etadot = dgtheta * .5 / k2;
    yprime[1] = etadot;


    //  8*pi*G*(rho+P)*sigma*a**2. 
   
    dgshear = Gpi8()*a2*( 4.0/3.0*  (tau2rho_g(tau) * shearg + tau2rho_nu(tau)*shearr) 
       +   rho_nu0()*nuNR()/(a2*a2)  * shearnu); 

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
    if (!coupled) {   //  Treat baryons and photons as uncoupled.
      thetabdot = -adotoa * thetab + k2 * cs2 * deltab + pb43 * drag;
    } else {
      //  Treat baryons and photons as tightly coupled. 
      //  Zeroth-order approximation to baryon velocity. 
      thetabdot = (-adotoa * thetab + k2 * cs2 * deltab 
		   + k2 * pb43 * (deltag * .25 - shearg)) / (pb43 + 1.);
      //  (\ddot a)/a. 
      adotdota = (adotoa * adotoa - gpres) * .5;

     //  First-order approximation to baryon-photon slip, thetabdot-thetagdot. 
      slip = pb43 * 2. / (pb43 + 1.) * adotoa * (thetab - thetag) + 
	1. / opac(tau) * (-adotdota * thetab - adotoa * k2 * .5 * deltag 
		     + k2 * (cs2 * deltabdot - deltagdot * .25)) / (pb43 + 1.);
      //  First-order approximation to baryon velocity. 
      thetabdot += pb43 / (pb43 + 1.) * slip;
    }
    yprime[5] = thetabdot;

    //  Photon total intensity and polarization equations of motion. 
    yprime[6] = deltagdot;

    //cs2 =  tau2cs2(tau);
    



    thetagdot = (-thetabdot - adotoa * thetab + k2 * cs2 *
	     deltab) / pb43 + k2 * (deltag * .25 - shearg);
    yprime[7] = thetagdot;

    if (!coupled) {           //  Treat baryons and photons as uncoupled. 
	yprime[8] = 8.0/15.0 * thetag - 0.6* k * y[9]
		 - opac(tau) * y[8] + 4.0/15.0 * hdot + 8.0/5.0*etadot + opac(tau) * .1 * polter;
	//  Polarization equations for l = 0, 1, 2. 
	yprime[lmaxg + 7] = -k * y[ lmaxg + 8] - opac(tau) * y[ lmaxg + 7] + opac(tau) * .5 * polter;
	yprime[lmaxg + 8] = k / 3. * (y[lmaxg + 7] - y[ lmaxg + 9] * 2.) - opac(tau) * y[ lmaxg + 8];
	yprime[lmaxg + 9] = k * (y[ lmaxg + 8] * .4 - y[ lmaxg + 10] * .6) - opac(tau) * y[ lmaxg + 9] + opac(tau) * .1 * polter;
	for (l = 3; l < lmaxg; ++l) {
	    yprime[l + 6] = k * denl[l] * (l * y[ l + 5] - (l + 1) * y[ l + 7]) - opac(tau) * y[ l + 6];
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

	//yprime[l + 6] = yprime[lmaxg + 7 + l] = 0;

      }
    }
    
     //  Truncate moment expansion. 
    // 	yprime(8+lmax)=ak*lmax*y(7+lmax)/(2*lmax+1)-opac*y(8+lmax) 
    // 	yprime(9+2*lmax)=ak*lmax*y(8+2*lmax)/(2*lmax+1)-opac*y(9+2*lmax) 
    yprime[lmaxg + 6] = k * y[lmaxg + 5] - (lmaxg + 1) / tau * y[lmaxg + 6] - opac(tau) * y[lmaxg + 6];
    yprime[lmaxg2 + 7] = k * y[lmaxg2 + 6] - (lmaxg + 1) / tau * y[lmaxg2 + 7] - opac(tau) * y[lmaxg2 + 7];

    //  Massless neutrino equations of motion. 
    if (cosmos->nuR()>0) {
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


  if (cosmos->haveMassiveNeutrinos()) {
    MassiveNeutrinos::propagateMassiveNeutrinoMomentsSynchronous(y, yprime, a, tau,
                                   neutrinoMass_kbT, k, hdot, etadot);
  }
  Psi = -1.5 * dgshear / k2 - adotoa + eta - adotoa / ( 2 * k2 ) * ( hdot + 6 * etadot );
}



void QuintSynchronous::initialScalarPerturbations(const ControlPanel &control, const double tau)
{ 
     double  a, c, h;
     int l;
     double a2, f1, delta0;
     double deltac, deltab, deltag, thetab, thetac, deltan;
     double deltar, thetag, thetan, shearr, thetar,ahdot,eta;
     double psi, kt2;
     double DelQ=0,DelQDot=0,Epsilon; // quintessence perturbation

    //  Initial conditions. 

    a = tau2a(tau);
    a2 = a * a;
    k2 = k*k;
#ifdef OUTPUTFILE
    ofs = new ofstream("colddark.dat");
#endif
  
    double sum = tau2rho(tau) + tau2p(tau); // COM s = grho + gpres
    double  fracnu = (tau2rho_nu(tau) + tau2rho_nuNR(tau) + tau2p_nu(tau) + tau2p_nuNR(tau) )  / sum;
  
    // Use yrad=rho_matter/rho_rad to correct initial conditions for 
    // matter+radiation.
    double yrad = tau2rho_m(tau) /  ( tau2rho_g(tau) + tau2rho_nu(tau) + tau2rho_nuNR(tau));
    double Omega_nu =  cosmos->rho2omega(tau2rho_nu(tau)+tau2rho_nuNR(tau),tau);
      
    //  Choose one of the following four cases for initial conditions, or 
    //  add your own.  Comment out the other cases. 
  
    switch (control.initialCond) {
    case ControlPanel::adiabatic: 
      //-----------------------------------------------------------
      //  First case. 
      //  Isentropic ("adiabatic") initial conditions. 
      psi = -1 / (4*Omega_nu + 15.0)*10;      
      c = (fracnu * 4. + 15.) / 20. * psi;   // 4*fracnu is approx 1.5 so 15.0 still rules 
      kt2 = pow(k*tau,2);
      h = c * kt2 * (1. - yrad * .2);  // 1.0 rules over yrad*2 (=1e-4)  observe:   h is propto (k tau)^2
      eta = c * 2. - (fracnu * 4. + 5.) / 6. / (fracnu * 4. + 15.) * c * kt2 * (1. - yrad / 3.);
      f1 = (fracnu * 4. + 23.) / (fracnu * 4. + 15.);
      deltac = h * -.5;
      deltag = h * -2.0/3.0 * (1. - kt2 / 36.);
      deltab = deltac;  
      deltar = h * -2.0/3.0 * (1. - kt2 / 36. * f1);   
      thetac = 0.;
      thetag = -c / 18. * kt2 * kt2 / tau;
      thetab = thetag;
      thetar = f1 * thetag;
      shearr = k2 * 4.0/15.0 / (Gpi8()*a2*sum) * psi * (yrad * 7.0/36.0 + 1.);


      quint->touch(tau);
      // initial conditions from Baccigalupi et. al. 
      // unfortunatly, inconsistency suspected (will soon be resolved)
      Epsilon = quintcosmos->H_0_cpm() * sqrt(cosmos->omega_relativistic()) * tau;
      DelQ = -Epsilon*h*tau *quintcosmos->tau2phidot(tau) /(20.0 *a);
      DelQDot = 4.0/tau * DelQ;      

      DelQ = 0;
      DelQDot = 0;


      ahdot = c * 2. * k2 * tau * a * (1. - yrad * .3);
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
    
    // QUINTESSENCE (appended to the back of the array)
    y[qidx]    = DelQ;
    y[qidx+1] = DelQDot;
    // y[qidx+2] = ahdot;

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
    y[(lmaxg << 1) + 8] = deltar;
    y[(lmaxg << 1) + 9] = thetar;
    y[(lmaxg << 1) + 10] = shearr * 2.;

    for (l = 3; l <=lmaxnr ; ++l) y[(lmaxg << 1) + 8 + l] = 0.;

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

void QuintSynchronous::getReady(const ControlPanel& control) {
  quint = quintcosmos->quintessence(); 
  Synchronous::getReady(control);  // now call parent objects get ready
  qidx = nvar-1;
}


void QuintSynchronous::calcPerturbationNr(const ControlPanel &control) {
  nqmax=MassiveNeutrinos::qGridSize();
  //   Calculate number of equations 
  if (control.scalar) {
    nvar = ((lmaxg + 1) << 1) + 5;
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
  } else nvar = 0;

  if (control.tensor) nvart = (lmaxt << 1) + 4; else nvart = 0;

  nvar +=2; // two more for quintessence

  y = new double[nvar+1];
  yprime = new double[nvar+1];
  yt = new double[nvart+1];
  ytprime = new double[nvart+1];
}
