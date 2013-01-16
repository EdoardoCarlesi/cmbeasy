#include "perturbation.h"

#ifdef MACOSX_PANTHER
extern "C" int isnan(double); 
extern "C" int isinf(double); 
#endif 

Perturbation::Perturbation(Cosmos* c)
  : SeedFlag(false), hnext(0), hnext_tensor(0),
    ForceTightCouplingOverrulePhoton(false),
    ForceTightCouplingOverruleBaryon(false),
    inform(false), cosmos(c), y(0), yt(0),
    yprime(0), ytprime(0), nqmax(0), lmaxnu(0),
    warnOnEarlySource(true), k(0)
{
  GPI8 = c->Gpi8();

  for (int j = 1; j <= LMX0; ++j)  denl[j] = 1. / (double) ((j << 1) + 1);
}

Perturbation::~Perturbation() {
  delete[] y;
  delete[] yprime;
  delete[] yt;
  delete[] ytprime;
}

/*!
  Propagates the scalar perturbation. If the switch of the tight coupling regime to the full system
  happens within tau and tauend, it will split the integration into two pieces. One which
  is definitly only in the isTightCoupling() true region and one in which it is definitly false.
  The slight jump is of the order 10^-10 and doesn't cause any harm.

  It does so for criterias both for photons and baryons. The jump is determined beforehand
  in tauTightCoupling()

  Please note that this routine is supplemented by a much more elegant version for 
  the up to date speedyinvariant classes


*/
void Perturbation::propagateScalar(double *tau, const double tauend, const double precision) {
  Miscmath::keyword method= Miscmath::rkutta;
  bool original = false;
  k2 = k*k;
  if (hnext == 0) hnext = (*tau)*1e-1;
  if (ForceTightCouplingOverrulePhoton) throw Bad_Error("Perturbation::propagateScalar()  Safety alert: ForceTightCouplingOverrule should be switched off except\nfor propagateScalar(). Something else seems to set it. Please correct.");
  if (ForceTightCouplingOverruleBaryon) throw Bad_Error("Perturbation::propagateScalar()  Safety alert: ForceTightCouplingOverrule should be switched off except\nfor propagateScalar(). Something else seems to set it. Please correct.");
  // first, let us find out whether our interval contains a jump in the tight coupling for photons
  MustCoupleTightlyPhoton = isTightCoupling(*tau,photon);
  if (MustCoupleTightlyPhoton && (! isTightCoupling(tauend,photon)) ) { // if there is a jump in tightCoupling of photons
    pair<double,double> TauJump = tauTightCoupling(*tau, tauend,photon);
    //    cout << "PROPAGATESCALAR photons: k: " << k << "  " << TauJump.first <<   "    ::::    " << TauJump.second << "   === " << TauJump.first - TauJump.second << endl;
    ForceTightCouplingOverrulePhoton = true; // we couldn't do it before, cause if the above would have used isTightCoupling which would have been overruled
    Miscmath::odeint(y,nvar,*tau,TauJump.first,precision,hnext,0, (moDerivs)&Perturbation::fderivsWrapper, *this,original,method);
    *tau = TauJump.second;   // so we make a tiny jump, just 10^-10
    hnext = 1e-1;  // small intitial step  
    MustCoupleTightlyPhoton = false; // from now on, don't couple photons
  }
  
  // now, let us find out whether the remaining interval contains a jump in the tight coupling for baryons
  MustCoupleTightlyBaryon = isTightCoupling(*tau,baryon);
  if (MustCoupleTightlyBaryon && (! isTightCoupling(tauend,baryon)) ) { // if there is a jump in tightCoupling of baryons
    pair<double,double> TauJump = tauTightCoupling(*tau, tauend,baryon);
    //    cout << "PROPAGATESCALAR baryons: k: " << k << "  " << TauJump.first <<   "    ::::    " << TauJump.second << "   === " << TauJump.first - TauJump.second << endl;
    ForceTightCouplingOverruleBaryon = true; // we couldn't do it before, cause if the above would have used isTightCoupling which would have been overruled
    Miscmath::odeint(y,nvar,*tau,TauJump.first,precision,hnext,0, (moDerivs)&Perturbation::fderivsWrapper, *this,original,method);
    *tau = TauJump.second;   // so we make a tiny jump, just 10^-10
    hnext = 1e-1;  // small intitial step  
    MustCoupleTightlyBaryon = false; // from now on, don't couple baryons
  }
  
  // So no matter what, there will be something left to do...
    ForceTightCouplingOverrulePhoton = true; 
  ForceTightCouplingOverruleBaryon = true;  // what ever is set by us, the perturbations should take
  hnext = Miscmath::odeint(y,nvar,*tau,tauend,precision,hnext,0, (moDerivs)&Perturbation::fderivsWrapper, *this,original,method);
  *tau = tauend;
  ForceTightCouplingOverrulePhoton = false; 
  ForceTightCouplingOverruleBaryon = false;  // release it

  for (int i=1; i<=nvar; i++) {
    if (isnan(y[i])) {
      cout << endl << "y["<<i<<"] isnan at time: " << *tau << "  [ " << cosmos->tau2a(*tau) << " ]"<<endl;
      throw Bad_Error("isnan !");
    }
  }
}


//X void Perturbation::setRegimes(double tau) {
//X   k2 = k*k;
//X   if (ForceTightCouplingOverrulePhoton) throw Bad_Error("Perturbation::propagateScalar()  Safety alert: ForceTightCouplingOverrule should be switched off except\nfor propagateScalar(). Something else seems to set it. Please correct.");
//X   if (ForceTightCouplingOverruleBaryon) throw Bad_Error("Perturbation::propagateScalar()  Safety alert: ForceTightCouplingOverrule should be switched off except\nfor propagateScalar(). Something else seems to set it. Please correct.");
//X   // first, let us find out whether our interval contains a jump in the tight coupling for photons
//X   MustCoupleTightlyPhoton = isTightCoupling(tau,photon);
//X   MustCoupleTightlyBaryon = isTightCoupling(tau,baryon);
//X }



void Perturbation::propagateTensor(double *tau, const double tauend, const double precision) { 
  k2 = k*k;
  if (hnext_tensor == 0) hnext_tensor = (*tau)*1e-1;
  Miscmath::odeint(yt,nvart,*tau,tauend,precision,hnext_tensor,0, (moDerivs)&Perturbation::fderivsTensorWrapper, *this,false);
  *tau = tauend;
}


void Perturbation::initDlnf() {
 // dlnf_0(q) / dlnq  where f(q) = 1 / (exp(q) + 1)
  double dq = 1.;
  for (int i = 0; i < NQMAX0; ++i) {
    double q = i * dq + .5;
    dlfdlq[i] = -q / (exp(-q)  + 1.);
  }
}

/*! 
  Make sure to call this if you reimplement your
  own getReady() 
*/
void Perturbation::getReady(const ControlPanel& control) {
  setSeedFlag(control.initialCond == ControlPanel::isoSeed);
  initDlnf();
  calcPerturbationNr(control);
}

/*!
  determine whether we are in tight coupling or not.

 If ForceTightCouplingOverrule is true then we return MustTightlyCouple.
 This avoids unwanted jumps in fderivs during integration. 
 Please note that we always integrate up to the change in the
 tightcoupling regime. Then we make a little jump (in propagateScalar),
 then we start again in the other regime.
*/
   
bool Perturbation::isTightCoupling(const double tau, Species species) {
#define EP0 1.0e-2
  if (tau > 200) return false;
//X   double ep;
//X   if (k > epsw) ep = EP0; else ep = 0.5*EP0;
    
  //    double tcp1 = k / opac(tau);
  // double tcp2 = 1. / (opac(tau) * tau);
  
  double epsilon = max(k / opac(tau) , 1.0/(opac(tau)*tau));
  double threshold = 5e-3;

  switch (species) {
  case photon:
    if (ForceTightCouplingOverrulePhoton) return MustCoupleTightlyPhoton;
    //    if (tcp1 > ep || tcp2 > ep && tcp1 > 1e-4)  return false;
    if (epsilon > threshold && tau2a(tau) > 1.35e-5) return false;
    return true;
    break;
  case baryon:
    if (ForceTightCouplingOverruleBaryon) return MustCoupleTightlyBaryon;
    //epsilon *=0.5;
    epsilon /= (1.0 + cosmos->tau2R(tau));
    if (epsilon > threshold && tau2a(tau) > 1.35e-5) return false;
    //if (tcp1 > ep || tcp2 > ep && tcp1 > 1e-4)  return false; 
    return true;
    break;
  default:
    throw Bad_Error("Perturbation::isTightCoupling(): Neither photons nor baryons requested");
  }
}

/*!
  Returns a pair of double which sandwiches the switch of  isTightCoupling()
  from true to false. The first value is such that isTightCoupling() will return
  true, the second one will give false.
  This is used to integrate the perturbations in a better controlled way.
  The switch in the integrand usually is troublesome to any integration mechanism,
  for it is non-analytic. 
*/
pair<double, double> Perturbation::tauTightCoupling(double min, double max, Species species) {  
  double original = max - min;
  double step = original*0.5;
  double mid = min + step;

  for(;;) {
  //while (step/mid > 1e-16) {
    mid = min + step;
    //if (step/mid < 1e-2) break;
    /*
    cout << "TauTightCoupling: " << min << " ? ";
    if (isTightCoupling(min,species)) cout << " yes "; else cout << " no: ";
    cout << "  max: " << max << " ? ";
    if (isTightCoupling(max,species)) cout << " yes "; else cout << " no: ";
    cout << endl;
    */
    if (mid == min) break;  // meaning: as precise as double gives us
    if (isTightCoupling(mid,species)) min = mid; else max = mid;
    step *= 0.5;
  }
  return make_pair(min,max);
}


/*! 
  Since tensors are gauge invariant and independent of dark energy,
  we implement them here since v4.0 
  Of course, a re-implentation is always possible
*/
void Perturbation::initialTensorPerturbations() {	
  
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

void Perturbation::fderivsTensor(const double tau,   double *y, double *yprime) {
    double  psie, htpr, a;
    int l;
    double a2, htdpr, rhonu, ep, ht, adotoa;

    double cs2;
    double deltap0, deltat0, tcp, pnu;
    int ind1, ind2,ind3;
    double tcp1, tcp2;

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

    if (cosmos->nuNR() > 0) {
      static bool didWarn = false;
      if (!didWarn) {
        didWarn=true;
        cout << "warning: ignoring tensor contribution of massive neutrinos\n";
      }
      //MassiveNeutrinos::nu1(a, &rhonu, &pnu);
    }

    //  8*pi*G*rho*a**2 and 8*pi*G*P*a**2.
    //  Tensors
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

void Perturbation::tensorSources(double tau, double *dt, double *dte, double *dtb) {
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
