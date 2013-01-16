#include "cncosmos.h"
#include "exponential.h"  // should be removed once history is more generic
#include "massiveneutrinosnew.h"

CnCosmos::CnCosmos(Quintessence::Type quinttype)
   : QuintCosmos(quinttype)
{
  mNeutrinoCoupling = 0.;
}

CnCosmos::~CnCosmos()
{
}
 
/*! 
  Implementation of the evolution equations - input are the current values of the
  background variables in the double array y, and the double array dy needs
  to be filled with their respective derivatives at the time tau
 
  As always, output of tau is in Mpc, i.e. d/dtau is in Mpc^-1 
*/
void CnCosmos::propagateHistoryInTau(const double tau,const double* y,double* dy)
{
  double a = y[1];
  double rho_g = y[2];
  double rho_nu = y[6];
  double rho_b = y[3];
  double rho_c = y[4];

  double phi = y[8];
  quint->setQ(y[8]);
  quint->setQDot(y[9]);
  double rho_nuNR = 0.;
  double p_nuNR = 0.;
  if (haveMassiveNeutrinos()) {
    double rnu, pnu;
    MassiveNeutrinosNew::nu1(a *avarmnu(phi), &rnu, &pnu);
    rho_nuNR = rho_nu0() * nuNR() * rnu/(a*a*a*a);
    p_nuNR = rho_nu0() * nuNR() * pnu/(a*a*a*a);
 //TauNuMass = new Spline(Tau2A,"Nu_Mass",&SplineAnchor);
  //TauNuMass->arm();
  }

  //X   rho_nuNR = nuNR() * massiveNeutrinoRho(0, A);
  //X   p_nuNR = nuNR() * massiveNeutrinoPressure(0, A);
  const double rhoq = quint->rho(a);
  double totalRho = rho_g + rho_c + rho_b + rho_nu  + rho_nuNR + rho__v + rhoq;

 /* cout << "rho_g: " << rho_g << " rho_c: " << rho_c << " rho_b: " << rho_b << " rho_nu: " << rho_nu << endl;
  cout << "rho_nuNR: " << rho_nuNR << " rho__v: "<< rho__v << endl;
  cout << "rhoq: " << rhoq << " totalRho: " << totalRho << endl;
  cout << "omegaq: " << rhoq/totalRho << "  at a = " << a << endl;
  */

  double H = a* sqrt(1.0/3.0*totalRho)/M_p();

  if (isnan(H)) {
    cout << "isnan H: " << tau << "a: " << a << endl;
    cout << "quintrho: " << quint->rho(a) << endl;
    throw Bad_Error("quintcosmos isnan H  !");
  }

  dy[1] = H*a;      // da/dtau is H*a
  dy[2] = -4.0*H*rho_g;  // gammas
  dy[3] = -3.0*H*rho_b;  // baryons
  dy[4] = -3.0*H*rho_c;  // cold dark matter


  if (quint->needsPropagation()) {
    const double massNu_kbt = mass_nu_kbT(phi);
    const double dlnmdphi = massNu_kbt?dMass_nu_kbTdphi(phi)/massNu_kbt:0.;
    //cout << mass_nu_eV(phi) << "   ?   " << dMass_nu_kbTdphi(phi) << "  ?  " << dlnmdphi << endl;
    dy[8] = quint->qDot(a);  // quintessence
    dy[9] = -2*H*quint->qDot(a)-a*a*quint->Vprime(quint->q(a),a,H)-a*a*dlnmdphi*(rho_nuNR-3.*p_nuNR)/M_p();     // quintessence
     //cout << a <<  "  " << (phi/m_p()) << " " << mass_nu_kbt(phi) << " "  << dmass_nu_kbtdphi(phi) << endl;
     //cout <<  dlnmdphi << "  /  " << (rho_nuNR-3.*p_nuNR) << endl;
     //cout << "test: "  << " " << (phi/m_p()) << " " << mass_nu_kbt(-1.2528*m_p()) << " "  << dmass_nu_kbtdphi(-1.2528*m_p()) << endl;
  } else {
    dy[8] = 0;
    dy[9] = 0;
  }
  dy[5] = a;   // dt/dtau
  dy[6] = -4.0*H*rho_nu;  // massless neutrinos

  double r = 3.0/4.0 * rho_b/rho_g;
  double cs = 1.0/sqrt(3*(1+r));  // sound speed
  dy[7] = cs;  // integrate for sound - horizon
}

void CnCosmos::history(bool inform)
{
  if (validHistory()) throw Bad_Error("CnCosmos::history() I already have a history! Call reset() first");
  getReady();
  quint->prepare();
  rho__v = omega2rho(omega_v());

  cross_W25 = -4.6;

  double y[10];
  double ydot[10];
  for (int k = 0; k < 9; k++) { y[k]=0; ydot[k] = 0;}

  y[1] = initialScaleFactor();
  y[2] = initialRho_g(y[1]);       
  y[6] = initialRho_nu(y[1]);   
  y[3] = initialRho_b(y[1]); 
  y[4] = initialRho_cdm(y[1]);
  y[5] = 0;
  y[7] = 0;

  double atmp = y[1];
  double rhonuNR = 0.;
  double pnuNR = 0.;
  if (haveMassiveNeutrinos()) {
    MassiveNeutrinosNew::initnu1(amnu);
    double rnu, pnu;
#warning assuming nu mass at start of history is the uncoupled one (can't use phi here, because I need rhonuNR to put phi_initial on attractor below)
    const double massNukbT = mass_nu_kbT_uncoupled();
//X     const double massNukbT = avarmnu(q);
    MassiveNeutrinosNew::nu1(atmp*massNukbT, &rnu, &pnu);
    //TauNuMass->set();
    rhonuNR = rho_nu0() * nuNR() * rnu/(atmp*atmp*atmp*atmp);
    pnuNR = rho_nu0() * nuNR() * pnu/(atmp*atmp*atmp*atmp);
  }

//X   y[8] = q;
//X   y[9] = quint->qDot(y[1]);
#warning setting q'_initial from KG, but ignoring phi-contrib to rho and therefore H, and assuming beta=0
#warning fix varying beta
  double beta = mNeutrinoCoupling;
  Exponential *exp = dynamic_cast<Exponential*>(quintessence());
  if (beta != 0 && !exp)
    throw Bad_Error("CnCosmos::history - only implemented for exponential potential");
  if (exp)  {
    beta = 0;
    const double lambda = exp->lambda();
    const double Omega_potRDE = (4.-3.*beta*(lambda+beta))/(3.*pow(lambda-beta, 2));
    const double Omega_rRDE = (4.-beta*(lambda+beta))/(3.*pow(lambda-beta, 2));
//X     cout << "Omega_RDE = " << Omega_rRDE << " and Omega_potRDE= " << Omega_potRDE << endl;
    double totalRhotmp = y[2]+y[6]+y[3]+y[4];
    const double H = atmp* sqrt(1.0/3.0*totalRhotmp)/M_p();
    double q = -M_p()/lambda*log(Omega_potRDE*totalRhotmp/M_p(4));

    static const bool setQPrimeFromKleinGordon = false;
//X     static const bool setQPrimeFromKleinGordon = true;
    const double S = beta * (rhonuNR-3.*pnuNR)/M_p();
    double qprime =  -atmp*atmp*quint->Vprime(q,atmp,H) - atmp*atmp*S;
    qprime /= 2.*H;
    if (setQPrimeFromKleinGordon) {
      cout << "Setting initial q' from KG - equation at a=" << atmp << endl;
      exp->setInitialQ(q);
      exp->setInitialQDot(qprime);
      cout << "setting initial q to " << (q/M_p()) << endl;
      cout << "setting initial q' to " << (qprime/M_p()) << " from KG - equation" << endl;
      cout << "Omega_potRDE was: " << Omega_potRDE << endl;
    } else if (true) {
//X       cout << "not setting initial q' from KG - equation" << endl;
//X       cout << "would have been: " << q/M_p() << " and " << qprime/M_p() << endl;
//X       cout << "for beta=" << beta << " and lambda= " << lambda << endl;
    }

//X     cout  << "q: "; cin >> q; cout << endl;
//X     cout << "qprime: "; cin >> qprime;
//X     exp->setInitialQ(q*M_p());
//X     exp->setInitialQDot(qprime*M_p());
  }


  quint->setInitial(y[1]);  // Inform quintessence to set its conditions
  y[8] = quint->q(y[1]);
  y[9] = quint->qDot(y[1]);

  //cout << "setting initial q to " << (y[8]/M_p()) << endl;
  //cout << "setting initial q' to " << (y[9]/M_p()) <<  endl;

  if (quint != 0 && inform) quint->printStatus();
  // now we want to get da/dtau in order to estimate the beginning tau....
  propagateHistoryInTau(0,y,ydot);  

  double tau = y[1]/ydot[1];
  
  if (inform) cout << "from initial step, we have: " << tau << endl;

  double dtau = tau*1e-1;     // and dtau is one tenth of this..
  // zeroing tau__equ and tau__ls as they will now be (re)calculated
  tau__equ = 0;
  tau__ls = 0;
 
  
  //  ofstream monitor("monit.dat");

  double A=0;
  double lastA = 0;
  double totalRho,totalPressure;
  Background b;
  double hnext = dtau*1e-12; // 1e-3;
  while (y[1] < a_0()) {   
    hnext = Miscmath::odeint(y,9,tau,tau+dtau,1e-9,hnext,0,(moDerivs)&QuintCosmos::propagateHistoryInTauWrapper, *this);
    tau += dtau;
    
    propagateHistoryInTau(tau,y,ydot);     // We would like to have the tau derivatives for splines
    b = fillHistorySplines(tau,y,ydot);         // this also updates quintessence to the right q and qdot
    lastA = A;
    A = b.A;
    Tau2Rho_qpot->set(quint->V(quint->q(A), A));
    Tau2Rho_qkin->set(quint->rho_kin(A));
    Tau2Rho_q->set(quint->rho(A));
    Tau2Phi->set(quint->q(A));
    Tau2PhiDot->set(quint->qDot(A));
  //      TauNuMass->set(mass_nu_eV(y[8]));

    // total rho and total pressure
    totalRho = b.rho_g + b.rho_nu + b.rho_nuNR + b.rho_c + b.rho_b  + rho__v + quint->rho(A);
    totalPressure = 1.0/3.0*(b.rho_g + b.rho_nu) + b.p_nuNR - rho__v + quint->p(A);
   
    // A2W.set(A,quint->w(A));

    Tau2W_q->set(quint->w(A));  
    Tau2Omega_q->set(quint->rho(A) / totalRho);
    Tau2Omega_qw ->set(quint->rho(A) / totalRho*quint->w(A));

    Tau2Vprime->set(quint->Vprime(quint->q(A),A,ydot[1]/A));
    Tau2VprimeOverRho->set(quint->Vprime(quint->q(A),A,ydot[1]/A)/quint->rho(A));

    // set eqn of state etc as function of logA (for Jan), also note crossing of 
    // w_q = -0.25 if a > 0.01...
    LogA2Omega_q->set(log(A),quint->rho(A) / totalRho);
    LogA2Omega_qw ->set(quint->rho(A) / totalRho / quint->w(A));
    if (quint->w(A) < -0.25 && cross_W25 == -4.6 && log(A) > -4.6)  cross_W25 = log(A); 

    Tau2Rho->set(totalRho); 
    Tau2P->set(totalPressure); 
    Tau2P_q->set(quint->p(A));
    
    // 8*pi*G*a^2*rho
    Tau2GRho->set(Gpi8()*A*A*totalRho);
    
    A2Omega_q->set(quint->rho(A) / totalRho);
    A2Omega_qw ->set(quint->rho(A) / totalRho*quint->w(A));
    
    Z2W->set(quint->w(A));

    for (int i=1; i<=7; i++) {
      if (isnan(y[i])) {
	cout << endl << "quintcosmos:: y["<<i<<"] isnan at time: " << tau << "  [ " << A << " ]"<<endl;
	throw Bad_Error("quintcosmos isnan !");
      }
    }

  
    // determine tau__equ by comparing densities
    // this is rough and should is refined by really narrowing in on
    // the zeros of the equality later
    if (b.rho_b + b.rho_c > b.rho_g && tau__equ == 0) tau__equ = tau;
 
    // determine t_LS by noting the crossing of a above the a_ls level...
    // rough, as above,  is refined in a later step by looking for
    // the interpolation crossing of a(tau) >= a_ls()
    if (A>= a_ls() && tau__ls == 0) tau__ls = tau;
   
    if (A < 0.1*a_0()) dtau *= 1.01;
    if (fabs(1.-lastA/A) > 0.05) dtau *= 0.5;
//X     cout << "hist: " << A << "  " << tau << " dtau: " << dtau <<  endl;
    dtau = min(10.0,dtau);
  }
  

  Tau2A->arm(Spline::thoseReady);
  splineT2Tau->arm(); 
  //  splineA2Tau->dump("a2tau");
  splineA2Tau->arm(Spline::thoseReady);
  LogA2Omega_q->arm(Spline::thoseReady);
   
  Z2W->flip(); // we flip this cause we do not want 10^6 ... 0  but the other way around (right order)
  Z2iH_spline->flip();
  
  Z2iH_spline->arm(Spline::thoseReady);

  //calculate distance stuff
  double dzlum = distzmax / double(ndistpoints-1);
  for (double ztemp = 0; ztemp < distzmax; ztemp += dzlum)
    Z2dist->set(ztemp, Z2iH_spline->integrate(0,ztemp));
  Z2dist->arm(Spline::thoseReady);

  // if inform is true, dump various background quantities to files
  if (inform) {
    Z2iH_spline->dump("inverseH",false);
    Tau2A->dump("his");
    Tau2ADot->dump("adot");
    Tau2Rho->dump("rho");
    Tau2Rho_g->dump("rho_g");
    Tau2Rho_b->dump("rho_b");
    Tau2Rho_nu->dump("rho_nu");
    Tau2Rho_nuNR->dump("rho_nuNR");
    Tau2P_nuNR->dump("p_nuNR");
    Tau2w_nuNR->dump("w_nuNR");
    Tau2wdot_nuNR->dump("wdot_nuNR");
    Tau2Rho_cdm->dump("rho_cdm");
    Tau2Rho_m->dump("rho_m");
    
    Tau2Rho_qpot->dump("rho_qpot");
    Tau2Rho_qkin->dump("rho_qkin");
    Tau2Rho_q->dump("rho_q");
    Tau2P_q->dump("pres_q");
    Tau2P->dump("pres");
    Tau2W_q->dump("wq");
    Tau2Omega_qw->dump("woq");
    Tau2Omega_q->dump("omega_q");
  }
  findSpecialMoments();
 
  Tau2Rho->derive(*Tau2RhoDot);
  Tau2RhoDot->arm();   
  if (inform) Tau2RhoDot->dump("rhodot");
   
  Tau2Rho_g->derive(*Tau2RhoDot_g);
  Tau2RhoDot_g->arm();   

  Tau2Rho_nu->derive(*Tau2RhoDot_nu);
  Tau2RhoDot_nu->arm();
  
  Tau2Rho_nuNR->derive(*Tau2RhoDot_nuNR);
  Tau2RhoDot_nuNR->arm();

  
  Tau2w_nuNR->derive(*Tau2wdot_nuNR);
  Tau2wdot_nuNR->arm();

  Tau2Rho_m->derive(*Tau2RhoDot_m);
  Tau2RhoDot_m->arm();   

  Tau2GRho->derive(*Tau2GRhoDot);
  Tau2GRhoDot->arm();   

  Tau2Vprime->derive(*Tau2DotVprime);
  Tau2DotVprime->arm();
  
  Tau2VprimeOverRho->derive(*Tau2DotVprimeOverRho);
  Tau2DotVprimeOverRho->arm();
  
  Tau2W_q->derive(*Tau2WDot_q);
  Tau2WDot_q->arm();   
 
   Z2W->derive(*Z2dWdZ,Spline::all);

   if (inform) Z2dWdZ->dump("dwdz");
   
   for (int k = 0; k < Z2dWdZ->size(); k++) {
     double z = Z2dWdZ->x(k);
     double tau = z2tau(z);  // building blocks for the distinguishing of fast and slow quintessence w-change
     Z2OmegaWZ->set(pow(Z2dWdZ->y(k),2) * tau2rho_q(tau)/tau2rho(tau) *z *z);
   } 

   Z2OmegaWZ->arm();
   
   ValidOmega = true; // new  
   quint->postPrepare(); // maybe some things need to be done with the knowledge of history ?
   if (quint->needsToUpdatePhi()) { // uups, is Phi not the real one ? (only for Arbitrary, so far)
     //     Tau2Phi->dump("phibefore");
     Tau2Phi->disarm();
     Tau2PhiDot->disarm();
     for (int i=0; i < Tau2Phi->size(); i++) {
       double tau=Tau2Phi->x(i);
       double a = tau2a(tau);
       Tau2Phi->setY(i,quint->q(a));
       Tau2PhiDot->setY(i,quint->qDot(a));
     } 
     Tau2Phi->arm();
     Tau2PhiDot->arm();
    
   }
   if (inform)
     printStatus();
}

Background CnCosmos::fillHistorySplines(double tau, const double *y, const double *ydot, bool writeSplines)
{
  Background b;
  b.A = y[1];
  b.rho_g = y[2];
  b.rho_nu = y[6];
  b.rho_b = y[3];
  b.rho_c = y[4];
  b.t = y[5];
  b.rs = y[7];
  const double phi = y[8];
  double A = y[1];
  
  //X   b.rho_nuNR = nuNR() * massiveNeutrinoRho(0, A);
  //X   b.p_nuNR = nuNR() * massiveNeutrinoPressure(0, A);
  
  double plain =0;
  double rnu=0, pnu=0;
  if (haveMassiveNeutrinos() !=0.0) {
    MassiveNeutrinosNew::nu1(A *avarmnu(phi), &rnu, &pnu) ;
    //cout << A << " - " << avarmnu(phi) << endl;
    plain = rho_nu0() * nuNR() /(A*A*A*A);
  }
  
  b.rho_nuNR = rnu * plain;
  b.p_nuNR = pnu*plain;
  
  if (writeSplines)
  {
    double rho_m = b.rho_b + b.rho_c;
    // scale factor and its derivative
    Tau2A->set(tau,A);
    Tau2ADot->set(ydot[1]);
    Tau2ADotToA->set(ydot[1]/y[1]);
    
    // the photons and massless and massive (NonRelativistic) neutrinos
    Tau2Rho_g->set(b.rho_g);
    Tau2Rho_nu->set(b.rho_nu);
 
   Tau2Rho_nuNR->set(b.rho_nuNR);
  //au2Rho_nuNR->set();
	
    Tau2P_nuNR->set(b.p_nuNR);
    Tau2w_nuNR->set(b.p_nuNR/b.rho_nuNR);
    //TauNuMass->set(mass_nu_eV(y[8]),tau);

    // baryons, cdm and total matter
    Tau2Rho_b->set(b.rho_b);
    Tau2Rho_m->set(rho_m);
    Tau2Rho_cdm->set(b.rho_c);
    
    // 4/3 of the ratio of gammas over baryons (important for soundspeed)
    Tau2R->set(4.0/3.0*b.rho_g/b.rho_b);
    
    // getting from a <-> tau , t <-> tau
    splineA2Tau->set(A,tau);
    splineT2Tau->set(b.t,tau);
    
    // t(tau), soundhorizon(tau) and soundspeed^2(tau)
    Tau2T->set(b.t);
    Tau2Rs->set(b.rs);
    Tau2Cs2->set(ydot[7]*ydot[7]);
    Tau2Cs->set(ydot[7]);
    
    Z2iH_spline->set(a2z(A), A*A / ydot[1]);
  }
  
  return b;
}

void CnCosmos::printStatus(ostream &o)
{
	double tt1, tt2;

	tt1 = 3.e+16*pow(mpc2s(),-1)*tau2t(1744.3475);
	tt2 = 3.e+16*pow(mpc2s(),-1)*tau2t(1744.1854);

		cout << " tt1: " << tt1 << " tt2: " << tt2 << " tt1-tt2: " << tt1-tt2 << endl;

  if (quint) quint->printStatus();
  o.setf(ios::scientific);
  o << endl << endl;
  o << " ==============================================================="  << endl;
  o << " ======================= CnCosmos Status =======================" << endl;
  o << " ==============================================================="  << endl;
  o << " h : " << h() << endl;
  o << " Hubble's constant: " << H_0_cpm() << " Mpc^-1" << endl;
  o << " Hence hubble time is : " << 1.0 / (   cpm2sInv()*tau2Hubble(tau_0()) ) << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " omega_phot: " << omega_g() <<      "    rho_g0:    " << rho_g0() << endl;
  o << " omega_nu  : " << omega_nu() <<     "    rho_nu0:   " << rho_nu0() << endl;
  o << " omega_nuNR: " << omega_nuNR() <<   "    rho_nuNR0: " << rho_nuNR0() << endl;
  o << " omega_cdm : " << omega_cdm() <<    "    rho_cdm0:  " << rho_cdm0()  << endl;
  o << " omega_b   : " << omega_b() <<      "    rho_b0:    " << rho_b0() << endl;
  o << " omega_quint: " << omega_q() <<      "    rho_q0:    " << tau2rho_q(tau_0()) << endl;
  o << " omega_kurv.   : " << omega_k() << endl;
  o << " --------------------------" << endl;
  o << " omega_0:      " << omega_0() << "  (the user-defined value)" << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " omega_nuNR h^2: " << omega_nuNR()*h2() << endl;
  o << " omega_cdm  h^2: " << omega_cdm()*h2() << endl;
  o << " omega_b    h^2: " << omega_b()*h2() <<  endl;
  o << " omega_m    h^2: " <<  omega_m()*h2() << endl;
  o << " -------------------------------------------------------------" << endl;
  o << " uncoupled sum of neutrino masses (in eV) : " << nuNR()*mass_nu_eV_uncoupled() << endl;
  o << " neutrino coupling = " << mNeutrinoCoupling << endl;

  if (validHistory()) {
    o << " today's sum of neutrino masses (in eV) : " << nuNR()*mass_nu_eV(tau2phi(tau_0())) << endl;
    o << " -------------------------------------------------------------"  << endl;
    o << " ---------------conformal times & scale factors  -------------" << endl;
    o << " -------------------------------------------------------------"  << endl;
    o << " equality::        tau = "  << tau_equ() << "   a = " << a_equ() << endl;  
    o << " last scattering:: tau = "  << tau_ls() << "   a = " << a_ls()  << endl;  
    o << " today::           tau = "  << tau_0() << "   a = " << a_0() << endl;  
    o << " z_ls: " << z_ls() << endl;
   
    o << " -------------------------------------------------------------"  << endl;
    o << " optdlss: " << optDistanceLss() << endl; 
    o << " reionizationFraction " << reionizationFraction() << endl;
    o << " -------------------------------------------------------------"  << endl;
    o << " sound horizons: " << endl;
    o << " -------------------------------------------------------------"  << endl;
    o << " equality: " << tau2rs(tau_equ()) << " Mpc" << endl;
    o << " ls      : " << tau2rs(tau_ls()) << " Mpc" << endl;
    o << " -------------------------------------------------------------"  << endl;
 
    o << " time equ:   " << mpc2year()*tau2t(tau_equ())  << " years" << endl;
    o << " time ls:    " << mpc2year()*tau2t(tau_ls()) << " years" << endl;
    o << " time today: " << mpc2year()*tau2t(tau_0()) << " years" << endl;
    o << " ------------------------------------------------------------"  << endl;
   
    o << "::::  w_0: " <<  tau2w_q(tau_0()) << endl;
    o << "::: Omega_q average till ls: " << tau2AvOmega_q(tau_ls()) << endl;
    o << "::: Omega_q average till 0: " << Tau2Omega_q->average(tau_0()) << endl;
    o << "::: Omega_q structure f. : "<<  omesf(1.0/3.0) << endl; 
    o << endl << endl;
  }
}

void CnCosmos::getReady() {
  //     Initialize neutrino mass and storage.
  if (nuNR() == 0. || mass_nu_eV_uncoupled() == 0.) {
    amnu = 0.;
  } else {
    amnu = mass_nu_kbT_uncoupled();
  }
//X   cout << "Set up uncoupled neutrino mass: " << amnu << " - in eV: " << mass_nu_eV_uncoupled() << endl;
}

bool CnCosmos::massiveNeutrinosAreRelativistic(const double tau)
{
  return (mass_nu_cpm(Tau2Phi->fastY(tau)) < kbT(T_nu()/tau2a(tau)));
}

