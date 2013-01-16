#ifndef cOSMOS_H
#define cOSMOS_H 

#include "global.h"
#include "spline.h"
#include "mass_function.h"
#include <iostream>

#include "anchor.h"
#include "safevector.h"
#include "controlpanel.h"
#include "cl.h"
#include "gauge.h"
#include "basecosmos.h"

#define MAX_SOUNDHORIZON 10000
#define MAX_FAKTOR 10000



#define NTHERMO 10000  //! From CMBFAST

class SplineWeb;

/*!
  Small struct for passing data of background 
  quantities. This is only used by history() and
  fillHistorySplines(). The reason for it to exist is
  that a class derived from Cosmos does not need
  to replicate the code for filling some basic
  splines with background data. 
 */
struct Background {
  double A, rho_g, rho_nu, rho_b;
  double rho_c, t, rs;
  double rho_nuNR, p_nuNR;
};


class Recombination;
class Perturbation;

/*!
  General class for cosmological background evolution. Also
  features conversion and many convenience functions. 

  \code All quantities are  in powers of Mpc. \endcode 

  A cosmos consists of some parameters, like Omega_baryon, the Hubble parameter 
  today and so on and so on.

  The background evolution is calculated by calling history(). 

  Access to the background quantities is particularily simpe by  using functions that
  are named  
  \code
  y = x2y(x);
  \endcode

  where some quantity y is calculated as a function of x. 

  An example may clarify the use:
  \code
  #include "cosmos.h"

  Cosmos mine;                // creates a cosmos called mine
  mine.setOmega_v(0.7);   // sets vaccum Omega today to 0.7
  mine.setOmega_b(0.05); // sets baryon Omega to 0.05
  mine.setOmega_cdm(0.25); // sets cold dark matter to 0.25
  mine.seth(0.65);            // and a small Hubble h of 0.65, i.e. H = 100*h km/(s*Mpc)

  mine.history();               //  calculates the background evolution

  double lastscattering = mine.tau_ls(); // special moments like equality, ls 
  double today = mine.tau_0();   // and today (_0) are readily accessible

  double a = tau2a(tau_ls());    // If you'd like to have the scale-factor at ls
  double a = a_ls();                 // that's the same

  double zls = a2z(a);             // this is the redshift of lastscattering
  double zls = z_ls();               // the same :-)

  double age = tau2t(tau_0());          // this converts conformal time tau to usual time t

  cout << "The age of the universe in seconds: " << mpc2s()*age << endl;
  cout << "and in years: " << mpc2year()*age << endl;
  \endcode

  As you can see, convenience functions tau2XX exist for many quantities.
  They are calculated in history() and stored in Splines. The access is particularily
  fast, and also cached, so if you access, say tau2a(tau) twice, you will get the
  second call almost for free :-)

  The spline behind tau2a() is the mother (see class Spline) of all tau2XX()
  Splines. So if you do not mind, call this function first, if you have to call
  tau2XX for some quantities of the same tau, as this speeds up the calculation
  of all other tau2XX Splines.

  For examples, see fderivs() in invariant.cc.
*/




class Cosmos : public baseCosmos {
 public:  // cmbfast neutrino stuff
   double amnu;

 protected:
  static const double distzmax;// Maximum redshift luminosity distance calculated to
  static const int ndistpoints;// Number of points in luminosity distance

  Spline* Z2iH_spline; //!< used for calculating luminosity distance  
  double mInitialScaleFactor; //!< at which scale factor to start the integration of the background equations

 public:

  Spline* Z2H_spline; 
  Spline *splineA2Tau;
  Spline *splineT2Tau;

  double Omega_b;
  double Omega_cdm;
  double Omega_vacuum;

  double Omega_nu;
  double Omega_nuNR;
  double Omega_curv;
  double Omega_g;
  double Rho_crit;   //!< Today's Critical energy density set by findSpecialMoments()
  double h_h;              //!< hubbles small h
  double Y_P;              //!< primordial helium mass fraction
  double T_cmB;
  // double f_Nu; out of order until further restructuring

  double tau__equ, tau__ls, tau__0, tau__min, tau__max; //!< history() initializes these special values of the conformal time tau. Conformal time in [Mpc]

  double t__equ, t__ls, t__0, t__min, t__max; //!< history() initializes these special values of the usual time t. time in [Mpc]
  double a__equ;
  double z__ls;
  double rho__v; //!< vaccum energy set by history() excl. Quintessence

  double NrRelativisticNu;    //!< Nr of relativstic neutrinos
  double mNrNonRelativisticNu;
  double mMassNu_eV; //!< mass of a massive nu

  double Optdlss; // optdlss
  double ReionizationFraction;  // rif
  double ReionizationZ;   // zri
  double ReionizationZStop; //

#ifndef PRERELEASE
  double DarkmatterDipole;  // [units] = e*cm
  double DarkMatterMass; // in GeV
#endif

  bool ValidThermo; //!< true if thermal history has been calculated
  bool pkBaryon;
  Spline *splineCs2;
  Recombination *recombination; //!< This is recfast and co.
  vector<double> RecombinationParameter; //!< only needed for some varying alpha theories

#define CACHESPLINE(a) Spline* a; double a##X, a##Y      /* a define for the spline variables for convenience: one refernce and two caching double variables */


 CACHESPLINE(Visibility);
 CACHESPLINE(DVisibility);
 CACHESPLINE(DDVisibility);
 CACHESPLINE(Tau2A);
 CACHESPLINE(Tau2ADot);
 CACHESPLINE(Tau2ADotToA);
 CACHESPLINE(Opac);
 CACHESPLINE(D2Kappa);
 CACHESPLINE(Expmmu);
 CACHESPLINE(FreeElectronFraction);

 CACHESPLINE(Tau2Rho);
 CACHESPLINE(Tau2Rho_g);
 CACHESPLINE(Tau2Rho_nu);
 CACHESPLINE(Tau2Rho_nuNR);
 CACHESPLINE(Tau2Rho_m);
 CACHESPLINE(Tau2Rho_b);
 CACHESPLINE(Tau2Rho_cdm);

 CACHESPLINE(Tau2RhoDot);
 CACHESPLINE(Tau2RhoDot_b); 
 CACHESPLINE(Tau2RhoDot_cdm); 
 CACHESPLINE(Tau2RhoDot_m); 
 CACHESPLINE(Tau2RhoDot_g);
 CACHESPLINE(Tau2RhoDot_nu);
 CACHESPLINE(Tau2RhoDot_nuNR);


 CACHESPLINE(Tau2GRho);
 CACHESPLINE(Tau2GRhoDot);

 CACHESPLINE(Tau2P_nuNR);
 CACHESPLINE(Tau2w_nuNR);
 CACHESPLINE(Tau2wdot_nuNR);

 CACHESPLINE(Tau2P);

 CACHESPLINE(Tau2R);
 CACHESPLINE(Tau2T);
 CACHESPLINE(Tau2Rs);
 CACHESPLINE(Tau2Cs2);

#undef CACHESPLINE

 Spline* Tau2Cs; //!< Only needed for av. cs calculation, use tau2cs2() to get the soundspeed^2 as cached spline !

 Spline* Z2dist; //!< distance measure spline

 public:
 SafeVector<double> InitialPower; //!< The spectral index of initial perturbations. You may specify a whole bunch and it is *public*
 SafeVector<double> InitialPower_dnsdlnk; //!< d n_s / d lnk, the running spectral index,if you like one. Will be set to zero by default and in setInitialPower()
 SafeVector<double> InitialTensorPower; //!< The corresponding initial spectrum for tensor fluctuations
 // SafeVector<double> TensorRatio; //!< The ratio of tensor to scalar quadrupole   OBSOLTET

 vector<double> sigma8; //!< for each member of IntiailPower, this is designed to hold the sigma8 value. CobeNormalize() fill this
 SafeVector<double> PowerNormalization; //!< for each  member of IntiailPower, store CobeNormalize() CMB normaliziation of the power-spectrum 

 private:
 SplineWeb*  Power_CDM; //!< The splineweb containing the cdm powerspectrum 
 SplineWeb*  Power_b; //!< The splineweb containing the baryon powerspectrum
 SplineWeb*  Power_g; //!< The splineweb containing the photon powerspectrum
 SplineWeb*  Power_nu; //!< The splineweb containing the neutrino powerspectrum
 SplineWeb*  Power_nuNr; //!< The splineweb containing the massive neutrino powerspectrum
 SplineWeb*  Power_psi;

 protected:
 std::vector<double> mExtraTimeSteps;

 public:
 Cosmos();  //!< standard constructor
 virtual ~Cosmos() {};  
 virtual void reset(); //!< resets() cosmos. Frees all Splines and SplineWebs and calls initSplines()
 virtual void initSplines();  //!< Allocates Splines and SplineWebs
 Anchor SplineAnchor; //!< this is the cosmos' Anchor it takes care of the splines and is asked by reset() to get rid of the all splines to avoid memory leaks
 Anchor ThermoAnchor; //!< convenience anchor. currently only used by recombination. 
 double BeginRecombination; //!< Tau of the start of recombination

 SafeVector<double> TimeStep; // !< stores the conformal times at which sources should be computed

 void setInitialScaleFactor(const double a) { mInitialScaleFactor=a; }
 double initialScaleFactor() const { return mInitialScaleFactor; }

  double cpm2msun();
 double mpc2s(int n =1) const; //!< convert Mpc to second 
 double mpc2year() const { return  mpc2s()/(365*24*3600);}
 double cpm2sInv(int n=1) const ; //!< convert Mpc^-1 to second^-1
  double s2mpc(double s) const { return cpm2sInv()*s;} //!< convert second to Mpc
  double sInv2cpm(double s) const { return mpc2s()*s; }
  double mpc2veG(double m) const { return m*1.5637e38;}
  double veG2mpc(double v) const { return v/1.5637e38; }  //!< Gev^-1 -> mpc conversion factor
  double Gev2cpm(int n=1) const;
  double cpm2Gev(int n=1) const;

  double sInv2Gev(double s) const { return cpm2Gev()*sInv2cpm(s); }
  double Gev2sInv(double g) const { return g*cpm2sInv()*Gev2cpm(); } //!< Gev -> s^-1

  double tau(double t); //!< conformal time at given time t  
  double tau_equ() const {return tau__equ;} //!< conformal time at equality
  double tau_ls() const {return tau__ls;}  //!< conformal time at last scattering
  double tau_0()  const { return tau__0;} //!< conformal time today
  double tau_min() const {  return Tau2A->start(); } //!< minimal conformal time simulated
  double tau_max() const { return Tau2A->stop();  }   //!< maximum conformal time simulated
  double tau2t(const double); //!< return usual time t given conformal time tau
  //double tau2tRough(const double); //!< 
  double tau2z(const double);  //!< redshift as a function of conformal time
  double tau2a(const double); //!< this gives a(tau);
  double tau2adot(const double); //!< this gives da(tau)/dtau;

  double tau2adotoa(const double) ; //!< return \dot a / a

  double tau2H(const double); //<! return the conformal Hubble parameter da/dtau / a 
  double tau2Hubble(const double); //!< return the usual Hubble parameter da/dt / a = da/dtau / a^2

  double tau2rho_g(const double); //!< this gives rho_photons(tau);
  double tau2rho_nu(const double); //!< this gives rho_nu(tau) (massless);
  double tau2rho_nuNR(const double); //!< this gives rho_nuNR(tau);
  double tau2rho_b(const double); //!< this gives rho_b(tau);
  double tau2rho_cdm(const double); //!< this gives rho_cdm(tau);
  double tau2rho(const double); //!< total rho
  double tau2rho_v() const { return rho__v; } //!< return vacuum energy density
  double tau2rho_relativistic( const double); //!< return relativistic rho (rough in case of massive nu)
  double tau2rhodot(const double); //!< d/dtau of total rho
  double tau2rhodot_g(const double); //!< d/dtau of  rho gamma
  double tau2rhodot_cdm(const double); //!< d/dtau of  rho cold dark matter
  double tau2rhodot_b(const double); //!< d/dtau of  rho baryon
  double tau2rhodot_m(const double); //!< d/dtau of  rho matter
  double tau2rhodot_nu(const double); //!< d/dtau of  rho relat. neutrinos
  double tau2rhodot_nuNR(const double); //!< d/dtau of  rho relat. neutrinos

  double tau2grho(const double); //!< d/dtau of total rho
  double tau2grhodot(const double); //!< d/dtau of total rho

  double tau2R(const double); //!< 4.0/3.0 *rho_r / rho_b
  double tau2rho_m(const double); //!< this gives rho_m(tau);
  double tau2rs(const double); //!< this gives sound horizon(tau);
  double tau2cs2(const double); //!< this gives soundspeed^2(tau);

  //! Return physical size of horizon at time tau 
  double tau2horizon(const double tau) { return tau2a(tau)*tau; }
  //! Return physical size of horizon seen TODAY of horizon at tau (this is what you usually mean by 'horizon'
  double tau2horizon_seenToday(const double tau) const { return tau;}
  //! Return comoving k of horizon size, i.e.  2*pi/ (2*d_horizon()) = k
  double tau2k_horizon(const double tau) const { return 2.0*M_PI/tau;};
  //! Return time tau at which mode k will cross inside horizon
  double k_horizon2tau(const double k) const { return 2.0*M_PI/k; }

  double tau2p_nu(const double tau) { return tau2rho_nu(tau)/3.0; }  //!< pressure of massless nu's
  double tau2p_nuNR(const double); //!< pressure of massive nu's
  double tau2p_g(const double tau) { return tau2rho_g(tau)/3.0; } //!< photon pressure
  double tau2p_v() { return -tau2rho_v(); } //!< vacuum pressure
  double tau2p(const double);  // !< total pressure
  double tau2w_nuNR(const double tau); //!< equation of state for massive neutrinos
  double tau2wDot_nuNR(const double tau); //!< tau-derivative of the eqs of massive neutrinos

  virtual double tau2AvOmega_q(const double tau) { return 0;}

  double tau2w_nu(const double tau) {
    throw Bad_Error("Cosmos::tau2w_nu - unimplemented");
    if (nuNR() == 0) return 1.0/3.0; //!< for definitness, 0/0 = 1.0/3.0 :-)
    return tau2rho_nuNR(tau)/tau2p_nuNR(tau);
  }

  double a_0() const { return 1.0; } //!< scale factor today
  double a_equ() const { return a__equ; } //!< scale factor at equality
  double a_ls() { return z2a( z_ls() ); } //!< scale factor at last scattering
  double a_min() const { return splineA2Tau->start();}  //!< minimal a
  double a_max() const { return splineA2Tau->stop();} //!< maximum a
  double a(double t);  //!< scale factor as a function of time
  double a2t(double a); //!< time as a function of a
  double a2tau(double);  //!< conformal time as a function of a

  void tauBounds(const double, string="unknown") const;   //!< checks whether tau is within simulated bounds.
  void aBounds(double,string = "unknown");
  void tBounds(double,string = "unknown");

  double t_equ() const { return t__equ; }  //!< usual time at equality
  double t_ls() const { return t__ls; }  //!< time at last scattering 
  double t_0() const {return t__0; } //!< time today, i.e. age of the universe
  double t_min() const { return splineT2Tau->start(); } 
  double t_max() const { return splineT2Tau->stop(); }  
  double t2z(double); //!< redshift as a function of time 
  double t2a(double); //!< scale factor as a function of time
  double t2tau(double);  //!< conformal time as a function of time 
 
  double z_equ(); //!< Matter-Radiation equaltiy redshift
  //! Redshift of last scattering. If thermal history has been calculated it uses exact value, otherwise estimate
  double z_ls() { 
    if (validThermo()) return z__ls;
    return 1000*pow(omega_b(),-0.027/(1+0.11*log(omega_b())));}  //!< last scattering redshift
  double z2tau(double); //!< conformal time as a function of redshift
  double z2t(double); //!< time as a function of redshift

  // the visibilities etc, also cached splines
  double visibility(double);   //!< visibility function, see for instace Seljak & Zaldarriaga
  double dvisibility(double);  //!< derivative of visibility w.r.t. conformal time 
  double ddvisibility(double); //!< second deriv. w.r.t conformal time 
  double opac(double);  
  double d2kappa(double);  
  double expmmu(double); //!< exp(-kappa(tau)-kappa(tau_0) )

  double physical2comovingK(const double k, const double tau); //!< given a physical wavenumber k, return comoving wavenumber 


  double rho2omega(const double rho, double tau = -1); //!<   Omega_0 for  energy density rho 
  double omega2rho(const double omega);//!< return rho_0 given some Omega_0

  bool ValidOmega; //!< If this is true then history() has been calculated
  bool validOmega()  const { return ValidOmega; }  

	void setPkBaryon(bool p){pkBaryon=p;}
	bool pkWithBaryons(){return pkBaryon;}

  double omega_b(bool fromHistory=true); //!< omega Baryon
  double omega_cdm(bool=true);  //!< omega cold dark matter  if bool is true, try to take it from actual history (if already calculated)
  double omega_c() { return omega_cdm(); }  //!< Synonym for omega_cdm() Cold dark Matter
  double omega_v(); //!< omega vacuum excluding quintessence

  double omega_g(); //!< Gammas
  double omega_nu(); //!< relativistic neutrinos
  double omega_nuNR(bool fromHistory=true);
  double omega_relativistic() { return omega_g() + omega_nu(); }  // Total of massless species, disregards massive nu's

  double omega_m(bool fromHistory) { return omega_cdm(fromHistory) + omega_b(fromHistory); }   //!< Total matter, cold + baryonic
  double omega_m() { return omega_cdm() + omega_b(); }   //!< Total matter, cold + baryonic
  virtual double omega_0() { return omega_m() + omega_g() + omega_nu() + omega_nuNR() +  omega_v(); } //!< Return total omega today
  double rho_crit() const { return Rho_crit; } //!< today's critical energy density, only possible if validHistory() 

  double omegaH2_b() { return omega_b()*h2(); } //!< Omega_b * h^2  
  double omegaH2_cdm() { return omega_cdm()*h2();}  //!< Omega_cdm * h^2
  void setOmega_b(double x) { Omega_b = x; }  //!< set  Omega_baryon
  void setOmega_cdm(double x) { Omega_cdm = x;} //!< set Omega_cdm
  //! set omega_cdm such that a flat universe comes out
  virtual void setOmega_cdm_flat() {
     double omegaCdm = 1.-omega_v()-omega_b(false)-omega_g();
     if (nuR()>0)
        omegaCdm-=omega_nu();
     if (nuNR()>0)
       omegaCdm-=omega_nuNR(false);
     if (omegaCdm<0.)
       throw Bad_Error("Cosmos::setOmega_cdm_flat() - would need negative Omega_cdm");
     setOmega_cdm(omegaCdm);
  }
  virtual void setOmega_vacuum_flat() {
     double omegaVacuum = 1.-omega_cdm()-omega_b(false)-omega_g();
     if (nuR()>0)
        omegaVacuum-=omega_nu();
     if (nuNR()>0)
       omegaVacuum-=omega_nuNR(false);
     if (omegaVacuum<0.)
       throw Bad_Error("Cosmos::setOmega_vacuum_flat() - would need negative Omega_cdm");
     setOmega_vacuum(omegaVacuum);
  }

  enum AdjustNeutrinoMass { DoAdjustNeutrinoMass, DoNotAdjustNeutrinoMass }; // when setting Omega for massive nu
  void setOmega_nuNR(double x, AdjustNeutrinoMass adjustMass=DoAdjustNeutrinoMass) {
    Omega_nuNR = x;
    if(adjustMass == DoAdjustNeutrinoMass) {
      setNeutrinoMassFromOmegaNuNR();
    }
  }

  void setOmega_nu(double x) { Omega_nu = x;} //!< set Omega_nu
  void setOmega_vacuum(double x) { Omega_vacuum = x;} 

  void setOmega_g(double x) { Omega_g = x; }   //!< Photon density

  void setOmegaH2_b(double x) { Omega_b = x/h2(); }  //!< Omega * h^2 of baryons
  void setOmegaH2_cdm(double x) { Omega_cdm = x/h2();} //!< Omega * h^2 of cdm
  void setOmegaH2_vacuum(double x) { Omega_vacuum = x/h2();} //!< Omega * h^2 of vacuum  
  void setOmegaH2_nu(double x) { Omega_nu = x/h2();} //!< Omega * h^2 of neutrinos
  void setOmegaH2_nuNR(double x, AdjustNeutrinoMass a=DoAdjustNeutrinoMass) { setOmega_nuNR(x/h2(), a); } //!< Omega * h^2 of massive neutrinos
  void setOmegaH2_g(double x) { Omega_g = x/h2(); }   //!< Omega *h^2 of photons

  double rstar() { return tau2rho_relativistic(tau_ls()) / tau2rho_m(tau_ls()); } //!< relativistic / matter at ls 
  double l_A() { return M_PI/tau2rs(tau_ls()) * (tau_0() - tau_ls()); }; //!< return the acoustic scale 

  double Z2iH(double z) { return Z2iH_spline->fastY(z); }  //!< Translates from redshift to inverse Hubble parameter

  //Distances
  double propermotionDistance(double z) const; //!< return proper motion distance
  double propermotionDistance(double z0, double z1) const; //!< return proper motion distance
  double angulardiameterDistance(double z) const; //!< return angular diameter distance
  double luminosityDistance(double z) const; //!< return luminosity distance, ignoring difference between heliocentric and CMB frame redshifts
  double luminosityDistance(double zhel, double zcmb) const; //!< return luminosity distance

  //Distance derivatives.  In a flat cosmology, the derivate of the proper
  // motion distance is just Z2iH, the inverse Hubble parameter
  double propermotionDistanceDeriv(double z) const { return Z2iH_spline->fastY(z); } //!< return derivate of proper motion distance with respect to redshift
  double angulardiameterDistanceDeriv(double z) const { 
    double opz=1+z; return Z2iH_spline->fastY(z)/opz - propermotionDistance(z)/(opz*opz); } //!< return derivate of angular diameter distance with respect to redshift
  double luminosityDistanceDeriv(double z) const { return propermotionDistance(z) + (1+z)*Z2iH_spline->fastY(z); } //!< return derivative of luminosity distance with respect to redshift


  inline double nuR() { return NrRelativisticNu;} //!< nr of relativisic nu species
  inline double nuNR() const {return mNrNonRelativisticNu;} //!< nr of non-relativstic nu species

  //! set nr of relativistic nu's
  void setNuR(double x) {
    NrRelativisticNu = x;
    if (NrRelativisticNu==0)
      setOmega_nu(0.);
  }

  //! set nr of non-relativistic nu's
  void setNuNR(double x) { mNrNonRelativisticNu = x; }

  virtual double mass_nu_eV() const {return mMassNu_eV;} //!< return mass of massive neutrino in eV
  virtual double mass_nu_cpm() const {return mMassNu_eV*1e-9*Gev2cpm();} //!< return mass of massive neutrino in mpc^-1
  virtual double mass_nu_kbT() const {return Gev2cpm()*mMassNu_eV*1e-9/kbT(T_nu());} //!< return mass of massive neutrino in units of k_b*T
  virtual void setMass_nu_eV(const double x) {
    mMassNu_eV = x;
    setOmegaNuNRFromNeutrinoMass();
  }

  void setOmegaNuNRFromNeutrinoMass();
  void setNeutrinoMassFromOmegaNuNR();

  virtual bool massiveNeutrinosAreRelativistic(const double); //!< true, if mass_nu() < k_B * T
  //! true if we include massive neutrinos
  bool haveMassiveNeutrinos() const { return /*const_cast<Cosmos*>(this)->omega_nuNR() > 0 &&*/ nuNR() > 0; }

  double massiveNeutrinoPressure(int nr, double a);
  double massiveNeutrinoRho(int nr, double a);

  double h() const { return h_h; } //!< Small Hubble h today
  //! Setting h also entails setting a different Omega_gamma and Omega_nu
  void seth(double h) {
    h_h = h;
    setOmegaH2_nu(1.695e-5);   // peakock, P 664
    setOmegaH2_g(2.488e-5);
    setOmegaNuNRFromNeutrinoMass();
  }
  double h2() const { return h_h*h_h;} //!< h^2
  double H_0() const { return 3.242e-18*h(); } //!< Hubbles constant in s^-1
  double H_0_Gev() { return sInv2Gev(H_0()); }
  double H_0_cpm() { return sInv2cpm(H_0()); } //!< Hubble constant in inverse Mpc

   double Y_p() { return Y_P;}  //!< Helium abundance
   double Y_he() { return Y_P;} //!< Synonym for Y_p()
   void setY_he(double x) { Y_P = x; }

  double rho_g0() const { return M_p()*M_p()* pow(T_cmb(),4) * 1.4952e-13;}   //!< Todays Photon energy density from temperature of microwave Background
  double rho_nu0() const { return rho_g0()*7./8.*pow(4./11., 4./3.); }
  double rho_nuNR0();
  double rho_0() const { return  M_p()*M_p()*h2()* 3.3379e-7;  } ; // today's critical energy density
  // double rho_m0() {return M_p()*M_p()*h2()* 3.3379e-7;  } ;  //!< Today's matter energy density divided by omega_m(), i.e you have to multiply with omega_m() to really get the energy density. Sorry for this slight misnomer, but most of the time, this and not omega_m()* ... is needed.

  double rho_b0() { return rho_0() * omega_b(); } //!< today's energy baryon density in Mpc^-4
  double rho_cdm0() { return rho_0() * omega_cdm(); } //!< today's cold dark matter energy density in Mpc^-4  

  double initialRho_g(const double a) { return  rho_g0()*pow(a/a_0(),-4); } //!< Initital photon density for a given a. For a = a0, this is of course rho_g0 
  //! Initital massless neutrino density for a given a. This is the sum for _all_ relativistic species of neutrinos
  double initialRho_nu(const double a) { return  nuR()*rho_nu0()*pow(a/a_0(),-4); }
  double initialRho_b(const double a) { return   rho_b0()*pow(a/a_0(),-3); }
  double initialRho_cdm(const double a) { return    rho_cdm0()*pow(a/a_0(),-3); }
  virtual double initialRho_nuNR(const double a, double* p_nuNR=0);

  virtual void setExtraBackgroundTimesteps(const std::vector<double>& extraSteps);  //!< list of tau-values for extra timesteps
  /*! call this in reimplementations of history() to make sure that tau+dtau
   *  will not jump over the next entry in mExtraTimeSteps.
   */
  void checkForExtraTimestep(double tau, double& dtau);
  virtual void propagateHistoryInTau(const double , const double*,double *); //!< Return array ziel[] with the tau derivatives of a, rho_r etc for the background equations
  //! Wrapps virtual propageHistoryInTau for use with odeint()
  void  propagateHistoryInTauWrapper(const double tau, const double* y,double *ydot) { propagateHistoryInTau(tau,y,ydot);}

  //! fills some history splines.
  virtual Background fillHistorySplines(double tau, const double *y, const double *ydot, bool writeSplines=true);

  virtual void history(bool inform=false);
  bool validHistory() const { return ValidHistory; }  //!< history() sets this flag  
  bool validThermo() const { return ValidThermo; }  //!< thermo() sets this flag  
  virtual void findSpecialMoments(); //!< get tau__ls, etc. to hight accuracy
  double findEquality(const double); //!< needed by findSpecialMoments()


  /*! CDM - Powerspectrum SplineWeb 
    
    On calling reset(), the initSplines() function will
    built up SplineWebs for the power spectra of cdm, baryons,
    photons, neutrinos (and the things you may be interested in)

    The SplineWebs are anchored in SplineAnchor
    and hence a new call to reset() or the destruction
    of cosmos will delete these SplineWebs. 
    
    In short: Do not use the old pointers
    after killing or reset-ing  cosmos.
  */
  SplineWeb* power_cdm() const { return Power_CDM;} //<! Return a pointer to the cdm-power web 
  SplineWeb* power_baryon() const { return Power_b;} 
  SplineWeb* power_gamma() const { return Power_g;}
  SplineWeb* power_nu() const { return Power_nu; }
  SplineWeb* power_nuNR() const { return Power_nuNr; }
  SplineWeb* power_psi() const { return Power_psi; }

  void setPower_cdm(string s) { Power_CDM = new SplineWeb(s, &SplineAnchor,50,500); } // size is good guess.
  void setPower_baryon(string s) { Power_b = new SplineWeb(s, &SplineAnchor,50,500); }
  void setPower_gamma(string s) { Power_g= new SplineWeb(s, &SplineAnchor,50,500); }
  void setPower_nu(string s) { Power_nu = new SplineWeb(s, &SplineAnchor,50,500); }
  void setPower_nuNR(string s) { Power_nuNr = new SplineWeb(s, &SplineAnchor,50,500); }
  void setPower_psi(string s) { Power_psi = new SplineWeb(s, &SplineAnchor,50,500); }

  virtual std::vector<SplineWeb*> powerSplineWebs() const {
    std::vector<SplineWeb*> vec;
    vec.push_back(power_cdm()); vec.push_back(power_baryon());
    vec.push_back(power_gamma()); vec.push_back(power_nu());
    vec.push_back(power_nuNR()); vec.push_back(power_psi());
    return vec;
  }

  MassFunction* dumpMassFunction(string name, double z, bool print, double norm); //!< Write mass function to disk in file name

  Spline* dumpNumberDensity(string name, double m); //!< Write number density
 // Returns comoving volume element dV
  double comovingVolume(double z0, double z1);
	
  double integrateOverComovingVolume(string name, double m, double z0, double z1, double t1, double t2, double p1, double p2);

// Integrates over a whole redshift range
  double integrateComovingVolume(double z0, double z1);

  Spline* allMatterPowerSpline(double tau); 

  void dumpPower(int n,string name, SplineWeb* w =0, double z=0); //!< Write power spectrum to disk in file name
  Spline* createPower(int n,string name, SplineWeb* w =0, Anchor *a=0, double tau=0); //!< return Spline as a section of power Spline Web
  void fillPower(SplineWeb* power,double k, double tau, double delta); //!< set new point at (k,tau) of value delta for power web *power. This is a low-level routine
  void fillPower(Perturbation *pert, double tau); //!< Asks perturbation to return all delta's needed to fill the power-webs available. 
  void fillPsiWeb(SplineWeb* power,double k, double tau, double delta);//!< like fillPower, but stores "delta" without k factors as function of k, not k/h

  enum GrowthFactorType {Baryons, Cdm, Weighted};

  Spline* growthFactor(const double k=0.05, std::string name="", GrowthFactorType species = Cdm, double* norm = 0, Anchor* a=0);
  
  Spline* growthIndex(Spline* growth_factor=0, Anchor *a=0, string name="");
  Spline* growth2Z;
  Spline* growthDt2Z;
  Spline* growthDt2t;
  Spline* growthT2t;
  Spline* growthOverA2Z;

  double ageYr;
  Spline *a2year;
  Spline *year2a;

  void dumpGrowthFactor(double);
  void printOutputList(char*, double, double, double);

  double Sigma8_z(double z);

  void setInitialPower(double min=1.0, double max=1.0, int steps=1);
  double Pk2Delta2(Spline *P, double k);
  //  Spline* createDelta2(string name, SplineWeb* power, double tau=0, Anchor* a=0); //!< create delta^2 of given species
  Spline * createLensingPowerLimber( int n, string name, int lstart, int lend, Anchor* a=0 ); //!< returns the Cl-spectrum of the gravitational potential times l^4, used by AllSkyLensing

  void setRecombinationParameter(vector<double>& v) { RecombinationParameter=v;}
  Recombination* getRecombination() { return recombination;}

  double T_cmb() const { return T_cmB;}  //!< Temperature of Microwavebackground
  void setT_cmb(double x) {T_cmB = x; }
  double T_nu() const { return T_cmb() * pow(4./11., 1./3.); }

  double twooverh3() const;

  double M_p() const { return baseCosmos::M_p(); }
  double M_p(int n) const;

  static double Gpi8() { const double Mp=baseCosmos::M_p(); return 1.0/(Mp*Mp); } //!< 8*pi*G in Mpc ^2
  double kbT(const double T=1) const {
    static const double k = 8.617343e-5; // Boltzmann constant eV K^−1
    return k*1e-9*Gev2cpm()*T;
  } //!< k_boltzmann times Temperature (in K) return value in Mpc^-1

  double omega0H2() { return omega_0() * h2(); }
  double omegabH2() { return omega_b() * h2(); }
  double k_equ() { 

    return 7.46e-2*pow(T_cmb()/2.7,-2)*omega0H2();

    return sInv2cpm(H_0())* sqrt(2*omega_0()*a_0()); }

  virtual void printStatus(ostream& o=cout);
  virtual void printStatus(const char*);

  void setReionizationFraction(double f) { ReionizationFraction = f;}
  double reionizationFraction() { return ReionizationFraction;}
  void setReionizationZ(double z) { ReionizationZ = z;}
  double reionizationZ() { return ReionizationZ;}
  void setReionizationZStop(double z) { ReionizationZStop = z;}
  double reionizationZStop() { return ReionizationZStop; }


  double x_e(const double);
  
  /*! set optical distance to last scattering surface. If 
    this is set to zero, it also sets the reionization fraction 
    to zero if it is non-zero, it will set the reionization fraction to unity*/
  void setOptDistanceLss(double o) { Optdlss = o; if (o) setReionizationFraction(1); else setReionizationFraction(0); } 
  double optDistanceLss() { return Optdlss;}
  
  
  double reionizationTau() { return z2tau(reionizationZ()); }
  double reionizationTauStop() { return z2tau(reionizationZStop()); }



  virtual void getReady(); //!< set up the mass of the neutrino amnu
  //  void onlyThermo(); //!< convenience wrapper around finithermo()


  // ======== CMBFAST 

  double dtauda(double a);


  void ionize(double *tempb, double *a, double * adot, double *dtau, double *xe);
  void ionhe(double *tempb, double *a, double *x0, double *x1, double *x2);
  void recint(double *a, double *xe);
  void reiopar();
  
  void cobeNormalize(const ControlPanel&,CL&, const vector<int> &lval);


  //! Return initial power-law power spectrum propto k^n
  static double powerSpectrum(const double k, const double n,const double dnsdlnk=0.0);
  double scalarSpectrum(const double k, unsigned int q); //!< return scalar power spectrum with spectral index Intialower[q]
  double tensorSpectrum(const double k, int n,const ControlPanel&); //!< return the initial tensor spectrum for initial value number n, i.e. InitialTensorPower[n]
  // ------------ thermal history and timestep determination ---------
  // void finithermo(const double, const double, const double, double *taurend, double *dlntau0, int *n1, int *nstep, bool);

 void thermo(double MaxK, ControlPanel&);
 void integrateDotKappa(const double z, double *kappa, double *dotkappa);
 void determineTimeSteps();
 double determineStartOfThermo(double k); //!< Determine minimal tau for thermal history
 virtual double tauStartforKMode(double k); //!< at which time to start computing perturbations for wavelength k

};



/*!
  The Caching mechanism used here is quite simple.
  For some spline, say Tau2A, there are two double
  variables created: Tau2AX and Tau2AY. 
  Now, If Tau2AX equals the argument of
  tau2a(tau) then tau2a(tau) will return Tau2AY.
  If the argument is not the same, it will call
  the spline interpolation, set the Tau2AY value
  to this and also set the Tau2AX to the argument
  tau. 
  
  This increases the speed of functions
  such as fderivs() by about 20%. This and only this is the
  reason why we don't just allocate all possible
  variables in each function, just to avoid calling
  tau2a() many times in one function call: Its cached,
  its extremly fast and its clean, cause we don't need
  any more variables in each function.

  Drawback: If the Spline changes (say by disarm(),
  modification and arm(), or just because you
  called reset() on cosmos, you will get the
  *old* and invalid interpolation value, *if* you
  call e.g. tau2a(tau) where tau is the last
  argument you used for tau2a() *before* you
  modified the Tau2A spline (which as mentioned
  happens on reset()ing the cosmos.)

  However, reset() will set the X-value to some
  insane number and you can be sure to get an
  uncached value just after calling reset().

*/
//#define CACHE(myspline)  if (myspline##X == tau) return myspline##Y; myspline##X = tau; return myspline##Y = (*myspline)(tau)
#define CACHE(myspline)  return (*myspline)(tau)



//                           tau -> 
inline double Cosmos::tau2a(const double tau) { CACHE(Tau2A); }
inline double Cosmos::tau2adot(const double tau)  { CACHE(Tau2ADot); }
inline double Cosmos::tau2adotoa(const double tau)  { CACHE(Tau2ADotToA); }

inline double Cosmos::tau2z(const double tau)  { return a2z(tau2a(tau)); }

inline double Cosmos::tau2rho_g(const double tau)  { CACHE(Tau2Rho_g); }
inline double Cosmos::tau2rho_nu(const double tau)  { CACHE(Tau2Rho_nu); }
inline double Cosmos::tau2rho_nuNR(const double tau)  { CACHE(Tau2Rho_nuNR); }
inline double Cosmos::tau2rho_b(const double tau)  { CACHE(Tau2Rho_b); }
inline double Cosmos::tau2rho_cdm(const double tau)  { CACHE(Tau2Rho_cdm);}
inline double Cosmos::tau2rho_m(const double tau)  { CACHE(Tau2Rho_m); }
inline double Cosmos::tau2rho(const double tau)  { CACHE(Tau2Rho); }


inline double Cosmos::tau2rhodot(const double tau)  { CACHE(Tau2RhoDot); }
inline double Cosmos::tau2rhodot_b(const double tau)  { CACHE(Tau2RhoDot_b); }
inline double Cosmos::tau2rhodot_cdm(const double tau)  { CACHE(Tau2RhoDot_cdm); }
inline double Cosmos::tau2rhodot_m(const double tau)  { CACHE(Tau2RhoDot_m); }
inline double Cosmos::tau2rhodot_g(const double tau)  { CACHE(Tau2RhoDot_g); }
inline double Cosmos::tau2rhodot_nu(const double tau)  { CACHE(Tau2RhoDot_nu); }
inline double Cosmos::tau2rhodot_nuNR(const double tau)  { CACHE(Tau2RhoDot_nuNR); }



inline double Cosmos::tau2grho(const double tau)  { CACHE(Tau2GRho); }

inline double Cosmos::tau2grhodot(const double tau)  { CACHE(Tau2GRhoDot); }



inline double Cosmos::tau2p_nuNR(const double tau)  { CACHE(Tau2P_nuNR); }
inline double Cosmos::tau2p(const double tau)  { CACHE(Tau2P); }

inline double Cosmos::tau2w_nuNR(const double tau)  { CACHE(Tau2w_nuNR); }
inline double Cosmos::tau2wDot_nuNR(const double tau)  { CACHE(Tau2wdot_nuNR); }

inline double Cosmos::tau2R(const double tau) { CACHE(Tau2R);}
inline double Cosmos::tau2t(const double tau) { CACHE(Tau2T);}
inline double Cosmos::tau2rs(const double tau) {  CACHE(Tau2Rs);}
inline double Cosmos::tau2cs2(const double tau) {  CACHE(Tau2Cs2);}


//                           a  ->
inline double Cosmos::a2tau(double a) {  return (*splineA2Tau)(a);}
inline double Cosmos::a2t(double a) {  return tau2t(a2tau(a));}

//                           t -> 
inline double Cosmos::t2tau(double t)  {return (*splineT2Tau)(t);}
inline double Cosmos::t2a(double t) {return tau2a(t2tau(t)); }
inline double Cosmos::t2z(double t) { return tau2z(t2tau(t)); }

//                           z -> 
inline double Cosmos::z2tau(double z) {  return a2tau(z2a(z)); }
inline double Cosmos::z2t(double z) { return tau2t(z2tau(z));}

//   visibilities etc
inline double Cosmos::visibility(const double tau) { CACHE(Visibility);}
inline double Cosmos::dvisibility(const double tau) { CACHE(DVisibility);}
inline double Cosmos::ddvisibility(const double tau) { CACHE(DDVisibility);}
inline double Cosmos::opac(const double tau) { CACHE(Opac);}
inline double Cosmos::d2kappa(const double tau) { CACHE(D2Kappa);}
inline double Cosmos::expmmu(const double tau) { CACHE(Expmmu);}

#undef CACHE


/*! 
  The first of various bound checkers. Checks validHistory() and
  wheter the requested x value (here, it is tau) is within the tabulated range 
*/
inline  void Cosmos::tauBounds(double tau, string s) const {
  if (validHistory()) {
    if (tau >= tau_min() && tau <= tau_max()) return;
    throw Bad_Error("Cosmos::" + s + "() tau is out of computed range");
  }
  throw Bad_Error("Cosmos::" + s + "() no Valid history. Please calculated one first");
}

inline void Cosmos::aBounds(double a, string s) {
  if (validHistory()) {
    if (a >= a_min() && a <= a_max() ) return;
    throw Bad_Error("Cosmos::" + s + "() a is out of computed range");
  }
  throw Bad_Error("Cosmos::" + s + "() no Valid history. Please calculated one first");
}

inline void Cosmos::tBounds(double t, string s) {
  if (validHistory()) {
    if (t >= t_min() && t <= t_max() ) return;
    throw Bad_Error("Cosmos::" + s + "() a is out of computed range");
  }
  throw Bad_Error("Cosmos::" + s + "() no Valid history. Please calculated one first");
}

#endif

