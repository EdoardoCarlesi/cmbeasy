#ifndef cncosmos_H
#define cncosmos_H 

#include "quintcosmos.h"
//#include <math>
#include "global.h"

/*!
  CnCosmos implements the background for mass-varying neutrino scenarios,
  i.e. neutrinos coupled to dark energy.
*/
class CnCosmos : public QuintCosmos
{
  public: 
   CnCosmos(Quintessence::Type = Quintessence::none); ~CnCosmos();
  
    virtual void propagateHistoryInTau(const double , const double*,double *);
    //!< here is where the background evolution equations are defined void
    void history(bool inform=false); //!< re-implemented from Cosmos
  
    void printStatus(ostream &o=cout);
  
    // reimplementations of methods relating to massive neutrino handling
    virtual void getReady(); //!< set up the mass of the neutrino amnu virtual
    bool massiveNeutrinosAreRelativistic(const double tau); Background
	    fillHistorySplines(double tau, const double *y, const double *ydot,
			    bool writeSplines=true);
  
    //double mass_nu_eV() const { throw Bad_Error("CnCosmos::mass_nu_eV() -
    //call the phi-dependent routine instead."); } //!< coupled neutrinos - use
    //mass_nu_xx(phi) instead
    double mass_nu_eV_uncoupled() const { return mMassNu_eV; }
    //double mass_nu_cpm() const { throw Bad_Error("CnCosmos::mass_nu_cpm() -
    //call the phi-dependent routine instead."); } //!< coupled neutrinos - use
    //mass_nu_xx(phi) instead double mass_nu_kbT() const { throw
    //Bad_Error("CnCosmos::mass_nu_kbT() - call the phi-dependent routine
    //instead."); } //!< coupled neutrinos - use mass_nu_xx(phi) instead
    double mass_nu_kbT_uncoupled() const { return
	    Gev2cpm()*mass_nu_eV_uncoupled()*1e-9/kbT(T_nu()); } 
    double mass_nu_eV(const double phi) const { //cout << "phi: " << phi << "
		    //Mp()  " << M_p() << " massnuev exp:  " <<
		//	    exp(mNeutrinoCoupling*phi/(M_p())) << endl; 
	    /* added a phi*phi dependence on neutrino mass!! */
return mMassNu_eV*exp(timeVaryingNuCoupling(phi)*phi/(M_p())); } //!< return
//mass of massive neutrino in eV                
	    //eturn mMassNu_eV*exp(mNeutrinoCoupling*phi/(M_p())); } //!<
	    //return mass of massive neutrino in eV
    double mass_nu_cpm(const double phi) const { return
	    mass_nu_eV(phi)*1e-9*Gev2cpm(); } //!< return mass of massive
	//    neutrino in mpc^-1 
	double mass_nu_kbT(const double phi) const {
		    return Gev2cpm()*mass_nu_eV(phi)*1e-9/kbT(T_nu()); } //!<
		    //return mass of massive neutrino in units of k_b*T 
		    
		    double dMass_nu_kbTdphi(const double phi) const { 

      //    cout << " coupling: " << mNeutrinoCoupling << " mass nu: " <<
      //    mass_nu_eV(phi) << " massNu Uncoupled " << mMassNu_eV  << " kbT "
      //    << kbT(T_nu())  <<endl;
return mNeutrinoCoupling*mass_nu_kbT(phi); } //!< d m_nu(phi)/dphi double

double avarmnu(const double phi) const { return mass_nu_kbT(phi); }
  
    void setNeutrinoCoupling(const double beta) { mNeutrinoCoupling = beta; }
    void setNeutrinoCouplingFromOmega_quint() {
//X       Exponential* exp = dynamic_cast<Exponential*>(quintessence()); X
//if (!exp)
//X         throw Bad_Error("CnCosmos::setNeutrinoCouplingFromOmega_quint() - no exponential quintessence set.");
//X 
//X       const double lambda = exp->lambda();
//X       const double omega_q(false/*from history*/))
//X       const double beta = 1.-
//X       mNeutrinoCoupling = beta;
    }
    double neutrinoCoupling() const { return mNeutrinoCoupling; }

    void dumpNuMass() {
	//TauNuMass->dump("nu_mass");
	Tau2Rho_nuNR->dump("nu_mass-rho_nuNR");	
	}

    double timeVaryingNuCoupling(double phi) const { 
	    
	    double tDependentCoupling;
	    tDependentCoupling = mNeutrinoCoupling*(exp(abs(phi/M_p())));	    
//	    return tDependentCoupling;  
    return mNeutrinoCoupling; // returns constant coupling
    }

    double dBetadphi(const double phi) const { return 0; }

    void setMass_nu_eV_uncoupled(const double m) {
      setMass_nu_eV(m);
    }

    void setMass_nu_kbT_uncoupled(const double m) {
      mMassNu_eV=m/(Gev2cpm()*1e-9/kbT(T_nu()));
      setOmegaNuNRFromNeutrinoMass();
    }

  private:
   double mNeutrinoCoupling;
	Spline* TauNuMass;
   };

#endif  // cncosmos_H

