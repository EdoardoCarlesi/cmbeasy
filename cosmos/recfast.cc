#include "recfast.h"

#include "cosmos.h"

#include <fstream>


void Recfast::compute()
{
  //  ofstream helium("helium_old.dat");
  C  = 2.99792458e8;
  k_B = 1.380658e-23;
  h_P = 6.6260755e-34;
  static const double m_e = 9.1093897e-31;
  //av. H atom; note: neglecting deuterium, making an O(e-5) effect
  m_H = 1.673725e-27;
  not4=3.9715; // mass He/H atom
  static const double sigma = 6.6524616e-29;
  static const double a = 7.565914e-16;
  static const double G = 6.6742e-11;    // new value
  // Fundamental constants in SI units
  // ("not4" pointed out by Gary Steigman)

  static const double L_H_ion=1.096787737e7;     // level for H ion. (in m^-1)
  static const double L_H_alpha=8.225916453e6;   // averaged over 2 levels
  static const double L_He1_ion=1.98310772e7;    // from Drake (1993)
  static const double L_He2_ion=4.389088863e7;   // from JPhysChemRefData (1987)
  static const double L_He_2s=1.66277434e7;      // from Drake (1993)
  /* class member */  L_He_2p=1.71134891e7;      // from Drake (1993)
  // photon rates and atomic levels in SI units

  double zinitial = 1e4;
  double z = zinitial;
  double zfinal = 0;

  // this will be the redshift at which ion()
  //  switches from tight for T_b to non-tight
  SwitchTbEvolution = 0;

  OmegaB = cosmos.omega_b();
  OmegaT = cosmos.omega_m();       // total matter
  Tnow = cosmos.T_cmb();
  double Yp = cosmos.Y_he();

  double mu_H = 1.0/(1.0-Yp);       // Mass per H atom
  mu_T = not4/(not4-(not4-1.)*Yp);  // !Mass per atom
  fHe = Yp/(not4*(1.0-Yp));        // n_He_tot / n_H_tot

  H = cosmos.H_0();                 // in seconds (that's originals recfast's HO)
  //double bigH = 100e3/(1e6*3.0856775807e16);
  //H = cosmos.h()*bigH;

  unsigned int dHdzSize = cosmos.splineA2Tau->size();
  mdHdz = new Spline(dHdzSize, "Recfast::mdHdz", &convenient);
  Spline hzSpline(dHdzSize, "Recfast::compute::hzSpline");
  for (unsigned int i=0; i<dHdzSize; ++i) {
    double z = cosmos.a2z(cosmos.splineA2Tau->x(i));
    double hubble = cosmos.tau2Hubble(cosmos.splineA2Tau->y(i));
    hzSpline.set(z, hubble);
  }
  hzSpline.flip();
  hzSpline.arm();
  hzSpline.derive(*mdHdz);
  mdHdz->arm();

  //  today's number density of hydrogen
  Nnow = 3.*H*H*OmegaB/(8.0*M_PI*G*mu_H*m_H);
  // UNUSED double n = Nnow * pow(1 +z,3);
  // UNUSED double fnu = (cosmos.nuR()*7.0/8.0)*pow(4.0/11.0, 4.0/3.0);

  // set up some constants so they don't have to be calculated later

  const double Lalpha = 1.0 / L_H_alpha;
  const double Lalpha_He = 1.0 / L_He_2p;
  const double DeltaB = h_P*C*(L_H_ion-L_H_alpha);
  CDB = DeltaB/k_B;
  double DeltaB_He = h_P*C*(L_He1_ion-L_He_2s); // 2s, not 2p
  CDB_He = DeltaB_He/k_B;
  CB1 = h_P*C*L_H_ion/k_B;
  CB1_He1 = h_P*C*L_He1_ion/k_B;                // ionization for HeI
  CB1_He2 = h_P*C*L_He2_ion/k_B;                // ionization for HeII
  CR = 2.0*M_PI*(m_e/h_P)*(k_B/h_P);
  CK = pow(Lalpha,3)/(8.0*M_PI);
  CK_He = pow(Lalpha_He,3)/(8.0*M_PI);
  CL = C*h_P/(k_B*Lalpha);
  CL_He = C*h_P/(k_B/L_He_2s);                  // comes from det.bal. of 2s-1s
  CT = (8.0/3.0)*(sigma/(m_e*C))*a;
  Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B;

  double y[4]; // 4 dim for 3 variables for odeint(); 

  //    Matter departs from radiation when t(Th) > H_frac * t(H)
  //    choose some safely small number
  H_frac = 1e-3;

  //    Fudge factor to approximate for low z out of equilibrium effect
  unsigned int Hswitch = mControl.recfastHSwitch;
  if (Hswitch==0) {
    fu=1.14;
  } else {
    fu=1.105;
  }

  b_He = 0.86;


  //   Set initial matter temperature
  y[3] = Tnow*(1.0+z);      //Initial rad. & mat. temperature
  //Tmat = y[3];

  getInit(z); // set x_H0, x_He0 and x0
  y[1] = x_H0;
  y[2] = x_He0;

  // OK that's the initial conditions
  int Nz=10000;
  //int Nz=1000;

  double hnext = 1e-1; // initial stepsize guess for odeint
  for (int i = 1; i <= Nz;i++) {
    //     calculate the start and end redshift for the interval at each z
    //     or just at each z
    double zstart = zinitial + (i-1)*(zfinal-zinitial)/(double)Nz;
    double zend   = zinitial + i*(zfinal-zinitial)/(double)Nz;

    //     Use Saha to get x_e, using the equation for x_e for ionized helium
    //     and for neutral helium.
    //     Everything ionized above z=8000.  First ionization over by z=5000.
    //     Assume He all singly ionized down to z=3500, then use He Saha until
    //     He is 99% singly ionized, and *then* switch to joint H/He recombination.

    double dy[4];  // for getting the accurate d (ln T_baryon / dz) used only at late times
    double dTbTb = 1.0 / (1.0 + z); //   d (ln T_baryon) / dz
    double rhs;
    z = zend;
    if (zend > 8000) {
      x_H0 = 1;
      x_He0 = 1;
      x0 = 1.0 + 2.0*fHe;
      y[1] = x_H0;
      y[2] = x_He0;
      y[3] = Tnow*(1.0+z);
    } else if (z > 5000) {
      x_H0 = 1;
      x_He0 = 1;
      rhs = exp( 1.5 * log(CR*Tnow/(1.0+z)) - CB1_He2/(Tnow*(1.0+z)) ) / Nnow;
      rhs *= 1.0;       // ratio of g's is 1 for He++ <-> He+
      x0 = 0.5 * ( sqrt( pow( rhs-1.0 -fHe,2.0)  + 4.0*(1.0 +2.0*fHe)*rhs) - (rhs-1.0-fHe) );
      y[1] = x_H0;
      y[2] = x_He0;
      y[3] = Tnow*(1.0+z);
    } else if (z > 3500) {
      x_H0 = 1.0;
      x_He0 = 1.0;
      x0 = x_H0 + fHe*x_He0;
      y[1] = x_H0;
      y[2] = x_He0;
      y[3] = Tnow*(1.0 +z);
    } else if (y[2] > 0.99) {
      x_H0 = 1.0;
      rhs = exp(1.5 * log(CR*Tnow/(1.0+z)) - CB1_He1/(Tnow*(1.0 +z)) ) / Nnow;
      rhs *= 4.0;      //ratio of g's is 4 for He+ <-> He0
      x_He0 = 0.5 * (sqrt( pow(rhs-1.0,2.0) + 4.0*(1.0+fHe)*rhs ) - (rhs-1.0));
      x0 = x_He0;
      x_He0 = (x0 - 1.0)/fHe;
      y[1] = x_H0;
      y[2] = x_He0;
      y[3] = Tnow*(1.0+z);
      //      helium << z << " " << x_He0 << endl;
    } else if (y[1] > 0.99) {
      rhs = exp( 1.5 * log(CR*Tnow/(1.0+z)) - CB1/(Tnow*(1.0+z)) ) / Nnow;
      x_H0 = 0.5   * (sqrt( rhs*rhs+4.0*rhs ) - rhs );

      hnext = Miscmath::odeint(y, 3, zstart, zend, 1e-9, hnext, 0,
                               (moDerivs)&Recfast::ion, *this, false);
      y[1] = x_H0;
      x0 = y[1] + fHe*y[2];
      ion(zend,y,dy); // to get derivative of Tb
      dTbTb = dy[3] / y[3]; 
      //      helium << z << " " << y[2] << endl;

    } else {
      hnext = Miscmath::odeint(y, 3, zstart, zend, 1e-9, hnext, 0,
                               (moDerivs)&Recfast::ion, *this, false, Miscmath::bstoer);
      x0 = y[1] + fHe*y[2];
      ion(zend,y,dy); // to get derivative of Tb
      if (SwitchTbEvolution > 0) {  // so we have switched already
        if (SwitchTbEvolution - zend > 5) { // drop the first z's
          dTbTb = dy[3] / y[3]; 
        } // i.e. do nothing if it just switched, keep canonical dTbTb value
      }
    }
    noteX->set(-zend, x0);  // set the spline
    noteTb->set(y[3]); // y[3] holds baryon temperature
    noteTbDeriv->set( -(1+zend) * dTbTb); //  d ln Tb / d ln a
    // eqn (68) of Ma & Bertschinger, the + instead of minus comes from
    // d ln Tb_ dln a = -(1+z) d ln Tb / dz
    double cs2 = y[3] * BoltzmannOverWeightAndC2() * (1.0 +  1.0/3.0 *(1+zend)*dTbTb);
    noteCs2->set(cs2);
  }  // end the loop over the z's
  noteX->arm(Spline::all);
  //  noteX->dump("recombin");
}



void Recfast::getInit(double z)
{
  double rhs;
  if (z > 8000) {
    x_H0 = 1.0;
    x_He0 = 1.0;
    x0 = 1.0+2.0*fHe;
  } else if (z > 3500) {
    x_H0 = 1.0;
    x_He0 = 1.0;
    rhs = exp( 1.5 * log(CR*Tnow/(1.0+z)) - CB1_He2/(Tnow*(1.0+z)) ) / Nnow;
    rhs *= 1.0;		// ratio of g's is 1 for He++ <-> He+
    x0 = 0.5 * ( sqrt( pow(rhs-1.0-fHe,2) + 4.0*(1.0+2.0*fHe)*rhs) - (rhs-1.0-fHe) );
  } else if (z > 2000) {
    x_H0 = 1.0;
    rhs = exp( 1.5 * log(CR*Tnow/(1.0+z)) - CB1_He1/(Tnow*(1.0+z)) ) / Nnow;
    rhs *= 4.0;		//ratio of g's is 4 for He+ <-> He0
    x_He0 = 0.5*( sqrt( pow(rhs-1.0,2) + 4.0*(1.0+fHe)*rhs ) - (rhs-1.0));
    x0 = x_He0;
    x_He0 = (x0 - 1.0)/fHe;
  } else {
    rhs = exp( 1.5 * log(CR*Tnow/(1.0+z)) - CB1/(Tnow*(1.0+z)) ) / Nnow;
    x_H0 = 0.5 * (sqrt(rhs*rhs +4.0*rhs ) - rhs );
    x_He0 = 0;
    x0 = x_H0;
  }
}



void Recfast::ion(const double z, const double *y, double *f)
{
  double CfHe_t;
  static const double Lambda = 8.2245809;
  static const double Lambda_He = 51.3;

  //      the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen	
  static const double a_PPB = 4.309;
  static const double b_PPB = -0.6166;
  static const double c_PPB = 0.6703;
  static const double d_PPB = 0.5300;


  static const double A2P_s=1.798287e9;           //  Morton, Wu & Drake (2006)
  static const double A2P_t=177.58e0;             //  Lach & Pachuski (2001)
  static const double L_He_2Pt=1.690871466e7;     //  Drake & Morton (2007)
  static const double L_He_2st=1.5985597526e7;    //  Drake & Morton (2007)
  static const double L_He2St_ion=3.8454693845e6; //  Drake & Morton (2007)
  static const double sigma_He_2Ps=1.436289e-22;  //  Hummer & Storey (1998)
  static const double sigma_He_2Pt=1.484872e-22;  //  Hummer & Storey (1998)
  // Atomic data for HeI

  static const double AGauss1=-0.14;               //  Amplitude of 1st Gaussian
  static const double AGauss2=0.05;               //  Amplitude of 2nd Gaussian
  static const double zGauss1=7.28;               //  ln(1+z) of 1st Gaussian
  static const double zGauss2=6.75;               //  ln(1+z) of 2nd Gaussian
  static const double wGauss1=0.18;               //  Width of 1st Gaussian
  static const double wGauss2=0.33;               //  Width of 2nd Gaussian
  //  Gaussian fits for extra H physics (fit by Adam Moss)


  // the Verner and Ferland type fitting parameters for Helium
  //  fixed to match those in the SSS papers, and now correct
  double	a_VF = pow(10.0,(-16.744));
  double	b_VF = 0.711;
  double	T_0 = pow(10.0,0.477121); // 3K
  double	T_1 = pow(10.0,5.114);
  // fitting parameters for HeI triplets
  // (matches Hummer's table with <1% error for 10^2.8 < T/K < 10^4)
  double a_trip = pow(10.0, -16.306);
  double b_trip = 0.761;

  double	x_H = y[1];
  double	x_He = y[2];
  double	x = x_H + fHe * x_He;
  double	Tmat = y[3];


  double n = Nnow * pow(1.0 +z,3);
  double n_He = fHe * Nnow * pow(1.0 +z,3);
  double Trad = Tnow * (1.0 +z);

  double Hz = cosmos.tau2Hubble(cosmos.z2tau(z))*cosmos.cpm2sInv(); // general H(z)
  double dHdz = mdHdz->fastY(z)*cosmos.cpm2sInv();
  /* --  we don't need any of this, since we take the real values from cosmos
  //H = cosmos.H_0();
  double H0=cosmos.H_0();
  double bigH = 100.0e3/(1.0e6*3.0856775807e16);
  H0 = cosmos.h()*bigH;
  static const double G = 6.6742e-11;    // new value
  double a=7.565914e-16;
  double fnu = (cosmos.nuR()*7.0/8.0)*pow(4.0/11.0, 4.0/3.0);
  fnu = (21.0/8.0)*pow(4.0/11.0,4.0/3.0);
  double z_eq = (3.0*pow(H0*C,2)/(8.0*M_PI*G*a*(1.0+fnu)*pow(Tnow,4)))*OmegaT;
  Hz = H0 * sqrt(pow((1.0+z),4)/(1.0+z_eq)*OmegaT + OmegaT*pow(1.0+z,3)
       + cosmos.omega_v());
  dHdz = (pow(H0,2)/2.0/Hz)*(4.0*pow(1.0+z,3)/(1.0+z_eq)*OmegaT
         + 3.*OmegaT*pow(1.0+z,2));
  // -- */

  // Get the radiative rates using PPQ fit (identical to Hummer's table)
  double Rdown=1e-19*a_PPB*pow(Tmat/1e4,b_PPB) /(1.0+c_PPB*pow(Tmat/1e4,d_PPB));
  double Rup = Rdown * pow(CR*Tmat,1.5)*exp(-CDB/Tmat);

  // calculate He using a fit to a Verner & Ferland type formula
  double sq_0 = sqrt(Tmat/T_0);
  double sq_1 = sqrt(Tmat/T_1);
  //  typo here corrected by Wayne Hu and Savita Gahlaut
  double Rdown_He = a_VF/(sq_0*pow(1.0 +sq_0,1.0 -b_VF));
  Rdown_He = Rdown_He/pow(1.0+sq_1,1.0+b_VF);
  double Rup_He = Rdown_He*pow(CR*Tmat,1.5)*exp(-CDB_He/Tmat);
  Rup_He *= 4.0;     // statistical weights factor for HeI
  //  avoid overflow (pointed out by Jacques Roland)
  double He_Boltz;
  if (Bfact/Tmat > 680.0) He_Boltz = exp(680.0);
  else He_Boltz = exp(Bfact/Tmat);

  double K = CK/Hz;		// Peebles coefficient K=lambda_a^3/8piH

  // now deal with H and its fudges
  unsigned int Hswitch = mControl.recfastHSwitch;
  if (Hswitch==0) {
    K = CK/Hz;  // Peebles coefficient K=lambda_a^3/8piH
  } else {
    // fit a double Gaussian correction function
    K = CK/Hz*(1.0
          +AGauss1*exp(-pow((log(1.0+z)-zGauss1)/wGauss1, 2.))
          +AGauss2*exp(-pow((log(1.0+z)-zGauss2)/wGauss2, 2.)));
  }

  // add the HeI part, using same T_0 and T_1 values
  double Rdown_trip = a_trip/(sq_0*pow((1.0+sq_0),(1.0-b_trip)));
  Rdown_trip = Rdown_trip/(pow((1.0+sq_1),(1.0+b_trip)));
  double Rup_trip = Rdown_trip*exp(-h_P*C*L_He2St_ion/(k_B*Tmat));
  Rup_trip = Rup_trip*(pow((CR*Tmat),(1.50)))*(4.0/3.0);
  // last factor here is the statistical weight

  unsigned int Heswitch = mControl.recfastHeliumRecombination;
  HeliumRecombinationSwitch Heflag;
  // try to avoid "NaN" when x_He gets too small
  if ((x_He<5e-9) || (x_He>0.980)) {
    Heflag = 0;
  } else {
    Heflag = Heswitch;
  }

  double K_He;
  if (Heflag==0) {                    // use Peebles coeff. for He
    K_He = CK_He/Hz;
  } else {           // for Heflag>0  // use Sobolev escape probability
    const double tauHe_s = A2P_s*CK_He*3.0*n_He*(1.0-x_He)/Hz;
    const double pHe_s = (1.0 - exp(-tauHe_s))/tauHe_s;
    K_He = 1.0/(A2P_s*pHe_s*3.0*n_He*(1.0-x_He));
    //  smoother criterion here from Antony Lewis & Chad Fendt
      if (((Heflag==2) || (Heflag>=5)) && (x_H<0.9999999)) {
        // use fitting formula for continuum opacity of H
        // first get the Doppler width parameter
        double Doppler = 2.0*k_B*Tmat/(m_H*not4*C*C);
        Doppler = C*L_He_2p*sqrt(Doppler);
        double gamma_2Ps = 3.0*A2P_s*fHe*(1.0-x_He)*C*C
                    /(sqrt(M_PI)*sigma_He_2Ps*8.0*M_PI*Doppler*(1.0-x_H))
                    /(pow(C*L_He_2p, 2.0));
        double pb = 0.36;    // value from KIV (2007)
        double qb = b_He;
        // calculate AHcon, the value of A*p_(con,H) for H continuum opacity
        double AHcon = A2P_s/(1.0+pb*(pow(gamma_2Ps, qb)));
        K_He=1.0/((A2P_s*pHe_s+AHcon)*3.0*n_He*(1.0-x_He));
      }
    if (Heflag>=3) {     // include triplet effects
      double tauHe_t = A2P_t*n_He*(1.0-x_He)*3.0;
      tauHe_t = tauHe_t /(8.0*M_PI*Hz*pow(L_He_2Pt, 3.0));
      const double pHe_t = (1.0 - exp(-tauHe_t))/tauHe_t;
      const double CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B;
      if ((Heflag==3) || (Heflag==5) || (x_H>0.99999)) {
        // no H cont. effect
        CfHe_t = A2P_t*pHe_t*exp(-CL_PSt/Tmat);
        CfHe_t = CfHe_t/(Rup_trip+CfHe_t);     // "C" factor for triplets
      } else {   // include H cont. effect
        double Doppler = 2.0*k_B*Tmat/(m_H*not4*C*C);
        Doppler = C*L_He_2Pt*sqrt(Doppler);
        double gamma_2Pt = 3.0*A2P_t*fHe*(1.0-x_He)*C*C
                    /(sqrt(M_PI)*sigma_He_2Pt*8.0*M_PI*Doppler*(1.0-x_H))
                    /(pow((C*L_He_2Pt),2.0));
        // use the fitting parameters from KIV (2007) in this case
        double pb = 0.66;
        double qb = 0.9;
        double AHcon = A2P_t/(1.0+pb*pow(gamma_2Pt, qb))/3.0;
        CfHe_t = (A2P_t*pHe_t+AHcon)*exp(-CL_PSt/Tmat);
        CfHe_t = CfHe_t/(Rup_trip+CfHe_t);  // "C" factor for triplets
      }
    }
  }

  //	Estimates of Thomson scattering time and Hubble time
  double timeTh=(1.0/(CT*pow(Trad,4)))*(1.0 +x+fHe)/x;	// Thomson time
  // better for hubble time :
  //double timeH=2.0/(3.0*H*pow(1.0+z,1.5)); // Hubble time
  double timeH = 2.0/ (3.0* cosmos.tau2Hubble(cosmos.z2tau(z)) * cosmos.cpm2sInv() ); // in s^-1

  //	calculate the derivatives
  //	turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
  //	(clunky, but seems to work)
  if (x_H > 0.99) {                  	// don't change at all
    f[1] = 0;
  } else if (x_H > 0.985) {          // use Saha rate for Hydrogen
    f[1] = (x*x_H*n*Rdown - Rup*(1.0-x_H)*exp(-CL/Tmat))/(Hz*(1.0+z));
    /*
    // cout << "z: " << z << "  y[1] : " << y[1] << "  f[1]: " << f[1] << endl;
    // for interest, calculate the correction factor compared to Saha
    // (without the fudge)
    const double factor=(1.0 + K*Lambda*n*(1.0-x_H))
           /(Hz*(1.0+z)*(1.0+K*Lambda*n*(1.0-x)
           +K*Rup*n*(1.0-x)));
    */
  } else {                             // use full rate for H
    f[1] = ((x*x_H*n*Rdown - Rup*(1.0-x_H)*exp(-CL/Tmat)) *(1.0 + K*Lambda*n*(1.0-x_H)));
    f[1] /=  Hz*(1.0+z)*(1.0/fu+K*Lambda*n*(1.0-x_H)/fu +K*Rup*n*(1.0-x_H));  // was (1-x) typo thank to Jens Chluba
    // compare this typo also to correct version in RecfastAlpha
  }
  //   turn off the He once it is small
  if (x_He < 1e-15) f[2]=0;
  else {
    f[2] = ((x*x_He*n*Rdown_He - Rup_He*(1.0-x_He)*exp(-CL_He/Tmat)) *(1.0+ K_He*Lambda_He*n_He*(1.0-x_He)*He_Boltz));
    f[2] /= Hz*(1.0+z) * (1.0 + K_He*(Lambda_He+Rup_He)*n_He*(1.0-x_He) *He_Boltz); 
    // Modification to HeI recombination including channel via triplets
    if (Heflag>=3) {
      f[2] = f[2]+ (x*x_He*n*Rdown_trip
             - (1.0-x_He)*3.0*Rup_trip*exp(-h_P*C*L_He_2st/(k_B*Tmat)))
             *CfHe_t/(Hz*(1.0+z));
    }
  }
  //    follow the matter temperature once it has a chance of diverging
  // M.Doran:: I'm not happy with this, as it still bears the possibility
  // of switching over a little too late or risking oscillations in T_b when 
  // switching on the coupling, such that it can run to it's "true" value
  // However, we drop the first 4-5 z's after the switch and therefore
  // we will most propably not see many oscillations. In addition,
  // at the time this is happening, the baryon sound speed is already
  // negligible and we don't have to worry for the soundspeed.
  // Yet,  I still would prefer a "perfect"  T_matter evolution...
  // G. Robbers 01/2009: comment superseeded by inclusion of epsilon term?                      
  if (timeTh < H_frac*timeH) {
    // f[3]=Tmat/(1.0+z); // Tmat follows Trad
    // additional term to smooth transition to Tmat evolution,
    // (suggested by Adam Moss)
    const double epsilon = Hz*(1.0+x+fHe)/(CT*pow(Trad,3)*x);
    f[3] = Tnow
           + epsilon*((1.0+fHe)/(1.0+fHe+x))*((f[1]+fHe*f[2])/x)
           - epsilon* dHdz/Hz + 3.0*epsilon/(1.0+z);       
  } else {
    if (SwitchTbEvolution == 0) SwitchTbEvolution = z;
    f[3]= CT * pow(Trad,4) * x / (1.0+x+fHe) * (Tmat-Trad) / (Hz*(1.0 +z)) + 2.0*Tmat/(1.0+z);
  }

  //  cout << "curioisity: " << z << "  :: " <<  (1.0 + K*Lambda*n*(1.0-x_H)) / (1.0/fu+K*Lambda*n*(1.0-x)/fu +K*Rup*n*(1.0-x)) << endl;


}          
