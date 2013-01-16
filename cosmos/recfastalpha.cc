#include "recfastalpha.h"
#include "cosmos.h"
#include <fstream>
#include "minmax.h"



void RecfastAlpha::compute() {
  // ofstream helium("helium.dat");
  Crossed98 = false;
  CrossedHelium = false;
  C  = 2.99792458e8;
  k_B = 1.380658e-23;
  static double h_P = 6.6260755e-34;
  static double m_e = 9.1093897e-31;
  m_H = 1.673725e-27;

  static double sigma = 6.6524616e-29;
  static double a = 7.565914e-16;
  
  static double G = 6.67259e-11;
  static double L_H_ion = 10967877.37;
  static double L_H_alpha = 8225916.453;
  static double L_He1_ion = 19831077.2;
  static double L_He2_ion = 43890888.63;
  static double L_He_2s = 16627743.4;
  static double L_He_2p = 17113489.1;

  double zinitial = 1e4;
  double z = zinitial;
  double zfinal = 0;

  SwitchTbEvolution = 0;  // this will be the redshift at which ion() switches from tight for T_b to non-tight

  OmegaB = cosmos.omega_b();
  OmegaT = cosmos.omega_m(); // total matter 
  //  double omegac = cosmos.omega_cdm();
  //double omegav = cosmos.omega_v();
  // double h = cosmos.h();
  // double tcmb = cosmos.T_cmb();
  Tnow = cosmos.T_cmb();
  double Yp = cosmos.Y_he();

  double mu_H = 1.0/(1.0-Yp);     //Mass per H atom
  mu_T = 4.0/(4.0-3.0*Yp); //Mass per atom
  fHe = Yp/(4.0*(1.0-Yp)); // n_He_tot / n_H_tot

  H = cosmos.H_0();  // in seconds (that's originals recfast's HO)
  
  //  today's number density of hydrogen
  Nnow = 3*H*H*OmegaB/(8.0*M_PI*G*mu_H*m_H);
	
  // set up some constants so they don't have to be calculated later
  
  double Lalpha = 1.0 / L_H_alpha;
  double Lalpha_He = 1.0 / L_He_2p;
  double DeltaB = h_P*C*(L_H_ion-L_H_alpha);
  CDB = DeltaB/k_B;
  double DeltaB_He = h_P*C*(L_He1_ion-L_He_2s); // 2s, not 2p
  CDB_He = DeltaB_He/k_B;
  CB1 = h_P*C*L_H_ion/k_B;
  CB1_He1 = h_P*C*L_He1_ion/k_B; // ionization for HeI
  CB1_He2 = h_P*C*L_He2_ion/k_B; // ionization for HeII
  CR = 2.0*M_PI*(m_e/h_P)*(k_B/h_P);
  CK = pow(Lalpha,3)/(8.0*M_PI);
  CK_He = pow(Lalpha_He,3)/(8.0*M_PI);
  CL = C*h_P/(k_B*Lalpha);
  CL_He = C*h_P/(k_B/L_He_2s); // comes from det.bal. of 2s-1s
  CT = (8.0/3.0)*(sigma/(m_e*C))*a;
  Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B;

  double y[4]; // 4 dim for 3 variables for odeint(); 

  //    Matter departs from radiation when t(Th) > H_frac * t(H)
  //    choose some safely small number
  H_frac = 1e-3;
      
  //    Fudge factor to approximate for low z out of equilibrium effect
  fu=1.14;
      
  //   Set initial matter temperature
  y[3] = Tnow*(1.0+z);      //Initial rad. & mat. temperature
  
  getInit(z); // set x_H0, x_He0 and x0
  y[1] = x_H0;
  y[2] = x_He0;
    
  //     OK that's the initial conditions, now start writing output file
  int Nz=10000;


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
    } else {
      x_H0 = 1;
      x_He0 = 1;
      rhs = exp( 1.5 * log(CR*Tnow/(1.0+z)) - CB1_He2*(1+2*delAlpha(z))/(Tnow*(1.0+z)) ) / Nnow;
      rhs *= 1.0;       // ratio of g's is 1 for He++ <-> He+
    
      if (rhs > 1e-7) {  // so if there is still He-II to recombine
	x0 = 0.5 * ( sqrt( pow( rhs-1.0 -fHe,2.0)  + 4.0*(1.0 +2.0*fHe)*rhs) - (rhs-1.0-fHe) );
	y[1] = x_H0;
	y[2] = x_He0;
	y[3] = Tnow*(1.0+z);
	//	cout << "z: " <<z << " rhs: " << rhs << " x0: " << x0 << endl;
      } else {  // now, let go and evolve everything, typically at redshifts 3500-5000 
	hnext = Miscmath::odeint(y,3,zstart,zend,1e-9,hnext,0, (moDerivs)&RecfastAlpha::ion,*this,false,Miscmath::rkutta);
	
	x0 = y[1] + fHe*y[2];
	ion(zend,y,dy); // to get derivative of Tb
	if (SwitchTbEvolution > 0) {  // so we have switched already 
	  if (SwitchTbEvolution - zend > 5) { // drop the first z's
	    dTbTb = dy[3] / y[3]; 
	  } // i.e. do nothing if it just switched, keep canonical dTbTb value
	}
      }
    }
    //Trad = Tnow * (1.0+zend);
    // Tmat = y[3];
    // x_H = y[1];
    //x_He = y[2];         
    //x = x0;
    noteX->set(-zend, x0);  // set the spline 
    noteTb->set(y[3]); // y[3] holds baryon temperature
    noteTbDeriv->set( -(1+zend) * dTbTb); //  d ln Tb / d ln a
    // eqn (68) of Ma & Bertschinger, the + instead of minus comes from 
    // d ln Tb_ dln a = -(1+z) d ln Tb / dz 
    double cs2 = y[3] * BoltzmannOverWeightAndC2() * (1.0 +  1.0/3.0 *(1+zend)*dTbTb);
    noteCs2->set(cs2);
  }  // end the loop over the z's
  noteX->arm(Spline::all);
  noteX->dump("recombin");  // if you would ike to stay informed
}


/*!
  (Re-)Implement your varying alpha behaviour here
  This one is just for fun and give 1% change at recombination
*/
double RecfastAlpha::fineStructureConstant(double z) {
  double rel = 0.01; 
  double Z = min(z,1100.0);
 
  return 1.0/137.0*( 1.0 + Z/1100.0 *rel);
} 

void RecfastAlpha::getInit(double z) {
  cout << "RecfastAlpha - WARNING: this class has not been updated to recfast 1.4.2 like the standard Recfast" << endl;
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
    double alpha = fineStructureConstant(z);
    double dalpha = fineStructureConstant(z) - fineStructureConstant(0);       
    rhs = exp( 1.5 * log(CR*Tnow/(1.0+z)) - CB1*(1.0 + 2.0*dalpha/alpha)/(Tnow*(1.0+z)) ) / Nnow;

    x_H0 = 0.5 * (sqrt(rhs*rhs +4.0*rhs ) - rhs );
    x_He0 = 0;
    x0 = x_H0;
  }
}


void RecfastAlpha::ion(const double z, const double *y, double *f) {
  bool print =false  ;

  static double Lambda = 8.2245809;
  static double Lambda_He = 51.3;

  //      the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen	
  static double a_PPB = 4.309;
  static double b_PPB = -0.6166;
  static double c_PPB = 0.6703;
  static double d_PPB = 0.5300;

  // the Verner and Ferland type fitting parameters for Helium
  //  fixed to match those in the SSS papers, and now correct
  double	a_VF = pow(10.0,(-16.744));
  double	b_VF = 0.711;
  double	T_0 = pow(10.0,0.477121); // 3K
  double	T_1 = pow(10.0,5.114); // 
    
  double	x_H = y[1];
  double	x_He = y[2]; 
  double	x = x_H + fHe * x_He;
  double	Tmat = y[3];
  
  double n = Nnow * pow(1.0 +z,3);
  double n_He = fHe * Nnow * pow(1.0 +z,3);
  double Trad = Tnow * (1.0 +z);

  //H = cosmos.H_0();
  double Hz = cosmos.tau2Hubble(cosmos.z2tau(z)) *cosmos.cpm2sInv()  ; // general H(z) 
  
  // Get the radiative rates using PPQ fit (identical to Hummer's table)
  double Rdown=1e-19*a_PPB*pow(Tmat/1e4,b_PPB) /(1.0+c_PPB*pow(Tmat/1e4,d_PPB));
  //  double Rup = Rdown * pow(CR*Tmat,1.5)*exp(-CDB/Tmat);
			      

  // calculate He using a fit to a Verner & Ferland type formula
  double sq_0 = sqrt(Tmat/T_0);
  double sq_1 = sqrt(Tmat/T_1);
  //  typo here corrected by Wayne Hu and Savita Gahlaut
  double Rdown_He = a_VF/(sq_0*pow(1.0 +sq_0,1.0 -b_VF));
  Rdown_He = Rdown_He/pow(1.0+sq_1,1.0+b_VF);
  double Rup_He = Rdown_He*pow(CR*Tmat,1.5)*exp(-CDB_He/Tmat);
  Rup_He *= 4.0;	 // statistical weights factor for HeI
  //	Avoid overflow (pointed out by Jacques Roland)
  double He_Boltz;
  if (Bfact/Tmat > 680.0) He_Boltz = exp(680.0);
  else He_Boltz = exp(Bfact/Tmat);
  
  double K = CK/Hz;		// Peebles coefficient K=lambda_a^3/8piH
  double K_He = CK_He/Hz;		// Peebles coefficient for Helium
  
  //	Estimates of Thomson scattering time and Hubble time
  double timeTh=(1.0/(CT*pow(Trad,4)))*(1.0 +x+fHe)/x;	// Thomson time
  // better for hubble time :
  //double timeH=2.0/(3.0*H*pow(1.0+z,1.5)); // Hubble time
  double timeH = 2.0/ (3.0* cosmos.tau2Hubble(cosmos.z2tau(z)) * cosmos.cpm2sInv() ); // in s^-1


   // 
  // From here on, we calculate some quantities needed for varying fine structure constant
  // 

  // calculate the effect of a changing fine structure constant
  double alpha = fineStructureConstant(z);
  double dalpha = fineStructureConstant(z) - fineStructureConstant(0); 
  double del = dalpha / alpha;

  // first, we calculate d ln Rdown / d T from the explicit formula above
  double t = Tmat * 1e-4; // the fit is expressed in terms of this quantity 
  double dlnRdowndT = 1e-4*(b_PPB/t - c_PPB * d_PPB * pow(t,d_PPB-2) / (1 + c_PPB*pow(t,d_PPB))); 

  // as the fitting formulae is for a table of Hummer and this table only goes down to T = 10 K,
  // I will cut the increasingly unphysical dlnRdowndT at about T = 10 K
  // if this is not done, then for really large changes in the fine structure constant,
  // you will end up with dln R / dT > 1, messing up peebles correction factor below
  // actually, as recombination is over, I would rather like to be on the safe side and
  // cut it out say at 50 K, cause then for all sensible models, dlnRdowndT is still < 1e-1
  // and we are fine...
  dlnRdowndT *= exp(- pow(50/ Tmat,3));

  // now we use the expression d / d alpha R = 2/alpha *(R - T* d R / d T) see M. Kaplinghat et al (1999), eqn (7)
  // to relate d / d alpha R to the T dependence
  double dlnRdown = 2.0*(1.0 - Tmat * dlnRdowndT)*del;  // d (ln Rdown)
  // get the effective Rup called effective beta in Kalplinghat et. at 
  //  (this is the full beta with the change in Rdown factored out and stuck in a changed correction factor)
  double RupEffective =  Rdown * pow(CR*Tmat,1.5)*exp(-CDB/Tmat - CB1/Tmat*(2*del + del*del) ); 
  
  // He - varying alpha

  double Rup_HeEffective = 4*Rdown_He*pow(CR*Tmat,1.5)*exp(-CDB_He/Tmat - CB1_He1/Tmat*(2*del + del*del) ); 
  
  //
  // End of the fine structure constant stuff
  // 


  //
  // M. Doran:  
  //
  // We use a different strategy to determine the onset of recombination than the original RECFAST
  // 
  // We do the following:
  // (1) The ratio r = Rup * exp(-CL/Tmat) / (n*Rdown) 
  //       Determines wether the strength to re-ionize is in principle stronger than the one to recombine
  //       
  // (2)  If the Peebles Correction factor C times the rate n*Rdown is >> Hubble 
  //       then our rates are sufficiently large to reach equillibrium
  //
  // (3)  So if r > 1 and C*n*Rdown >> H,
  //       we determine x_e by solving for stationary x_e evolution.
  //  
  // Please note that in (2) we need to know C and for this we need to know x_e. So in
  // (2), we *assume* rapid collisions, take step (3) and then go back and evaluate the rate vs Hubble
  
  double rup = RupEffective*(1.0 + del);  // for the correction factor, full changed Rup
  double k = K*(1.0 - 6*del);  // for the correction factor: change K
  double lambda = Lambda*(1.0 + 8*del); // for the correction factor: change Lambda
  
  double r  = RupEffective*exp(-CL/Tmat)/(n*Rdown);
  double dx = fHe*x_He;    
  double balance_x =  -(r+dx)/2 + sqrt(0.25*(r+dx)*(r+dx) + r);
  double one_minus_xp = 1.0 - balance_x;
  
 
  // here is the *good* formulae. This is much more precise than taking the root and subtracting
  // typical accuracy at r =10^4:  is  10^-13 !
  if (r > 1e3) {
    balance_x = 1.0 - (1.0+dx)/r + (2.0+3*dx + dx*dx)/(r*r) - (5.0+10*dx + 6*dx*dx)/(r*r*r) + (14 +35*dx + 30*dx*dx)/(r*r*r*r);
    one_minus_xp = (1.0+dx)/r - (2.0+3*dx + dx*dx)/(r*r) + (5.0 + 10*dx + 6*dx*dx)/(r*r*r) - (14 + 35*dx + 30*dx*dx)/(r*r*r*r);
  }
  
  double Cstar = (1.0 + k*lambda*n*one_minus_xp);
  Cstar /= Hz*(1.0+z)*(1.0/fu+k*lambda*n*one_minus_xp/fu +k*rup*n*one_minus_xp); // some fudge factors cancelled with nominator's 
  double Cprime = Cstar*(1.0 + dlnRdown); // absorb the change in Rdown just like Kaplinghat
  //  double absdown = Cprime*Rdown*n*Hz*(1+z);
  
  bool tight = true;
  double PotentialForceUp =  Cprime*RupEffective*exp(-CL/Tmat)*Hz*(1+z);
  if (PotentialForceUp / Hz < 10) tight = false;  // is up rate fast enough compared to hubble ?
  double TypicalBracket = max(n*Rdown, RupEffective*exp(-CL/Tmat)); // max typical bracket value  
  if (TypicalBracket*Cprime < 1) tight=false; // cause our reasoning is only ok, if the bracket needs to cancel...
  
  if (tight && !Crossed98) {
    if (print) cout << "@@ ";
    x_H = balance_x;
    f[1] = -(x_H-y[1]); // adjust
  } else {  // otherwise just use the full evolution equation, this time with the actual x_H  ...
    Crossed98 = true;
    if (print)    cout << "##";
    one_minus_xp = 1.0 - x_H;
    double Cstar = (1.0 + k*lambda*n*one_minus_xp);
    Cstar /= Hz*(1.0+z)*(1.0/fu+k*lambda*n*one_minus_xp/fu +k*rup*n*one_minus_xp); // some fudge factors cancelled with nominator's 
    double Cprime = Cstar*(1.0 + dlnRdown); // absorb the change in Rdown just like Kaplinghat
    f[1] = ((x_H+dx)*x_H*n*Rdown - RupEffective*one_minus_xp*exp(-CL/Tmat)) * Cprime;

    if (Cprime < 0) {
      cout << "dlnRdown: " << dlnRdown << endl;
      cout << "k*Lambda*n: " << k*lambda*n << "   one_minus_xp : " << one_minus_xp << endl;
      cout << "xp: " << x_H << "   k*rup: " << k*rup << endl;
      throw Bad_Error("RecfastAlpha:: Cprime <0, this is numerically unstable!!");
    }   
  }

  //   turn off the He once it is small
  if (x_He < 1e-15) f[2]=0; 
  else {
    // "das selbe in gruen" :: "the same procedure as above :-)"
    // We copy our treatment for Hydrogen, just slight changes, for instance fHe 
    // enters in "r"
    
    double rup = Rup_HeEffective*(1.0 + del);  // for the correction factor, full changed Rup
    double k = K_He*(1.0 - 6*del);  // for the correction factor: change K
    double lambda = Lambda_He*(1.0 + 8*del); // for the correction factor: change Lambda
    
    double r  = Rup_HeEffective*exp(-CL_He/Tmat)/(n*Rdown_He*fHe);
    double dx = x_H / fHe; 
    double balance_x =  -(r+dx)/2 + sqrt(0.25*(r+dx)*(r+dx) + r);  // equillibrium x 
    double one_minus_xhe = 1.0 - balance_x;
    
    // here is the *good* formulae. This is much more precise than taking the root and subtracting
    // typical accuracy at r =10^4:  is  10^-11 (even though dx is so bitchy for helium) !
    
    // careful: For helium this is much more difficult, as dx corresponds to x_H / fHe
    // and therefore dx is approx 12 and therefore >> 1 !!!

    if (r > 1e3) {
      balance_x = 1.0 - (1.0+dx)/r + (2.0+3*dx + dx*dx)/(r*r) - (5.0+10*dx + 6*dx*dx + dx*dx*dx)/(r*r*r) + (14 +35*dx + 30*dx*dx + 10*dx*dx*dx + dx*dx*dx*dx)/(r*r*r*r);
      one_minus_xhe = (1.0+dx)/r - (2.0+3*dx + dx*dx)/(r*r) + (5.0 + 10*dx + 6*dx*dx + dx*dx*dx)/(r*r*r) - (14 + 35*dx + 30*dx*dx +10*dx*dx*dx + dx*dx*dx*dx)/(r*r*r*r);
    }

    

    double dlnRdown_He = 0.0;
    
    double Cstar = 1.0+ k*lambda*n_He*one_minus_xhe*He_Boltz;
    Cstar /= Hz*(1.0+z) * (1.0 + k*(lambda +rup)*n_He*one_minus_xhe *He_Boltz);
    double Cprime = Cstar*(1.0 + dlnRdown_He); // absorb the change in Rdown just like Kaplinghat

    bool tight = true;
    double PotentialForceUp =  Cprime*Rup_HeEffective*exp(-CL_He/Tmat)*Hz*(1+z);
    if (PotentialForceUp / Hz < 10) tight = false;  // is up rate fast enough compared to hubble ?
    double TypicalBracket = max(n*Rdown_He, Rup_HeEffective*exp(-CL_He/Tmat)); // max typical bracket value  
    if (TypicalBracket*Cprime < 1) tight=false; // cause our reasoning is only ok, if the bracket needs to cancel...

    if (! tight  && !CrossedHelium && print) {
      if (PotentialForceUp/Hz < 10) cout << endl << "Force up triggered it"; else cout << endl <<  "Typical bracket triggered" << endl;
      cout << "at: y[2] : " << y[2] << " f[2]: " << f[2] << endl; 
    }

    //    if (one_minus_xhe < 1e-3) tight = true;  // for really small values always tight
    
    if (print && false) {
    cout << " absup/H: " << PotentialForceUp / Hz << " TB*C: " << TypicalBracket*Cprime << "  TB: " << TypicalBracket;
    cout << " r: " << r;
    }

    if (tight && !CrossedHelium) {
      x_He = balance_x;
      f[2] = -(x_He-y[2]); // adjust
      if (print) cout << " +++ "; // Cprime: " <<Cprime;
    } else {
      if (print) cout << " ---- " ;
      CrossedHelium = true;
      one_minus_xhe = 1.0 - x_He;
      
      double Cstar = 1.0+ k*lambda*n_He*one_minus_xhe*He_Boltz;
      Cstar /= Hz*(1.0+z) * (1.0 + k*(lambda +rup)*n_He*one_minus_xhe *He_Boltz);
      double Cprime = Cstar*(1.0 + dlnRdown_He); // absorb the change in Rdown just like Kaplinghat
      
      f[2] = ((x_H + fHe*x_He)  *x_He*n*Rdown_He - Rup_HeEffective*one_minus_xhe*exp(-CL_He/Tmat))*Cprime;
    }
  } 
    
  
  if (print) {
    cout << "  z: " << z << " y[1]: " << y[1] << " f[1]: " << f[1];
    cout << "  ||| y[2]: " << y[2] << "  f[2] " << f[2] << endl;
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
  if (timeTh < H_frac*timeH) 
    f[3]=Tmat/(1.0+z); // Tmat follows Trad
  else {
    if (SwitchTbEvolution == 0) SwitchTbEvolution = z;
    f[3]= CT * pow(Trad,4) * x / (1.0+x+fHe) * (Tmat-Trad) / (Hz*(1.0 +z)) + 2.0*Tmat/(1.0+z);
  }            
  
}          
