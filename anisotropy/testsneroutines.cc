#include "global.h"
#include "cosmos.h"
#include "distcosmos.h"
#include "distcelestine.h"
#include "analyzethis.h"

/*#include "snedata.h"*/

#include <iostream>
#include <cstdio>

using namespace std;

//Test program for baseCosmos type distances
int main(int argc, char* argv[]) {

  double h = 0.65;
  double ov = 0.8;

  Cosmos cosmos;

  cosmos.seth(h);
  cosmos.setOmega_vacuum(ov);
  cosmos.setOmega_cdm( 1.0-cosmos.omega_b() - 
		       cosmos.omega_v() - cosmos.omega_nu() -
		       cosmos.omega_g() );

  distCosmos lcosmos;
  lcosmos.setZMax(2.0);
  lcosmos.seth(h);
  lcosmos.setOmega_m(0.733);
  lcosmos.setOmega_v(1.23);
  lcosmos.setOmegaH2_rel( 0.0 );

  distCelestine ccosmos;
  ccosmos.setZMax(2.0);
  ccosmos.seth(h);
  ccosmos.setOmega_m(0.7);
  ccosmos.setOmega_v(0.3);
  ccosmos.setOmega_q(1.0);
  ccosmos.setOmegaH2_rel( 0.0 );
  ccosmos.setW0(-1);
  ccosmos.setWa(0.0);

  AnalyzeThis ai;

  try {
    printf("Cosmos\n");
    cosmos.history(false);
     
    cout << "For Universe with:\n";
    cout << " Hubble h: " << cosmos.h() << endl;
    cout << " Omega_rel: " << cosmos.omega_relativistic() << endl;
    cout << " Omega_matter: " << cosmos.omega_m() << endl;
    cout << " Omega_vacuum: " << cosmos.omega_v() << endl;

    printf(" Lum dist at z = 0.00: %10.5lf Mpc\n",
	   cosmos.luminosityDistance(0.00));
    printf(" Lum dist at z = 0.01: %10.5lf Mpc\n",
	   cosmos.luminosityDistance(0.01));
    printf(" Lum dist at z = 0.10: %10.5lf Mpc\n",
	   cosmos.luminosityDistance(0.1));
    printf(" Lum dist at z = 0.40: %10.5lf Mpc\n",
	   cosmos.luminosityDistance(0.4));
    printf(" Lum dist at z = 0.54: %10.5lf Mpc\n",
	   cosmos.luminosityDistance(0.54));
    printf(" Lum dist at z = 0.92: %10.5lf Mpc\n",
	   cosmos.luminosityDistance(0.92));
    printf(" Lum dist at z = 1.00: %10.5lf Mpc\n",
	   cosmos.luminosityDistance(1.0));
    printf(" Deriv Lum dist at z = 1.00: %10.5lf Mpc\n",
	   cosmos.luminosityDistanceDeriv(1.0));

    printf("distCosmos\n");
    lcosmos.setZMax(1.1);
    lcosmos.history(false);

    cout << "For Universe with:\n";
    cout << " Hubble h: " << lcosmos.h() << endl;
    cout << " Omega_rel: " << lcosmos.omega_relativistic() << endl;
    cout << " Omega_matter: " << lcosmos.omega_m() << endl;
    cout << " Omega_vacuum: " << lcosmos.omega_v() << endl;
     
    printf(" Lum dist at z = 0.00: %10.5lf Mpc\n",
	   lcosmos.luminosityDistance(0.00));
    printf(" Lum dist at z = 0.01: %10.5lf Mpc\n",
	   lcosmos.luminosityDistance(0.01));
    printf(" Lum dist at z = 0.10: %10.5lf Mpc\n",
	   lcosmos.luminosityDistance(0.1));
    printf(" Lum dist at z = 0.40: %10.5lf Mpc\n",
	   lcosmos.luminosityDistance(0.4));
    printf(" Lum dist at z = 0.54: %10.5lf Mpc\n",
	   lcosmos.luminosityDistance(0.54));
    printf(" Lum dist at z = 0.92: %10.5lf Mpc\n",
	   lcosmos.luminosityDistance(0.92));
    printf(" Lum dist at z = 1.00: %10.5lf Mpc\n",
	   lcosmos.luminosityDistance(1.0));
    printf(" Deriv Lum dist at z = 1.00: %10.5lf Mpc\n",
	   lcosmos.luminosityDistanceDeriv(1.0));

    if ( lcosmos.isUniverseRebounding() ) cout << "Rebounding Universe\n";

    printf("distCelestine\n");
    ccosmos.history(false);

    cout << "For Universe with:\n";
    cout << " Hubble h: " << ccosmos.h() << endl;
    cout << " Omega_rel: " << ccosmos.omega_relativistic() << endl;
    cout << " Omega_matter: " << ccosmos.omega_m() << endl;
    cout << " Omega_vacuum: " << ccosmos.omega_v() << endl;
    cout << " Omega_q: " << ccosmos.omega_q() << endl;
    cout << " w_0: " << ccosmos.w0() << " w_a: " << ccosmos.wa() << endl;
    
    printf(" Lum dist at z = 0.00: %10.5lf Mpc\n",
	   ccosmos.luminosityDistance(0.00));
    printf(" Lum dist at z = 0.01: %10.5lf Mpc\n",
	   ccosmos.luminosityDistance(0.01));
    printf(" Lum dist at z = 0.10: %10.5lf Mpc\n",
	   ccosmos.luminosityDistance(0.1));
    printf(" Lum dist at z = 0.40: %10.5lf Mpc\n",
	   ccosmos.luminosityDistance(0.4));
    printf(" Lum dist at z = 0.54: %10.5lf Mpc\n",
	   ccosmos.luminosityDistance(0.54));
    printf(" Lum dist at z = 0.92: %10.5lf Mpc\n",
	   ccosmos.luminosityDistance(0.92));
    printf(" Lum dist at z = 1.00: %10.5lf Mpc\n",
	   ccosmos.luminosityDistance(1.0));
    printf(" Deriv Lum dist at z = 1.00: %10.5lf Mpc\n",
	   ccosmos.luminosityDistanceDeriv(1.0));

    double alpha = 1.2773;

    /*
    cout << "Extended test\n";

    //Try getting a chisquare.  Here goes nothing...
    cout << "Trying Lowe Knop SNe without extinction correction\n";
    cout << "Chisquare: " << ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::Lowe, 
					    0.0, 0.17) 
	 << endl;
    cout << "Estimated scriptM: " << 
      ai.estimateScriptMKnop03(lcosmos, alpha, AnalyzeThis::Lowe, 0.0, 0.17)
	 << endl;

    cout << "Trying Riess04 gold sample:\n";
    cout << "Old chisquare: " << ai.Sn1aRiess04(lcosmos,true) << endl;
    cout << "New Chisquare: " << ai.SNIaRiess04(lcosmos,true) << endl;

    cout << "Trying Tonry03 172 sample:\n";
    cout << "Old chisquare: " << ai.Sn1aTonry(lcosmos,true,false) << endl;
    cout << "New chisquare: " << ai.SNIaTonry03(lcosmos) << endl;

    cout << "Trying Barris04 sample:\n";
    cout << "Old chisquare: " << ai.Sn1aTonry(lcosmos,true,true) << endl;
    cout << "New chisquare: " << ai.SNIaBarris04(lcosmos) << endl;

    cout << "Trying SDSS Baryon Acoustic Peak:\n";
    cout << "distCosmos Chisquare: " << ai.SDSS_BAP_A_chiSquared(lcosmos) 
	 << endl;
    cout << "Cosmos Chisquare: " << ai.SDSS_BAP_A_chiSquared(cosmos) 
	 << endl;
    */


    //Now try a variety of cosmological parameters on Knop03 Data
    cout << "========================================================" << endl;
    double beta = 0.0;
    cout << "Knop03 Low Extinction primary subset\n";
    cout << "Omega_m: " << lcosmos.omega_m() << " Omega_L: " <<
      lcosmos.omega_v() << " Alpha: " << alpha << " Beta: " << beta << endl;
    cout << "Lowe Knop chi^2, no extcorr: " << 
      ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::K03Lowe, beta,0.17) << endl;
    
    cout << "ScriptM estimate: " << 
      ai.estimateScriptMKnop03( lcosmos, alpha, AnalyzeThis::K03Lowe, beta, 
				0.17) << endl;

    lcosmos.setOmega_m(0.5);
    lcosmos.setOmega_v(0.9);
    lcosmos.reset(); lcosmos.history();
    cout << "Omega_m: " << lcosmos.omega_m() << " Omega_L: " <<
      lcosmos.omega_v() << " Alpha: " << alpha << " Beta: " << beta << endl;
    cout << "Lowe Knop chi^2, no extcorr: " << 
      ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::K03Lowe, beta,0.17) << endl;

    lcosmos.setOmega_m(0.6);
    lcosmos.setOmega_v(1.4);
    lcosmos.reset(); lcosmos.history();
    cout << "Omega_m: " << lcosmos.omega_m() << " Omega_L: " <<
      lcosmos.omega_v() << " Alpha: " << alpha  << " Beta: " << beta << endl;
    cout << "Lowe Knop chi^2, no extcorr: " << 
      ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::K03Lowe, beta,0.17) << endl;

    //With extinction correction
    beta = 4.1; //Extinction correction
    cout << "Omega_m: " << lcosmos.omega_m() << " Omega_L: " <<
      lcosmos.omega_v() << " Alpha: " << alpha  << " Beta: " << beta << endl;
    cout << "Lowe Knop chi^2, extcorr: " << 
      ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::K03Lowe, beta, 0.11) << endl;
    
    lcosmos.setOmega_m(0.733);
    lcosmos.setOmega_v(1.23);
    lcosmos.reset(); lcosmos.history();
    cout << "Omega_m: " << lcosmos.omega_m() << " Omega_L: " <<
      lcosmos.omega_v() << " Alpha: " << alpha  << " Beta: " << beta << endl;
    cout << "Lowe Knop chi^2, extcorr: " << 
      ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::K03Lowe, beta, 0.11) << endl;

    cout << "=======================================================" << endl;
    cout << "Fits to primary SN sample with extinction correction\n";
    cout << "Omega_m: " << lcosmos.omega_m() << " Omega_L: " <<
      lcosmos.omega_v() << " Alpha: " << alpha  << " Beta: " << beta << endl;
    cout << "Primary Knop chi^2, extcorr: " << 
      ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::K03Primary, beta, 0.11) << 
      endl;

    lcosmos.setOmega_m(1.6);
    lcosmos.setOmega_v(0.3);
    lcosmos.reset(); lcosmos.history();
    cout << "Omega_m: " << lcosmos.omega_m() << " Omega_L: " <<
      lcosmos.omega_v() << " Alpha: " << alpha  << " Beta: " << beta << endl;
    cout << "Primary Knop chi^2, extcorr: " << 
      ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::K03Primary, beta, 0.11) << 
      endl;

    lcosmos.setOmega_m(0.7);
    lcosmos.setOmega_v(1.3);
    lcosmos.reset(); lcosmos.history();
    cout << "Omega_m: " << lcosmos.omega_m() << " Omega_L: " <<
      lcosmos.omega_v() << " Alpha: " << alpha  << " Beta: " << beta << endl;
    cout << "Primary Knop chi^2, extcorr: " << 
      ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::K03Primary, beta, 0.11) << 
      endl;

    lcosmos.setOmega_m(1.4);
    lcosmos.setOmega_v(2.4);
    lcosmos.reset(); lcosmos.history();
    cout << "Omega_m: " << lcosmos.omega_m() << " Omega_L: " <<
      lcosmos.omega_v() << " Alpha: " << alpha  << " Beta: " << beta << endl;
    cout << "Primary Knop chi^2, extcorr: " << 
      ai.SNIaKnop03( lcosmos, alpha, AnalyzeThis::K03Primary, beta, 0.11) << 
      endl;

  } catch (Bad_Error x) {
    cout << "leandriver bad error:" << endl;
    cout << x.s << endl;
  }	      
  return 0;
}
