#ifndef MASSIVENEUTRINOS_H
#define MASSIVENEUTRINOS_H

#define NRHOPN 10000 //! From CMBFAST

namespace MassiveNeutrinos
{
   //! propagate massive neutrino moments in longitudinal gauge
   void propagateMassiveNeutrinoMomentsLongitudinal(const double *y, double *yprime,
                            const double a, const double tau,
                            const double nuMass, const double k,
                            const double PhiDot, const double Psi,
                            const double beta, const double X);
   void propagateMassiveNeutrinoMomentsSynchronous(const double *y, double *yprime,
                            const double a, const double tau,
                            const double nuMass, const double k,
                            const double hdot, const double etadot);
   void nu2(const double a, const double betanu,
            const double rhonu, const double pnu, const double deltaphi,
            double *drhonu, double *fnu, double *dpnu, double *shearnu,
            const double *psi0, const double * psi1, const double *psi2);
    double nuder(const double a, const double adotoa, const double beta,
                 const double rhonu, const double phidot,
                 const double *psi2, const double *psi2dot);

    void init();

    void setFirstIndex(unsigned int i);
    unsigned int qGridSize();
    unsigned int pertArraySize();

    void ninu1(const double a, double *rhonu, double *pnu,  const double amnu);
    void nu1(const double a, double *rhonu, double *pnu);


    struct Nu1d_ {
      double amin, amax,  dlna, r1[NRHOPN], p1[NRHOPN], dr1[NRHOPN], dp1[NRHOPN], ddr1[NRHOPN];
    };

    void initnu1(const double amnu);

    void setIninitalLongitudinalScalarPerturbations(double *y, const double a,
                                        const double mass_nu, const double deltan,
                                        const double Vn, const double Pi_nu);
}

#endif //MASSIVENEUTRINOS_H
