/* mass_function.h */
#ifndef massfunction_h
#define massfunction_h
#include "spline.h"
#include "anchor.h"

class MassFunction {
	private:
	Anchor SplineAnchor;
	Spline *Pk;
	Spline *Sigma;
	Spline *dLnSigma_dLnM;
	Spline *LnSigma_LnM;
	Spline *dLnSigma_M;
	
	double m_max, m_min;
	double rho_0;
	double z;

	bool normalize;

	protected: 
	Spline *Integral;
	Spline *Differential;

// Set Pk and z in the constructor
	public: 
        ~MassFunction(){};
	MassFunction(Spline* pk, double min, double max, double r);
	Anchor* anchor(){return &SplineAnchor;}

	double norm(double M);
	double sigma8();
	double M_r(double);
	double R_m(double);
	double W_k(double,double);

	virtual double param(int, double){};
	virtual double differential_number_density(double M){};
	virtual double differential_number_density(double M, double z){};

	double integral_number_density(double M);
	double integral_number_density(double z, double M);

	double get_integral(double m);

	void set_integral();
	void set_differential();	
	virtual void set_differential(double z){};	

	void print_differential(string name){Differential->dump("differential_"+name);}
	void print_integral(string name){Integral->dump("integral_"+name);}
	virtual void init_standard_param(){};
	
	void set_pk(Spline *p){Pk=p;}
	void set_rho0(double r){rho_0=r;}
	void set_z(double a){z=a;}
	
	void init(double norm);
	void initSplines();
	void fill();
	void integrate_sigma();

	void normalize_sigma8(double);
	void set_normalize(bool n){normalize=n;}
	bool normalization(){return normalize;}
	void set_m_max(double m){m_max=m;}
	void set_m_min(double m){m_min=m;}
	
	double sigma(double M){return Sigma->fastY(M);}
	double M_max(){return m_max;}
	double M_min(){return m_min;}
};

class Tinker : public MassFunction
{
	protected:
	// Tinker mass function parameters
	double p[4];

	public: 
        ~Tinker(){};
	// Set Pk and z in the constructor
	Tinker(Spline* pk, double min, double max, double r);

	virtual double param(int i, double z){};	
	double differential_number_density(double M);
	
	// Initialize parameters with the standard values
	void init_standard_param(); 
};

class TinkerZ : public Tinker
{
	// Set Pk and z in the constructor
	public: 
        ~TinkerZ(){};
	TinkerZ(Spline* pk, double min, double max, double r);

	double differential_number_density(double z, double M);
	void set_differential(double z);	
	double param(int i, double z);	
};

class ShethTormen : public MassFunction
{
	protected:
	// ShethTormen mass function parameters
	double p[3];

	// Set Pk and z in the constructor
	public: 
        ~ShethTormen(){};
	ShethTormen(Spline* pk, double min, double max, double r);

	virtual double param(int i, double z){};	
	double differential_number_density(double M);

	// Initialize parameters with the standard values
	void init_standard_param(); 
};
#endif
