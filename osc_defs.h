#ifndef OSC_DEFS
#define OSC_DEFS
#include <complex>

enum flavor {
	f_e = 0,
	f_m = 1,
	f_t = 2
};

// Hierarchy type
enum h_type {
	NH = 0,
	IH
};


// Describes a set of oscillation parameters, used to calculate oscillation probabilities
class parameters {
public:
	// Base hierarchy type (NH, IH)
	h_type type;
	////	Numerical data
	//	mass differences, in eV^2
	double dm2_21;
	double dm2_31;
	double dm2_32;
	
	//	mixing angles
	double t12;
	double t23;
	double t13;
	
	double d_cp;

	////	Calculated values	
	//	mass squared differences in a matrix form
	double dm2_mat[9];
	//	trigonometric
	double s12, s23, s13;
	double c12, c23, c13;

	complex<double> MNS[9];
	//	CP violating phase times i
	complex<double> id_cp;
	
	//	MNS sub-matrices
	complex<double> MNS_23[9];
	complex<double> MNS_13[9];
	complex<double> MNS_12[9];


	// Populate the data with a given hierarchy and best fit parameters (defined in pdg)
	void populate(h_type t, float d_cp = -1111.);
	//	does not change numerical data but calculates derived quantities (trigonometrics, MNS) 
	void populate_common();
	// Change the mass hierarchy
	void flip_hierarchy();

private:
	// Populate with PDG values
	void populate_nh(double d_cp);
	void populate_ih(double d_cp);
};

// Global PDG best fit values and uncertainties
parameters pdg[2];
parameters pdg_sd[2];


#endif
