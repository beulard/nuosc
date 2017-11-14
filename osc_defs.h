#ifndef OSC_DEFS
#define OSC_DEFS
#include <complex>

enum flavor {
	f_e = 0,
	f_m = 1,
	f_t = 2
};

enum h_type {
	NH,
	IH
};

struct hierarchy {
	// Base hierarchy type (NH, IH)
	h_type type;
	////	Numerical data
	//	mass differences, in eV^2
	double dm2_21;
	//	absolute value of dm^2 as given in pdg neutrino mixing paper
	double dm2;
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
};


//	helper function to populate hierarchy object
void populate(hierarchy* h, h_type t, float d_cp = -1111.);
//	does not change numerical data but calculates derived quantities (trigonometrics, MNS) 
void populate_common(hierarchy* h);

#endif
