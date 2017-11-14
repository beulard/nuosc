// change this line to change hierarchy everywhere
#include "osc_defs_IH.h"


// mass difference matrix, for indexing
const double dm2_mat[] = { 0.0,    -dm2_21, -dm2_31,
			   dm2_21,  0.0,    -dm2_32,
			   dm2_31,  dm2_32,  0.0  };

// trigonometric functions
const double s12 = TMath::Sin(t12);
const double s23 = TMath::Sin(t23);
const double s13 = TMath::Sin(t13);

const double c12 = TMath::Cos(t12);
const double c23 = TMath::Cos(t23);
const double c13 = TMath::Cos(t13);

//	already defined in osc_defs.h which is read when root starts
/*enum flavor {
	f_e = 0,
	f_m = 1,
	f_t = 2
};*/

//	strongly typed i * d_cp to avoid ambiguity when using exp
const complex<double> id_cp = 1.i * d_cp;

// global MNS matrix, defined in the main function nucp() 
complex<double> MNS[9];

//	rotation matrices that multiply (in this order) to give the MNS mixing matrix
const complex<double> MNS_23[] = 
{		1.,  0.,  0.,
		0.,  c23, s23,
		0., -s23, c23		};

const complex<double> MNS_13[] =
{		c13, 			   0., s13 * exp(-id_cp),
		0., 			   1., 0.,
		-s13 * exp(id_cp), 0., c13		};

const complex<double> MNS_12[] =
{		 c12, s12, 0.,
		-s12, c12, 0.,
		0.,   0.,  1.		};


