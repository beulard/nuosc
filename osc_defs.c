#include "osc_defs.h"
#include "TMath.h"
#include "numath.h"
#include "TString.h"

void populate_nh(hierarchy* h, float d_cp) {
	h->dm2_21 = 7.37e-5;
	h->dm2 = 2.50e-3;
	h->dm2_31 = h->dm2 + h->dm2_21 / 2.;
	h->dm2_32 = h->dm2 - h->dm2_21 / 2.;

	h->t12 = TMath::ASin(TMath::Sqrt(0.297));
	h->t23 = TMath::ASin(TMath::Sqrt(0.437));
	h->t13 = TMath::ASin(TMath::Sqrt(0.0214));

	// Default to best fit value
	if (d_cp < 0.) {
		h->d_cp = TMath::Pi() * 1.35;
	} else {
		h->d_cp = d_cp;
	}
}

void populate_ih(hierarchy* h, float d_cp) {
	h->dm2_21 = 7.37e-5;
	h->dm2 = -2.46e-3;
	h->dm2_31 = h->dm2 + h->dm2_21 / 2.;
	h->dm2_32 = h->dm2 - h->dm2_21 / 2.;

	h->t12 = TMath::ASin(TMath::Sqrt(0.297));
	h->t23 = TMath::ASin(TMath::Sqrt(0.569));
	h->t13 = TMath::ASin(TMath::Sqrt(0.0218));
	
	// Default to best fit value
	if (d_cp < 0.) {
		h->d_cp = TMath::Pi() * 1.32;
	} else {
		h->d_cp = d_cp;
	}
}

//	populates calculated values
void populate_common(hierarchy* h) {
	//	populate dm2 matrix
	double dm2_mat[] = { 0.0,      -h->dm2_21, -h->dm2_31,
			     h->dm2_21, 0.0,       -h->dm2_32,
			     h->dm2_31, h->dm2_32,  0.0  };

	memcpy(h->dm2_mat, dm2_mat, 9 * sizeof(double));

	h->s12 = TMath::Sin(h->t12);
	h->s23 = TMath::Sin(h->t23);
	h->s13 = TMath::Sin(h->t13);

	h->c12 = TMath::Cos(h->t12);
	h->c23 = TMath::Cos(h->t23);
	h->c13 = TMath::Cos(h->t13);

	h->id_cp = complex<double>(0.0, h->d_cp);

	complex<double> MNS_23[] = 
	{		1.,  0.,  0.,
			0.,  h->c23, h->s23,
			0., -h->s23, h->c23		};
	complex<double> MNS_13[] = 
	{		h->c13, 			   0., h->s13 * exp(-h->id_cp),
			0., 			   1., 0.,
			-h->s13 * exp(h->id_cp), 0., h->c13		};
	complex<double> MNS_12[] =
	{		 h->c12, h->s12, 0.,
			-h->s12, h->c12, 0.,
			0.,   0.,  1.		};
	
	memcpy(h->MNS_23, MNS_23, 9 * sizeof(complex<double>));
	memcpy(h->MNS_13, MNS_13, 9 * sizeof(complex<double>));
	memcpy(h->MNS_12, MNS_12, 9 * sizeof(complex<double>));
	
	complex<double> temp[9];
	complex<double> MNS[9];
	
	mat_mult(MNS_23, MNS_13, temp);
	mat_mult(temp, MNS_12, MNS);
	memcpy(h->MNS, MNS, 9 * sizeof(complex<double>));
}

void populate(hierarchy* h, h_type t, float d_cp) {
	h->type = t;
	if(t == NH)
		populate_nh(h, d_cp);
	else if(t == IH)
		populate_ih(h, d_cp);
	else
		return;
	populate_common(h);
}
