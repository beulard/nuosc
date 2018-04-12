#include "osc_defs.h"
#include "TMath.h"
#include "numath.h"
#include "TString.h"


void parameters::populate_nh(double d) {
	//h->dm2_21 = 7.37e-5;
	
	//h->dm2_21 = 7.5e-5;
	dm2_21 = pdg[NH].dm2_21;
	dm2_31 = pdg[NH].dm2_31;
	dm2_32 = pdg[NH].dm2_32;

	t12 = pdg[NH].t12;
	t13 = pdg[NH].t13;
	t23 = pdg[NH].t23;
	//h->dm2 = 2.50e-3;
	//h->dm2_31 = h->dm2 + h->dm2_21 / 2.;
	//h->dm2_32 = h->dm2 - h->dm2_21 / 2.;
	
	//h->dm2_31 = 2.457e-3;
	//h->dm2_32 = h->dm2_31 - h->dm2_21;

	//h->t12 = TMath::ASin(TMath::Sqrt(0.297));
	//h->t23 = TMath::ASin(TMath::Sqrt(0.437));
	//h->t13 = TMath::ASin(TMath::Sqrt(0.0214));
	//h->t12 = 0.5843;
	//h->t23 = 0.738;
	//h->t13 = 0.148;

	// Default to best fit value
	if (d < -1000.) {
		d_cp = pdg[NH].d_cp;
	} else {
		d_cp = d;
	}
}

void parameters::populate_ih(double d) {
	//h->dm2_21 = 7.37e-5;
	//h->dm2_21 = 7.5e-5;
	//h->dm2_31 = -2.449e-3;
	dm2_21 = pdg[IH].dm2_21;
	dm2_31 = pdg[IH].dm2_31;
	dm2_32 = pdg[IH].dm2_32;

	t12 = pdg[IH].t12;
	t13 = pdg[IH].t13;
	t23 = pdg[IH].t23;
	

	//h->dm2 = -2.46e-3;
	//h->dm2_31 = h->dm2 + h->dm2_21 / 2.;
	//h->dm2_32 = h->dm2 - h->dm2_21 / 2.;
	//h->dm2_32 = h->dm2_31 - h->dm2_21;

	//h->t12 = TMath::ASin(TMath::Sqrt(0.297));
	//h->t23 = TMath::ASin(TMath::Sqrt(0.569));
	//h->t13 = TMath::ASin(TMath::Sqrt(0.0218));
	//h->t12= 0.5843;
	//h->t23 = 0.864;
	//h->t13 = 0.148;
	
	// Default to best fit value
	if (d < -1000.) {
		d_cp = pdg[IH].d_cp;
	} else {
		d_cp = d;
	}
}

//	populates only calculated values
void parameters::populate_common() {
	//	populate dm2 matrix
	double tmp_dm2_mat[] = { 0.0,      -dm2_21, -dm2_31,
			     dm2_21, 0.0,       -dm2_32,
			     dm2_31, dm2_32,  0.0  };

	memcpy(dm2_mat, tmp_dm2_mat, 9 * sizeof(double));

	s12 = TMath::Sin(t12);
	s23 = TMath::Sin(t23);
	s13 = TMath::Sin(t13);

	c12 = TMath::Cos(t12);
	c23 = TMath::Cos(t23);
	c13 = TMath::Cos(t13);

	id_cp = complex<double>(0.0, d_cp);

	complex<double> tmp_MNS_23[] = 
	{		1.,  0.,  0.,
			0.,  c23, s23,
			0., -s23, c23		};
	complex<double> tmp_MNS_13[] = 
	{		c13, 			   0., s13 * exp(-id_cp),
			0., 			   1., 0.,
			-s13 * exp(id_cp), 0., c13		};
	complex<double> tmp_MNS_12[] =
	{		 c12, s12, 0.,
			-s12, c12, 0.,
			0.,   0.,  1.		};
	
	memcpy(MNS_23, tmp_MNS_23, 9 * sizeof(complex<double>));
	memcpy(MNS_13, tmp_MNS_13, 9 * sizeof(complex<double>));
	memcpy(MNS_12, tmp_MNS_12, 9 * sizeof(complex<double>));
	
	complex<double> temp[9];
	//complex<double> MNS[9];
	
	mat_mult(MNS_23, MNS_13, temp);
	mat_mult(temp, MNS_12, MNS);

	//memcpy(h->MNS, MNS, 9 * sizeof(complex<double>));
}

void parameters::populate(h_type t, float d_cp) {
	type = t;
	if(t == NH)
		populate_nh(d_cp);
	else if(t == IH)
		populate_ih(d_cp);
	else {
		Printf("Invalid hierarchy type");
		return;
	}
	populate_common();
}

void parameters::flip_hierarchy() {
	dm2_31 = -dm2_31;
	dm2_32 = dm2_31 - dm2_21;
	// flip the type
	type = (h_type)(1 - type);

	// And call populate_common to repopulate dm2 matrix 
	populate_common();
}

