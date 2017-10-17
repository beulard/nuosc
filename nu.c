#include "osc_defs_common.h"

//	baseline for experiment: DUNE has L = 1300km
const double L = 1300.0;

double P(double E) {
	double p = 0.;
	//	initial and final states
	int a = f_m;
	int b = f_e;

	for(int i=0; i<3; ++i) {
		p += pow(abs(MNS[a*3 + i] * std::conj(MNS[b*3 + i])), 2);		
	}

	for(int i=0; i<3; ++i) {
		for(int j=0; j<3; ++j) {
			if(j > i) {
				complex<double> x = -1.i * 2 * 1.2668 * dm2_mat[i*3 + j] * 500 / E;
				p += 2. * real(MNS[a*3 + i] * conj(MNS[a*3 + j]) * conj(MNS[b*3 + i]) * MNS[b*3 + j]
				    					    * exp(x));
			}
		}
	}
	return p;
}

double plot_P(double* x, double* par) {
	return P(x[0]);
}

void nu() {
	//	calculate MNS matrix by multiplying rotations
	complex<double> temp[9];
	mat_mult(MNS_23, MNS_13, temp);
	mat_mult(temp, MNS_12, MNS);


	// TODO WHATS GOING ON?
	TCanvas* c1 = new TCanvas();
	TF1* m_e = new TF1("", "P(x)", 0.0, 5.0);
	m_e->SetMaximum(1.0);
	m_e->SetMinimum(0.0);
	m_e->SetNpx(500000);
	m_e->Draw();
}
