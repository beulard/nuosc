#include "osc_defs_common.h"

//	baseline for experiment: DUNE has L = 1300km
const double L = 1300.0;


double P(double E) {
	double p = 0.;

	//	initial and final states
	const int a = f_m;
	const int b = f_e;

	for(int i=0; i<3; ++i) {
		p += pow(abs(MNS[a*3 + i] * std::conj(MNS[b*3 + i])), 2);
	}

	for(int i=0; i<3; ++i) {
		for(int j=0; j<3; ++j) {
			if(j > i) {
				p += 2. * std::real(MNS[a*3 + i] * std::conj(MNS[a*3 + j]) 
						* std::conj(MNS[b*3 + i]) * MNS[b*3 + j]
						* std::exp(complex<double>(-2.i * 1.2668 *  dm2_mat[i*3 + j] * L/E)));
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

	//	TODO restructure code to allow for comparison between
	//	TODO		NH and IH
	//	TODO	values of delta CP
	TCanvas* c1 = new TCanvas();
	TF1* m_e = new TF1("", plot_P, 0.5, 10.0);

	m_e->SetMaximum(.08);
	m_e->SetMinimum(0.0);
	m_e->SetNpx(10000);
	m_e->Draw();
}
