#include "helper.h"
#include "vector.h"
#include "numath.h"
#include "osc_defs_NH.h"
#include <complex>

//	calculate e^(i * d_cp) CP violation phase
std::complex<double> phase = std::exp(std::complex<double>(1.0i * d_cp));

//	complex MNS matrix
complex<double> MNS[] = 
{      c12 * c13, 							 s12 * c13, 			  		      s13 * phase,
	  -s12 * c23 - c12 * s23 * s13 * phase,  c12 * c23 - s12 * s23 * s13 * phase, s23 * c13,
	   s12 * s23 - c12 * s23 * s13 * phase, -c12 * s23 - s12 * c23 * s13 * phase, c23 * c13
};


//	probability of a neutrino initially in flavor a to transition to flavor b, as a function
//	of L/E (km/GeV)
double P(flavor a, flavor b, double x) {
	double p = 0.;

	//	formula found in Zuber p.192
	for(int i=0; i<3; ++i) {
		p += std::pow(std::abs(MNS[a*3 + i] * std::conj(MNS[b*3 + i])), 2);
	}

	for(int i=0; i<3; ++i) {
		for(int j=0; j<3; ++j) {
			if(j > i) {
				p += 2. * std::real(MNS[a*3 + i] * std::conj(MNS[a*3 + j]) 
					   			  * std::conj(MNS[b*3 + i]) * MNS[b*3 + j] 
					    		  * std::exp(std::complex<double>(-2.i * 1.2668 * dm2_mat[i*3 + j] * x)));
			}
		}
	}

	return p;
}

double plot_P(double* x, double* par) {
	return P((flavor)par[0], (flavor)par[1], x[0]);
}

//	main: in_f is the flavor at t=0
void nucp(int in_f = 0) {
	TCanvas* c1 = new TCanvas();
	if(in_f == f_e) {
		TF1* e_e = new TF1("", plot_P, 0., 35000., 2);
		e_e->SetParameters(f_e, f_e);
		TF1* e_m = new TF1("", plot_P, 0., 35000., 2);
		e_m->SetParameters(f_e, f_m);
		TF1* e_t = new TF1("", plot_P, 0., 35000., 2);
		e_t->SetParameters(f_e, f_t);
		
		e_e->SetRange(0., 35000.);
		e_e->SetNpx(10000);
		e_e->SetMinimum(0.0);
		e_e->SetMaximum(1.0);
		e_m->SetNpx(10000);
		e_t->SetNpx(10000);
		
		e_m->SetLineColor(4);
		e_t->SetLineColor(8);

		e_e->Draw();
		e_m->Draw("same");
		e_t->Draw("same");
	}
	if(in_f == f_m) {
		TF1* m_e = new TF1("", plot_P, 0., 35000., 2);
		m_e->SetParameters(f_m, f_e);
		TF1* m_m = new TF1("", plot_P, 0., 35000., 2);
		m_m->SetParameters(f_m, f_m);
		TF1* m_t = new TF1("", plot_P, 0., 35000., 2);
		m_t->SetParameters(f_m, f_t);
		
		m_e->SetRange(0., 35000.);
		m_e->SetNpx(10000);
		m_e->SetMinimum(0.0);
		m_e->SetMaximum(1.0);
		m_m->SetNpx(10000);
		m_t->SetNpx(10000);
		
		m_m->SetLineColor(4);
		m_t->SetLineColor(8);

		m_e->Draw();
		m_m->Draw("same");
		m_t->Draw("same");
	}
	if(in_f == f_t) {
		TF1* t_e = new TF1("", plot_P, 0., 35000., 2);
		t_e->SetParameters(f_t, f_e);
		TF1* t_m = new TF1("", plot_P, 0., 35000., 2);
		t_m->SetParameters(f_t, f_m);
		TF1* t_t = new TF1("", plot_P, 0., 35000., 2);
		t_t->SetParameters(f_t, f_t);
		
		t_e->SetRange(0., 35000.);
		t_e->SetNpx(10000);
		t_e->SetMinimum(0.0);
		t_e->SetMaximum(1.0);
		t_m->SetNpx(10000);
		t_t->SetNpx(10000);
		
		t_m->SetLineColor(4);
		t_t->SetLineColor(8);

		t_e->Draw();
		t_m->Draw("same");
		t_t->Draw("same");
	}
}
