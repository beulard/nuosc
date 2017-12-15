#include "helper.h"
#include "vector.h"
#include "numath.h"
#include "osc_defs.h"
#include <complex>

hierarchy* h = new hierarchy;

//	probability of a neutrino initially in flavor a to transition to flavor b, as a function
//	of L/E (km/GeV) expression from zuber p. 192
double P(flavor a, flavor b, double x) {
	double p = 0.;

	//	formula found in Zuber p.192
	for(int i=0; i<3; ++i) {
		p += pow(abs(h->MNS[a*3 + i] * std::conj(h->MNS[b*3 + i])), 2);
	}

	for(int i=0; i<3; ++i) {
		for(int j=0; j<3; ++j) {
			if(j > i) {
				p += 2. * real(h->MNS[a*3 + i] * conj(h->MNS[a*3 + j]) 
				  * conj(h->MNS[b*3 + i]) * h->MNS[b*3 + j] 
				  * exp(complex<double>(-2.i * 1.2668 * h->dm2_mat[i*3 + j] * x)));
			}
		}
	}

	return p;
}

// alternative expression for P (zuber p. 196)
/*double P(flavor a, flavor b, double x) {
	double p = (a == b ? 1.0 : 0.0);

	for (int i=0; i<3; ++i) {
		for (int j=0; j<3; ++j) {
			if (i > j) {
				complex<double> K = MNS[a*3 + i] * conj(MNS[b*3 + i]) 
					  * conj(MNS[a*3 + j]) * MNS[b*3 + j];
					
				p -= 4.0 * real(K) * sinsq(1.2668 * dm2_mat[i*3 + j] * x);
				p += 4.0 * imag(K) * TMath::Sin(1.2668 * dm2_mat[i*3 + j] * x)
					     * TMath::Cos(1.2668 * dm2_mat[i*3 + j] * x);
			}
		}
	}
	return p;
}*/

double plot_P(double* x, double* par) {
	return P((flavor)par[0], (flavor)par[1], x[0]);
}

//	main: in_f is the flavor at t=0
void nucp(int in_h = NH, int in_f = 0) {
	populate(h, (h_type)in_h);

	// Test if the probability adds up to 1 at every point.
	// It does!
	/*int N = 35000;
	for (int i=0; i<N; ++i) {
		// we want L to go from 15000 to 25000
		float L = 0 + (float)i / N * 35000.;
		float totP = P(f_e, f_e, L) + P(f_e, f_m, L) + P(f_e, f_t, L);
		Printf("%f", totP);
	}*/
	
	TCanvas* c1 = new TCanvas();
	c1->SetFillColor(ci[CI_BACKGROUND]);
	if(in_f == f_e) {
		TF1* e_e = new TF1("#nu_{e} survival", plot_P, 0., 35000., 2);
		e_e->SetParameters(f_e, f_e);
		TF1* e_m = new TF1("#nu_{#mu} appearance", plot_P, 0., 35000., 2);
		e_m->SetParameters(f_e, f_m);
		TF1* e_t = new TF1("#nu_{#tau} appearance", plot_P, 0., 35000., 2);
			
		e_t->SetParameters(f_e, f_t);
		
		e_e->SetRange(0., 35000.);
		e_e->SetNpx(10000);
		e_e->SetMinimum(0.0);
		e_e->SetMaximum(1.0);
		e_m->SetNpx(10000);
		e_t->SetNpx(10000);
		
		e_e->SetLineColor(ci[CI_E]);
		e_m->SetLineColor(ci[CI_MU]);
		e_t->SetLineColor(ci[CI_TAU]);
	//		e_e->SetRange(0., 5000.);
		e_e->Draw();
		e_m->Draw("same");
		e_t->Draw("same");
		c1->BuildLegend();
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
		
		m_e->SetLineColor(ci[CI_E]);
		m_m->SetLineColor(ci[CI_MU]);
		m_t->SetLineColor(ci[CI_TAU]);

		m_e->SetTitle("#nu_{#mu} transition");
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
		
		t_e->SetLineColor(ci[CI_E]);
		t_m->SetLineColor(ci[CI_MU]);
		t_t->SetLineColor(ci[CI_TAU]);

		t_e->SetTitle("#nu_{#tau} transition");
		t_e->Draw();
		t_m->Draw("same");
		t_t->Draw("same");
	}
}
