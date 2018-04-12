#include "helper.h"
#include "vector.h"
#include "numath.h"
#include "osc_defs.h"
#include <complex>
#include "spectrum.h"

parameters* p = new parameters;


double plot_P(double* x, double* par) {
	return P((flavor)par[0], (flavor)par[1], 1., x[0], p, false);
}

//	main: in_f is the flavor at t=0
void nucp(int in_h = NH, int in_f = 0) {
	p->populate(NH, 0);
	if(in_h == 1) {
		p->flip_hierarchy();
	}

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

		e_e->GetYaxis()->SetTitle("P(L, E)");
		e_e->GetXaxis()->SetTitle("L/E (km/GeV)");

		e_e->GetXaxis()->SetTitleSize(0.049);
		e_e->GetXaxis()->SetLabelSize(0.045);
		e_e->GetYaxis()->SetTitleSize(0.049);
		e_e->GetYaxis()->SetLabelSize(0.045);
		e_e->GetYaxis()->SetTitleOffset(0.9);

		//e_e->SetLineWidth(2);
		//e_m->SetLineWidth(2);
		//e_t->SetLineWidth(2);
	//		e_e->SetRange(0., 5000.);
		e_e->Draw();
		e_m->Draw("same");
		e_t->Draw("same");
		c1->BuildLegend();
		e_e->SetTitle("");
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
		c1->BuildLegend();
		m_e->SetTitle("");
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
