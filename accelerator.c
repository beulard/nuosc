#include "osc_defs_NH.h"

// approximation functions
double e_transition(flavor final, double x) {
	switch(final) {
		case f_m:
			return sinsq(2.0 * t13) * sinsq(t23) * sinsq(1.2668 * dm2_32 * x);
			break;
		case f_t:
			return sinsq(2.0 * t13) * cossq(t13) * sinsq(1.2668 * dm2_32 * x);
			break;
		default:
			return 0.0;
			break;
	}
}

double mu_transition(flavor final, double x) {
	switch(final) {
		case f_t:
			return sinsq(2.0 * t23) * cossq(t13) * cossq(t13) * sinsq(1.2668 * dm2_32 * x);
			break;
		case f_e:
			return sinsq(2.0 * t13) * sinsq(t23) * sinsq(1.2668 * dm2_32 * x);
			break;
		default:
			return 0.0;
			break;
	}
}

//	plotting functions
double plot_transition(double* x, double* par) {
	if(par[0] == 0) 
		return e_transition((flavor)par[1], x[0]);
	else
		return mu_transition((flavor)par[1], x[0]);
}


void accelerator() {

		TCanvas* c1 = new TCanvas("c1", "", 900, 500);
		c1->Divide(2, 1);
		c1->SetTitle("accelerator experiment");
		// plot electron approximation
		TF1* e_m2 = new TF1("", plot_transition, 0., 5000., 2);
		e_m2->SetParameters(0, f_m);
		TF1* e_t2 = new TF1("", plot_transition, 0., 5000., 2);
		e_t2->SetParameters(0, f_t);
		e_t2->SetMaximum(.1);
		e_t2->SetMinimum(0.0);
		e_t2->SetNpx(10000);
		e_m2->SetNpx(10000);
		e_t2->SetLineColor(8);
		e_m2->SetLineColor(4);
		e_t2->SetTitle("e transition");
		c1->cd(1);
		e_t2->GetXaxis()->SetTitle("L/E (km/GeV)");
		e_t2->GetYaxis()->SetTitle("");
		e_t2->Draw();
		e_m2->Draw("same");

		//	plot muon approx
		c1->cd(2);
		
		TF1* m_e2 = new TF1("", plot_transition, 0., 5000., 2);
		m_e2->SetParameters(0, f_e);
		TF1* m_t2 = new TF1("", plot_transition, 0., 5000., 2);
		m_t2->SetParameters(0, f_t);
		m_t2->SetMaximum(.1);
		m_t2->SetMinimum(0.0);
		m_t2->SetNpx(10000);
		m_e2->SetNpx(10000);
		m_t2->SetLineColor(8);
		m_t2->SetTitle("#mu transition");
		m_t2->GetXaxis()->SetTitle("L/E (km/GeV)");
		m_t2->GetYaxis()->SetTitle("");
		m_t2->Draw();
		m_e2->Draw("same");
}
