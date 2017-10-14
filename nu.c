#include "helper.h"
#include "vector.h"
#include "numath.h"
#include "osc_defs_NH.h"

//	(real) MNS matrix
const double MNS[] = {  c12 * c13, 					  s12 * c13, 				   s13,
					   -s12 * c23 - c12 * s23 * s13,  c12 * c23 - s12 * s23 * s13, s23 * c13,
					    s12 * s23 - c12 * s23 * s13, -c12 * s23 - s12 * c23 * s13, c23 * c13 };

//	probability of a neutrino initially in flavor a to transition to flavor b, as a function
//	of L/E (km/GeV)
double P(flavor a, flavor b, double x) {
	double delta = a == b ? 1.0 : 0.0;
	double p = delta;
	
	//	formula found in Zuber p.192
	for(int i=0; i < 3; ++i) {
		for(int j=0; j < 3; ++j) {
			if(i > j) {
				p -= 4.0 * MNS[a*3 + i] * MNS[b*3 + i] * MNS[a*3 + j] * MNS[b*3 + j] 
						 * sinsq(1.2668 * dm2_mat[i*3 + j] * x);
			}
		}
	}
	
	return p;
}

double plot_P(double* x, double* par) {
	return P((flavor)par[0], (flavor)par[1], x[0]);
}


// main
void nu(int in_f = 0) {
	TCanvas* c1 = new TCanvas();
	if(in_f == 0) {
		TF1* e_e = new TF1("e_e", plot_P, 0., 35000., 2);
		e_e->SetParameters(f_e, f_e);
		TF1* e_m = new TF1("e_m", plot_P, 0., 35000., 2);
		e_m->SetParameters(f_e, f_m);
		TF1* e_t = new TF1("e_t", plot_P, 0., 35000, 2);
		e_t->SetParameters(f_e, f_t);


		e_e->SetNpx(10000);
		e_e->SetMinimum(0.0);
		e_e->SetMaximum(1.0);
		// TODO PLOT LEGEND?
		e_e->SetTitle("#nu_{e} transition");
		e_t->SetNpx(10000);
		e_m->SetNpx(10000);

		e_t->SetLineColor(8);
		e_m->SetLineColor(4);

		e_e->GetXaxis()->SetTitle("L/E (km/GeV)");
		e_e->GetYaxis()->SetTitle("");

		e_e->Draw();
		e_m->Draw("same");
		e_t->Draw("same");

	}
	else if(in_f == 1) {
		TF1* m_e = new TF1("m_e", plot_P, 0., 35000., 2);
		m_e->SetParameters(f_m, f_e);
		TF1* m_m = new TF1("m_m", plot_P, 0., 35000., 2);
		m_m->SetParameters(f_m, f_m);
		TF1* m_t = new TF1("m_t", plot_P, 0., 35000, 2);
		m_t->SetParameters(f_m, f_t);
		m_e->SetNpx(10000);
		m_e->SetMinimum(0.0);
		m_e->SetMaximum(1.0);
		m_t->SetNpx(10000);
		m_m->SetNpx(10000);

		m_t->SetLineColor(8);
		m_m->SetLineColor(4);
	
		m_e->GetXaxis()->SetTitle("L/E (km/GeV)");
		m_e->GetYaxis()->SetTitle("");

		m_e->SetTitle("#nu_{#mu} transition");
		m_e->Draw();
		m_m->Draw("same");
		m_t->Draw("same");
		

	}
	else if(in_f == 2) {
		TF1* t_e = new TF1("t_e", plot_P, 0., 35000., 2);
		t_e->SetParameters(f_t, f_e);
		TF1* t_m = new TF1("t_m", plot_P, 0., 35000., 2);
		t_m->SetParameters(f_t, f_m);
		TF1* t_t = new TF1("m_t", plot_P, 0., 35000, 2);
		t_t->SetParameters(f_t, f_t);
		t_e->SetNpx(10000);
		t_e->SetMinimum(0.0);
		t_e->SetMaximum(1.0);
		t_t->SetNpx(10000);
		t_m->SetNpx(10000);

		t_t->SetLineColor(8);
		t_m->SetLineColor(4);

		t_e->GetXaxis()->SetTitle("L/E (km/GeV)");
		t_e->GetYaxis()->SetTitle("");

		t_e->SetTitle("#nu_{#tau} transition");
		t_e->Draw();
		t_m->Draw("same");
		t_t->Draw("same");
	}
}

