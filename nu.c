#include "helper.h"
#include "vector.h"
#include "numath.h"

using namespace TMath;

//	mass differences, in eV^2
const double dm2_21 = 7.37e-5;
//	absolute value of dm^2 as given in pdg neutrino mixing paper
const double dm2 = 2.50e-3;
const double dm2_31 = dm2 + dm2_21 / 2.;
const double dm2_32 = dm2 - dm2_21 / 2.;

// mass difference matrix, for indexing
const double dm2_mat[] = { 0.0,    dm2_21, dm2_31,
					   dm2_21, 0.0,    dm2_32,
					   dm2_31, dm2_32, 0.0  };

//	mixing angles
const double t12 = ASin(Sqrt(0.297));
const double t23 = ASin(Sqrt(0.437));
const double t13 = ASin(Sqrt(0.0214));

// trigonometric functions
const double s12 = Sin(t12);
const double s23 = Sin(t23);
const double s13 = Sin(t13);

const double c12 = Cos(t12);
const double c23 = Cos(t23);
const double c13 = Cos(t13);

enum flavor {
	f_e = 0,
	f_m = 1,
	f_t = 2
};


//	MNS matrix
const double MNS[] = {  c12 * c13, 					  s12 * c13, 				   s13,
					   -s12 * c23 - c12 * s23 * s13,  c12 * c23 - s12 * s23 * s13, s23 * c13,
					    s12 * s23 - c12 * s23 * s13, -c12 * s23 - s12 * c23 * s13, c23 * c13 };


//	probability of a neutrino initially in flavor a to transition to flavor b, as a function
//	of L/E (km/GeV)
double P(flavor a, flavor b, double x) {
	double delta = a == b ? 1.0 : 0.0;
	double p = delta;

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

double e_transition(flavor final, double x) {
	switch(final) {
		case f_t:
			return sinsq(2.0 * t13) * sinsq(t13) * sinsq(1.2668 * dm2_32 * x);
			break;
		case f_m:
			return sinsq(2.0 * t13) * cossq(t23) * sinsq(1.2668 * dm2_32 * x);
			break;
		default:
			return 0.0;
			break;
	}
}

double plot_transition(double* x, double* par) {
	if(par[0] == 0) 
		return e_transition((flavor)par[1], x[0]);
	else
		return mu_transition((flavor)par[1], x[0]);
}

double plotP(double* x, double* par) {
	return P((flavor)par[0], (flavor)par[1], x[0]);
}


// main
void nu(int in_f) {
	if(in_f == 0) {
		TCanvas* c1 = new TCanvas();
		TF1* e_e = new TF1("e_e", plotP, 0., 35000., 2);
		e_e->SetParameters(f_e, f_e);
		TF1* e_m = new TF1("e_m", plotP, 0., 35000., 2);
		e_m->SetParameters(f_e, f_m);
		TF1* e_t = new TF1("e_t", plotP, 0., 35000, 2);
		e_t->SetParameters(f_e, f_t);
		e_e->SetNpx(10000);
		e_e->SetMinimum(0.0);
		e_e->SetMaximum(1.0);
		e_e->SetTitle("electron neutrino transition");
		e_t->SetNpx(10000);
		e_m->SetNpx(10000);

		e_t->SetLineColor(8);
		e_m->SetLineColor(4);

		e_e->Draw();
		e_m->Draw("same");
		e_t->Draw("same");

		TCanvas* c2 = new TCanvas();
		// plot approximations
		TF1* e_m2 = new TF1("", plot_transition, 0., 5000., 2);
		e_m2->SetParameters(0, f_m);
		TF1* e_t2 = new TF1("", plot_transition, 0., 5000., 2);
		e_t2->SetParameters(0, f_t);
		e_t2->SetMaximum(.05);
		e_t2->SetMinimum(0.0);
		e_t2->SetNpx(10000);
		e_m2->SetNpx(10000);
		e_t2->SetLineColor(4);
		e_t2->Draw();
		e_m2->Draw("same");
	}
	else if(in_f == 1) {
		TF1* m_e = new TF1("m_e", plotP, 0., 35000., 2);
		m_e->SetParameters(f_m, f_e);
		TF1* m_m = new TF1("m_m", plotP, 0., 35000., 2);
		m_m->SetParameters(f_m, f_m);
		TF1* m_t = new TF1("m_t", plotP, 0., 35000, 2);
		m_t->SetParameters(f_m, f_t);
		m_e->SetNpx(10000);
		m_e->SetMinimum(0.0);
		m_e->SetMaximum(1.0);
		m_t->SetNpx(10000);
		m_m->SetNpx(10000);

		m_t->SetLineColor(8);
		m_m->SetLineColor(4);
	
		m_e->SetTitle("muon  neutrino transition");
		m_e->Draw();
		m_m->Draw("same");
		m_t->Draw("same");
		

		TCanvas* c2 = new TCanvas();
		// plot approximations
		TF1* m_e2 = new TF1("", plot_transition, 0., 5000., 2);
		m_e2->SetParameters(0, f_e);
		TF1* m_t2 = new TF1("", plot_transition, 0., 5000., 2);
		m_t2->SetParameters(0, f_t);
		m_t2->SetMaximum(.05);
		m_t2->SetMinimum(0.0);
		m_t2->SetNpx(10000);
		m_e2->SetNpx(10000);
		m_t2->SetLineColor(4);
		m_t2->Draw();
		m_e2->Draw("same");
	}
	else if(in_f == 2) {
		TF1* t_e = new TF1("t_e", plotP, 0., 35000., 2);
		t_e->SetParameters(f_t, f_e);
		TF1* t_m = new TF1("t_m", plotP, 0., 35000., 2);
		t_m->SetParameters(f_t, f_m);
		TF1* t_t = new TF1("m_t", plotP, 0., 35000, 2);
		t_t->SetParameters(f_t, f_t);
		t_e->SetNpx(10000);
		t_e->SetMinimum(0.0);
		t_e->SetMaximum(1.0);
		t_t->SetNpx(10000);
		t_m->SetNpx(10000);

		t_t->SetLineColor(8);
		t_m->SetLineColor(4);

		t_e->SetTitle("tau neutrino transition");
		t_e->Draw();
		t_m->Draw("same");
		t_t->Draw("same");
	}
}



