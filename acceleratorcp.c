#include "osc_defs_NH.h"
using namespace TMath;

//	cf pdg neutrino mixing intro for formulas
double acc_cp_P(double x) {
	double J = c13 * Sin(2 * t12) * Sin(2 * t13) * Sin(2 * t23)
				   * Sin(1.2668 * dm2_32 * x) * Sin(1.2668 * dm2_21 * x);
	
	double P1 = sinsq(t23) * sinsq(2 * t13) * sinsq(1.2668 * dm2_32 * x);
	double P2 = cossq(t23) * sinsq(2 * t13) * sinsq(1.2668 * dm2_21 * x);
	double P3 = J * Sin(d_cp) * Sin(1.2668 * dm2_32 * x);
	double P4 = J * Cos(d_cp) * Cos(1.2668 * dm2_32 * x);

	return P1 + P2 - P3 + P4;
}

double acc_cp_plot_P(double* x, double* par) {
	return acc_cp_P(x[0]);
}

void acceleratorcp() {
	TCanvas* c1 = new TCanvas();
	TF1* m_e = new TF1("", acc_cp_plot_P, 0., 35000.);
	m_e->SetTitle("accelerator #nu_{#mu} #rightarrow #nu_{e} (with CP, low mass scale)");
	m_e->SetNpx(10000);
	m_e->SetMaximum(1.0);
	m_e->SetMinimum(0.0);
	m_e->GetXaxis()->SetTitle("L/E (km/GeV)");
	m_e->Draw();
}


