// This file enables us to compare transition probability in the presence
// or absence of matter effects on the oscillated neutrinos.
// We will plot on the same canvas the transition probability from
// muon neutrino to electron neutrino with and without matter effects,
// as a function of energy. Hence we need the functions defined in spectrum.c
// i.e. P and P_me
#include "spectrum.h"

void matter() {
	// We have to pick a hierarchy that effectively shows the difference.
	hierarchy h;
	//double d_cp = -pi/2.;
	populate(&h, NH);

	// number of samples
	const int N = 100000;
	double x[N];
	double y[N];
	double y_me[N];
	
	for (int i=0; i<N; ++i) {
		// Range x (E) from 0.1 to 10 GeV
		x[i] = 1e-1 + (double)i / (double)N * (10. - 1e-1);

		y[i] = P(f_m, f_e, x[i], &h, false);
		y_me[i] = P_me(f_m, f_e, x[i], &h, false);
	}

	TCanvas* c = new TCanvas("c1", "c1", 700, 500);
	c->SetFillColor(ci[CI_BACKGROUND]);
	c->SetLogx();
	TGraph* gp = new TGraph(N, x, y);
	TGraph* gpme = new TGraph(N, x, y_me);

	gp->SetTitle("Vacuum");
	gpme->SetTitle("Earth mantle");
	gp->GetXaxis()->SetLimits(1e-1, 10);
	gp->SetLineWidth(3);
	gpme->SetLineWidth(3);
	gp->SetLineColor(ci[CI_4]);
	gpme->SetLineColor(ci[CI_3]);

	gp->GetXaxis()->SetTitle("E (GeV)");
	gp->GetYaxis()->SetTitle("P(#nu_{#mu} #rightarrow #nu_{e})");

	gp->SetMarkerColor(ci[CI_4]);
	gpme->SetMarkerColor(ci[CI_3]);
	gp->Draw("");
	gpme->Draw("same");

	c->BuildLegend();

}
