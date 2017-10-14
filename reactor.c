#include "osc_defs_NH.h"
using namespace TMath;

double P(double x) {
	double p;
	p = 1.0 - Power(cossq(t13), 2) * sinsq(2. * t12) * sinsq(1.2668 * dm2_21 * x)
			- cossq(t12) * sinsq(2. * t13) * sinsq(1.2668 * dm2_31 * x)
			- sinsq(t12) * sinsq(2. * t13) * sinsq(1.2668 * dm2_32 * x);
	return p;
}

double plot_P(double* x, double* par) {
	return P(x[0]);
}

void reactor() {
	TCanvas* c1 = new TCanvas();

	TF1* survival = new TF1("", plot_P, 0., 35000.);
	survival->SetTitle("#bar{#nu}_{e} #rightarrow #bar{#nu}_{e} (reactor)");
	survival->SetNpx(10000);
	survival->SetMaximum(1.0);
	survival->SetMinimum(0.0);

	survival->GetXaxis()->SetTitle("L/E (km/GeV)");
	survival->GetYaxis()->SetTitle("");

	survival->Draw();
}


