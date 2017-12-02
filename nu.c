#include "osc_defs.h"
#include "helper.h"


//	hierarchies to be used: combinations of normal and inverted 
//	those first two are the best fit delta hierarchies
hierarchy nhd, ihd;
//	four hierarchies, one for each delta_cp = 0, 90, 180, 270 degrees
hierarchy nh[3];
hierarchy ih[3];


double plot_P_nhd(double* x, double* par) {
	return P(f_m, f_e, x[0], &nhd, (bool)par[1]);
}

double plot_P_ihd(double* x, double* par) {
	return P(f_m, f_e, x[0], &ihd, (bool)par[1]);
}

double plot_P_nh(double* x, double* par) {
	return P_me(f_m, f_e, x[0], &nh[(int)par[0]], (bool)par[1]);
}

double plot_P_ih(double* x, double* par) {
	return P(f_m, f_e, x[0], &ih[(int)par[0]], (bool)par[1]);
}

void plot_P(bool anti);

void nu() {
	for(int i=0; i<3; ++i) {
		populate(&nh[i], NH);
		//	set delta to a multiple of 90 degrees and recalculate MNS
		nh[i].d_cp = TMath::Pi() / 2. * (i - 1);
		//nh[0].d_cp = -0.8 * pi;
		populate_common(&nh[i]);
	}
	
	//	styling of plot titles
	gStyle->SetTitleAlign(33);
	gStyle->SetTitleX(.999);
	gStyle->SetTitleY(.5);
	
	//	UNCOMMENT FOR BEST FIT DELTA
	/*populate(&nhd, NH);
	populate(&ihd, IH);
	TCanvas* c2 = new TCanvas("c2", "", 800, 400);
	TF1* realdnh = new TF1("", plot_P_nhd, 0.5, 4.5);
	TF1* realdih = new TF1("", plot_P_ihd, 0.5, 4.5);
	realdnh->SetMaximum(0.1);
	realdnh->SetNpx(10000);
	realdih->SetNpx(10000);
	realdnh->SetLineColor(4);
	realdih->SetLineColor(8);
	realdnh->SetTitle("best fit #delta_{CP}");
	realdnh->Draw();
	realdih->Draw("same");
	*/
	
	//	delta != 0 plots
	// Generate random name for Canvas to avoid overwriting when re-running
	char name[64];
	srand(time(NULL));
	int r = rand();
	snprintf(name, 64, "c%d", r);
	TCanvas* c1 = new TCanvas(name, "", 1200, 500);

	c1->Divide(2, 1);
	c1->GetPad(1)->SetLogx();
	c1->GetPad(2)->SetLogx();
	//c1->GetPad(1)->SetTickx(2);
	


	c1->cd(1);
	plot_P(false);


	c1->cd(2);
	plot_P(true);


	c1->GetPad(1)->BuildLegend();
	c1->GetPad(2)->BuildLegend();
	c1->GetPad(1)->SetTicks();
	c1->GetPad(2)->SetTicks();

}

void plot_P(bool anti) {
	TF1* f[3];
	Color_t cols[3] = { 4, 2, 3 };

	for(int i=0; i<3; ++i) {
		char title[64];
		snprintf(title, 64, "#delta_{CP} = %d deg", (i-1) * 90);
		f[i] = new TF1(title, plot_P_nh, 0.1, 10., 2);
		f[i]->SetParameter(0, (double)i);
		f[i]->SetParameter(1, (double)anti);

		f[i]->SetTitle("");
		f[i]->SetMaximum(.2);
		f[i]->SetMinimum(0.0);
		f[i]->SetNpx(10000);
		f[i]->GetXaxis()->SetTitle("E (GeV)");

		f[i]->SetLineWidth(0);

		f[i]->SetFillColor(cols[i]);
		f[i]->SetFillStyle(1001);
	}
	// need to change order in drawing depending on anti or not
	if (anti) {
		f[2]->Draw();
		f[1]->Draw("same");
		f[0]->Draw("same");
	}
	else {
		f[0]->Draw();
		f[1]->Draw("same");
		f[2]->Draw("same");
	}
}
