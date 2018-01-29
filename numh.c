#include "osc_defs.h"
#include "helper.h"


hierarchy nh;
hierarchy ih;


double plot_P(double* x, double* par) {
	hierarchy* h;
	if((int)par[0] == 0)
		h = &nh;
	else
		h = &ih;
	double p = P_me(f_m, f_e, x[0], h, false);
	// Fix plots overflowing out of their frame/axes
	//if (p > 0.2)
	//	p = 0.2;
	return p;
}


void plot_P();

void numh() {
	populate(&nh, NH);
	populate(&ih, IH);
	
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
	TCanvas* c1 = new TCanvas(name, "", 800, 500);
	c1->SetFillColor(ci[CI_BACKGROUND]);

	//c1->Divide(2, 1);
	//c1->GetPad(1)->SetLogx();
	//c1->GetPad(2)->SetLogx();
	c1->SetLogx();
	


	//c1->cd(1);
	plot_P();


	//c1->cd(2);
	//plot_P(true);


	//c1->GetPad(1)->BuildLegend();
	//c1->GetPad(2)->BuildLegend();
	//c1->GetPad(1)->SetTicks();
	//c1->GetPad(2)->SetTicks();
	
	c1->BuildLegend();
	c1->SetTicks();

}

void plot_P() {
	TF1* f[2];
	//Color_t cols[3] = { 4, 2, 3 };
	int cols[2] = { ci[CI_1], ci[CI_2] };

	for(int i=0; i<2; ++i) {
		char title[64];
		if (i==0)
			strcpy(title, "Normal hierarchy");
		else
			strcpy(title, "Inverted hierarchy");
		f[i] = new TF1(title, plot_P, 0.1, 10., 1);
		f[i]->SetParameter(0, (double)i);

		f[i]->SetTitle("");
		f[i]->SetMaximum(.2);
		f[i]->SetMinimum(0.0);
		f[i]->SetNpx(20000);
		f[i]->GetXaxis()->SetTitle("E (GeV)");
		f[i]->GetYaxis()->SetTitle("P(#nu_{#mu} #rightarrow #nu_{e})");

		f[i]->SetLineColor(cols[i]);
		f[i]->SetLineWidth(3);

		//f[i]->SetFillColor(cols[i]);
		//f[i]->SetFillStyle(1001);
	}
	f[0]->Draw("");
	f[1]->Draw("same");
}
