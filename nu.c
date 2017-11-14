#include "osc_defs.h"
#include "helper.h"

//	baseline for experiment: DUNE has L = 1300km
//const double L = 1300.0;

//	hierarchies to be used: combinations of normal and inverted 
//	those first two are the best fit delta hierarchies
hierarchy nhd, ihd;
//	four hierarchies, one for each delta_cp = 0, 90, 180, 270 degrees
hierarchy nh[4];
hierarchy ih[4];

double P(double E, hierarchy* h) {
	double p = 0.;

	//	initial and final states
	const int a = f_m;
	const int b = f_e;

	for(int i=0; i<3; ++i) {
		p += pow(abs(h->MNS[a*3 + i] * std::conj(h->MNS[b*3 + i])), 2);
	}

	//Printf("%f", E);
	for(int i=0; i<3; ++i) {
		for(int j=0; j<3; ++j) {
			if(j > i) {
				p += 2. * std::real(h->MNS[a*3 + i] * std::conj(h->MNS[a*3 + j]) 
						* std::conj(h->MNS[b*3 + i]) * h->MNS[b*3 + j]
						* std::exp(complex<double>(-2.i * 1.2668 *  h->dm2_mat[i*3 + j] * L/E)));
			}
		}
	}
	return p;
}

double plot_P_nhd(double* x, double* par) {
	return P(x[0], &nhd);
}

double plot_P_ihd(double* x, double* par) {
	return P(x[0], &ihd);
}

double plot_P_nh(double* x, double* par) {
	return P(x[0], &nh[(int)par[0]]);
}

double plot_P_ih(double* x, double* par) {
	return P(x[0], &ih[(int)par[0]]);
}

void nu() {
	for(int i=0; i<4; ++i) {
		populate(&nh[i], NH);
		populate(&ih[i], IH);
		//	set delta to a multiple of 90 degrees and recalculate MNS
		nh[i].d_cp = TMath::Pi() / 2. * i;
		ih[i].d_cp = TMath::Pi() / 2. * i;
		populate_common(&nh[i]);
		populate_common(&ih[i]);
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
	TCanvas* c1 = new TCanvas("c1", "", 800, 900);
	c1->Divide(1, 4);
	c1->GetPad(1)->SetTickx(2);
	TF1* me_nh[4];
	TF1* me_ih[4];

	for(int i=0; i<4; ++i) {
		c1->cd(i+1);
		if(i > 0)
			gPad->SetTopMargin(0.);
		if(i < 3)
			gPad->SetBottomMargin(0.);
		gPad->SetLeftMargin(.05);
		gPad->SetPad(0., (3.-i) / 4., 1., (4.-i) / 4.);

		me_nh[i] = new TF1("", plot_P_nh, 0.5, 4.5, 1);
		me_ih[i] = new TF1("", plot_P_ih, 0.5, 10.0, 1);
		me_nh[i]->SetParameter(0, (double)i);
		me_ih[i]->SetParameter(0, (double)i);

		me_nh[i]->SetMaximum(.1);
		me_nh[i]->SetMinimum(0.0);
		me_nh[i]->SetNpx(10000);
		me_nh[i]->GetXaxis()->SetTitle("E (GeV)");

		char title[64];
		snprintf(title, 64, "#delta_{CP} = %d deg", 90 * i);
		me_nh[i]->SetTitle(title);
		me_nh[i]->SetLineColor(4);
		me_nh[i]->Draw();
		me_ih[i]->SetNpx(10000);
		me_ih[i]->SetLineColor(8);
		
		me_ih[i]->Draw("same");
	}

}
