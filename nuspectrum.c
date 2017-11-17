#include "spectrum.h"
#include "helper.h"

// plotting and helper functions
void plot_initial(float* mu, float* antimu, float* e, float* antie);
void plot_normalized(float* mu, float* antimu, float* e, float* antie);
void plot_particles(float* mu, float* e);
void plot_antiparticles(float* antimu, float* antie);
void read_spectrum(initial_spectrum* s);


// Calculate (binned) chi squared of data set y against fit l.
// Formula from Barlow p. 105.
float chisq(const float* y, const float* l, int N) {
	float c2 = 0.;
	for (int i=0; i<N; ++i) {
		c2 += pow(y[i] - l[i], 2) / N;
	}
	return c2;
}

float mean_dc2(const float* nh, const float* ih, int N) {
	float r = 0.;
	for (int i=0; i<N; ++i) {
		r += pow(nh[i] - ih[i], 2) / ih[i];
	}
	return r;
}

void nuspectrum() {
	// initial spectrum in the LBNF's neutrino mode
	initial_spectrum nu;
	// and in the antineutrino mode
	initial_spectrum antinu;


	gStyle->SetOptStat(0);

	// best fit spectra, with best fit angles, delta, normal hierarchy
	// this will be used as our 'experimental' data when fitting with different parameters
	spectrum best_fit_nu;
	// Populate can now cause segfault if we forget to assign a spectrum's hierarchy to a
	// valid object. Careful!
	best_fit_nu.h = new hierarchy;
	populate(best_fit_nu.h, NH, 0.5 * TMath::Pi());

	// read data from root file
	read_spectrum(&nu);
	
	// switch particle/antiparticle spectrum around for the antinu mode
	for (int i=0; i<50; ++i) {
		antinu.mu[i] = nu.antimu[i];
		antinu.antimu[i] = nu.mu[i];
		antinu.e[i] = nu.antie[i];
		antinu.antie[i] = nu.e[i];
	}

	// perform average and plot it to see if it converges to mid point in the high
	// frequency region especially
	// apparently sampling 20 points from 0 to 10 looks pretty good
	int Nvals = 20;
	float xx[Nvals];
	float yy[Nvals];
	for (int i=0; i<Nvals; ++i) {
		// for each point we sample 20 times to get average
		int samples = 50;
		float avg_P = 0.;
		for (int j=0; j<samples; ++j) {
			float E = 1e-6 + (float)i / (float) Nvals * 10. + (float)j / (float)samples * 0.2;
			avg_P += P(f_m, f_m, E, best_fit_nu.h, false) / samples;
		}
		xx[i] = 10. * (float)i / (float)Nvals + 1e-6;
		yy[i] = avg_P;
	}
	Printf("t12 %f t23 %f t13 %f", best_fit_nu.h->t12, best_fit_nu.h->t23, best_fit_nu.h->t13);
	Printf("dm21 %e dm31 %e", best_fit_nu.h->dm2_mat[1 * 3 + 0], best_fit_nu.h->dm2_mat[2 * 3 + 0]);
	//TCanvas* c6 = new TCanvas();
	//TGraph* gP = new TGraph(Nvals, xx, yy);
	//gP->Draw();
	
	// plot the initial spectrum (un-normalized)
	//plot_initial(nu.mu, nu.antimu, nu.e, nu.antie);
	//plot_initial(antinu.mu, antinu.antimu, antinu.e, antinu.antie);
	
	// Do the best-fit oscillation
	oscillate(&nu, &best_fit_nu);

	// Draw best-fit parameter oscillated and normalized spectra
	/*TCanvas* c2 = new TCanvas();
	//c2->SetLogy();
	c2->SetGrid();
	TH1* hmu = new TH1F("#nu_{#mu}", "", 50, 0., 10.);
	TH1* hamu = new TH1F("#bar{#nu}_{#mu}", "", 50, 0., 10.);
	TH1* he = new TH1F("#nu_{e}", "", 50, 0., 10.);
	TH1* hae = new TH1F("#bar{#nu}_{e}", "", 50, 0., 10.);
	float norms[4] = {0};
	for (int i=0; i<50; ++i) {
		hmu->Fill(0.2 * i + 1e-2, best_fit_nu.mu[i]);
		norms[0] += 0.2 * best_fit_nu.mu[i];
		norms[1] += 0.2 * best_fit_nu.antimu[i];
		norms[2] += 0.2 * best_fit_nu.e[i];
		norms[3] += 0.2 * best_fit_nu.antie[i];
		hamu->Fill(0.2 * i + 1e-2, best_fit_nu.antimu[i]);
		he->Fill(0.2 * i + 1e-2, best_fit_nu.e[i]);
		hae->Fill(0.2 * i + 1e-2, best_fit_nu.antie[i]);
	}
	//Printf("%f %f %f %f", norms[0], norms[1], norms[2], norms[3]);
	hmu->SetLineWidth(2);
	he->SetLineWidth(2);
	hae->SetLineWidth(2);
	hamu->SetLineWidth(2);
	hmu->SetLineColor(1);
	he->SetLineColor(2);
	hae->SetLineColor(6);
	hamu->SetLineColor(4);
	hmu->SetMinimum(1e-1);
	hmu->Draw("HIST");
	hamu->Draw("HIST SAME");
	c2->BuildLegend();
		
	//TCanvas* cccc = new TCanvas();
	he->Draw("HIST SAME");

	//TCanvas* ccccc = new TCanvas();
	hae->Draw("HIST SAME");
	*/
	
	// As discussed in the CDR, we are exploring delta_CP(true) space and calculating
	// a chi squared for each. 
	const int N = 200;
	float d_cp[N];

	float dc2_mean[N] = {0};
	float dc2_cp[N] = {0};


	// Calculate mean delta chi squared at each d_cp
	hierarchy nh0, nhpi;
	populate(&nh0, NH, 0.);
	populate(&nhpi, NH, TMath::Pi());
	spectrum nh0s, nhpis;
	nh0s.h = &nh0;
	nhpis.h = &nhpi;
	oscillate(&nu, &nh0s);
	oscillate(&nu, &nhpis);

	hierarchy nh, ih;
	spectrum nhs, ihs;
	populate(&ih, IH);
	ihs.h = &ih;
	oscillate(&nu, &ihs);

	hierarchy bestnh;
	spectrum bestnhs;
	populate(&bestnh, NH, 0.55 * TMath::Pi());
	bestnhs.h = &bestnh;
	oscillate(&nu, &bestnhs);

	for (int i=0; i<N; ++i) {
		d_cp[i] = ((float)i / (float)(N-1) * 2. - 1.) * TMath::Pi();
		populate(&nh, NH, d_cp[i]);
		populate(&ih, IH, d_cp[i]);

		nhs.h = &nh;
		ihs.h = &ih;
		oscillate(&nu, &nhs);
		oscillate(&nu, &ihs);

		dc2_mean[i] = mean_dc2(nhs.e, ihs.e, 50);
		//dc2_mean[i] = mean_dc2(nhs.e, bestnhs.e, 50);
		// TODO plot  same for ih
		dc2_cp[i] = min(mean_dc2(nhs.e, nh0s.e, 50), mean_dc2(nhs.e, nhpis.e, 50));
		//dc2_cp[i] = mean_dc2(nhs.e, nhpis.e, 50);
	}

	TCanvas* c3 = new TCanvas();
	
	const int p = 50;
	populate(&nh, NH, 1. * TMath::Pi());
	populate(&ih, IH, 1. * TMath::Pi());
	nhs.h = &nh;
	ihs.h = &ih;
	oscillate(&nu, &nhs);
	oscillate(&nu, &ihs);

	TH1* hnh = new TH1F("hnh", "", 50, 0., 10.);
	TH1* hih = new TH1F("hih", "", 50, 0., 10.);
	for (int i=0; i<50; ++i) {
		hnh->Fill(i*0.2 + 1e-3, nhs.e[i]);
		hih->Fill(i*0.2 + 1e-3, ihs.e[i]);
	}

	hnh->SetLineWidth(2);
	hnh->SetLineColor(2);
	hih->SetLineWidth(2);
	hnh->Draw("hist");
	hih->Draw("hist same");

	float x[N];
	float y[N];
	float y2[N];
	for (int i=0; i<N; ++i) {
		x[i] = d_cp[i] / TMath::Pi();
		y[i] = sqrt(dc2_mean[i]);
		y2[i] = sqrt(dc2_cp[i]);
	}

	TCanvas* c4 = new TCanvas();
	TGraph* gdc2_mean = new TGraph(N, x, y);
	gdc2_mean->Draw();
	gdc2_mean->SetTitle("MH sensitivity");
	gdc2_mean->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mean->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_mean->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mean->SetMinimum(0);
	gdc2_mean->SetMaximum(25);

	TCanvas* c5 = new TCanvas();
	TGraph* gdc2_cp = new TGraph(N, x, y2);
	gdc2_cp->Draw("");
	gdc2_cp->SetMaximum(10);
	gdc2_cp->SetMinimum(0);
	gdc2_cp->GetXaxis()->SetLimits(-1, 1);


	/*memset(&nh, 0, sizeof(hierarchy));
	int dcpidx = 25;
	populate(&nh, IH, d_cp[dcpidx]); 
	Printf("dcp / pi = %f", d_cp[dcpidx] / TMath::Pi());
	spectrum tests;
	tests.h = &nh;
	oscillate(&nu, &tests);

	TCanvas* cc = new TCanvas();
	TH1* test0h = new TH1F("test0h", "", 50, 0., 10.);
	TH1* test5h = new TH1F("test5h", "", 50, 0., 10.);
	for(int i=0; i<50; ++i) {
		test0h->Fill(0.2 * i + 1e-3, best_fit_nu.e[i]);
		//Printf("test %f", tests.e[i]);
		test5h->Fill(0.2 * i + 1e-3, tests.e[i]);
		//Printf("dif %f", pow(best_fit_nu.e[i] - tests.e[i], 1));
	}
	//cc->SetLogy();
	test0h->SetLineWidth(2);
	test0h->SetLineColor(2);
	test5h->SetLineWidth(2);
	test0h->Draw("HIST");
	test5h->Draw("SAME HIST");
	*/
}

void read_spectrum(initial_spectrum* s) {
	// read from root file
	TFile* f = TFile::Open("data/spectrum.root");

	// get the tree
	TTree* t = (TTree*)f->Get("spectrum");

	// use tree readers for each branch
	TTreeReader r("spectrum", f);
	TTreeReaderValue<float> mu(r, "mu");
	TTreeReaderValue<float> antimu(r, "antimu");
	TTreeReaderValue<float> e(r, "e");
	TTreeReaderValue<float> antie(r, "antie");

	// fill up the arrays
	int i = 0;
	while(r.Next()) {
		s->mu[i] = *mu;
		s->antimu[i] = *antimu;
		s->e[i] = *e;
		s->antie[i] = *antie;	

		++i;
	}
}

void plot_initial(float* mu, float* antimu, float* e, float* antie) {
	TCanvas* c1 = new TCanvas("c1", "", 600, 600);
	// set log y axis and no statistics legend
	c1->SetLogy();	
	c1->SetTicks();
	gStyle->SetOptStat(0);

	TH1* h1 = new TH1F("#nu_{#mu}", "", 50, 0., 10.);
	TH1* h2 = new TH1F("#bar{#nu}_{#mu}", "", 50, 0., 10.);
	TH1* h3 = new TH1F("#nu_{e}", "", 50, 0., 10.);
	TH1* h4 = new TH1F("#bar{#nu}_{e}", "", 50, 0., 10.);
	
	// fill histograms manually with values from the spectrum
	for(int i=0; i<50; ++i) {
		h1->Fill(i * 0.2, mu[i]);
		h2->Fill(i * 0.2, antimu[i]);
		h3->Fill(i * 0.2, e[i]);
		h4->Fill(i * 0.2, antie[i]);
	}

	
	// plot the initial spectrum plots
	h1->SetLineWidth(2);
	h2->SetLineWidth(2);
	h3->SetLineWidth(2);
	h4->SetLineWidth(2);

	h1->GetXaxis()->SetTitle("Energy / GeV");
	h1->GetYaxis()->SetTitle("Flux");

	h1->SetLineColor(1);
	h1->SetMinimum(1e6);
	h1->SetMaximum(1e10);
	h1->Draw("HIST");
	h2->SetLineColor(4);	
	h2->Draw("SAME HIST");
	h3->SetLineColor(2);
	h3->Draw("SAME HIST");
	h4->SetLineColor(6);
	h4->Draw("SAME HIST");

	gPad->SetGrid();
	gPad->BuildLegend();

}

void plot_normalized(float* mu, float* antimu, float* e, float* antie) {
	
	TH1* h1 = new TH1F("h5", "", 50, 0., 10.);
	TH1* h2 = new TH1F("h6", "", 50, 0., 10.);
	TH1* h3 = new TH1F("h7", "", 50, 0., 10.);
	TH1* h4 = new TH1F("h8", "", 50, 0., 10.);
	
	// fill histograms manually with values from the spectrum
	for(int i=0; i<50; ++i) {
		h1->Fill(i * 0.2, mu[i]);
		h2->Fill(i * 0.2, antimu[i]);
		h3->Fill(i * 0.2, e[i]);
		h4->Fill(i * 0.2, antie[i]);
	}

	
	// plot the initial spectrum plots
	h1->SetLineColor(1);
	h1->SetMaximum(1.0);
	h1->SetMinimum(3e-3);
	h1->Draw("HIST");
	h2->SetLineColor(4);	
	h2->Draw("SAME HIST");
	h3->SetLineColor(2);
	h3->Draw("SAME HIST");
	h4->SetLineColor(6);
	h4->Draw("SAME HIST");
}

void plot_particles(float* mu, float* e) {
	TH1* h1 = new TH1F("Initial #nu_{#mu}", "", 50, 0., 10.);
	TH1* h2 = new TH1F("Initial #nu_{e}", "", 50, 0., 10.);
	
	// fill histograms manually with values from the spectrum
	for(int i=0; i<50; ++i) {
		h1->Fill(i * 0.2, mu[i]);
		h2->Fill(i * 0.2, e[i]);
	}

	
	// plot the initial spectrum plots
	h1->SetLineWidth(2);
	h2->SetLineWidth(2);

	h1->SetLineStyle(2);
	h2->SetLineStyle(2);

	h1->SetLineColor(1);
	h1->Draw("HIST SAME");
	h2->SetLineColor(2);	
	h2->Draw("SAME HIST");
}

void plot_antiparticles(float* antimu, float* antie) {
	TH1* h1 = new TH1F("Initial #bar{#nu}_{#mu}", "", 50, 0., 10.);
	TH1* h2 = new TH1F("Initial #bar{#nu}_{e}", "", 50, 0., 10.);
	
	// fill histograms manually with valuanties from the spectrum
	for(int i=0; i<50; ++i) {
		h1->Fill(i * 0.2, antimu[i]);
		h2->Fill(i * 0.2, antie[i]);
	}

	
	// plot thantie initial spectrum plots
	h1->SetLineWidth(2);
	h2->SetLineWidth(2);

	h1->SetLineStyle(2);
	h2->SetLineStyle(2);

	h1->SetLineColor(4);
	h1->Draw("SAME HIST");
	h2->SetLineColor(6);
	h2->Draw("SAME HIST");
}

