#include "spectrum.h"
#include "helper.h"

// plotting and helper functions
void plot_initial(float* mu, float* antimu, float* e, float* antie);
void plot_normalized(float* mu, float* antimu, float* e, float* antie);
void plot_particles(float* mu, float* e);
void plot_antiparticles(float* antimu, float* antie);
void read_spectrum(initial_spectrum* s);
void plot_dc2(float* d_cp, float* dc2_mh_n, float* dc2_mh_i, float* dc2_cp, float* dc2_cp_ih, int N);
void do_dc2(const initial_spectrum* s);


// Mean delta chi-squared between a test spectrum and a "true" spectrum,
// as given in 1210.3651 p.9
float mean_dc2(const float* test, const float* tru, int N = 37) {
	float r = 0.;
	for (int i=0; i<N; ++i) {
		r += pow(test[i] - tru[i], 2) / N;
	}
	return r;
}

void nuspectrum(int anti=0) {
	// initial spectrum in the LBNF's neutrino mode
	initial_spectrum nu;
	// and in the antineutrino mode
	initial_spectrum antinu;


	gStyle->SetOptStat(0);

	// best fit spectra, with best fit angles, delta, normal hierarchy
	//spectrum best_fit_nu;
	// Populate can now cause segfault if we forget to assign a spectrum's hierarchy to a
	// valid object. Careful!
	//best_fit_nu.h = new hierarchy;
	//populate(best_fit_nu.h, NH, 0.5 * pi);

	// read data from root file
	read_spectrum(&nu);
	
	// switch particle/antiparticle spectrum around for the antinu mode
	for (int i=0; i<50; ++i) {
		antinu.mu[i] = nu.antimu[i];
		antinu.antimu[i] = nu.mu[i];
		antinu.e[i] = nu.antie[i];
		antinu.antie[i] = nu.e[i];
	}

	
	// plot the initial spectrum (un-normalized)
	plot_initial(nu.mu, nu.antimu, nu.e, nu.antie);
	//plot_initial(antinu.mu, antinu.antimu, antinu.e, antinu.antie);
	
	// Do the best-fit oscillation
	//oscillate(&nu, &best_fit_nu);

	
	// Main functions
	if(anti)
		do_dc2(&antinu);
	else
		do_dc2(&nu);
}

void do_dc2(const initial_spectrum* s) {
	const int N = 61;
	float d_cp[N];

	// Delta chi squared for the mass hierarchy assuming true normal hierarchy
	float dc2_mh_n[N] = {0};
	float dc2_mh_i[N] = {0};
	// Delta chi squared for delta_CP
	float dc2_cp[N] = {0};
	float dc2_cp_ih[N] = {0};


	// Calculate mean delta chi squared at each d_cp
	hierarchy nh0, nhpi;
	populate(&nh0, NH, 0.);
	populate(&nhpi, NH, pi);
	spectrum nh0s, nhpis;
	nh0s.h = &nh0;
	nhpis.h = &nhpi;
	oscillate(s, &nh0s);
	oscillate(s, &nhpis);
	

	hierarchy ih0, ihpi;
	populate(&ih0, IH, 0.);
	populate(&ihpi, IH, pi);
	spectrum ih0s, ihpis;
	ih0s.h = &ih0;
	ihpis.h = &ihpi;
	oscillate(s, &ih0s);
	oscillate(s, &ihpis);

	hierarchy nh[N], ih[N];
	spectrum nhs[N], ihs[N];

	// We can calculate all the spectra beforehand to save cpu when taking chisquared
	for (int i=0; i<N; ++i) {
		d_cp[i] = ((float)i / (float)(N-1) * 2. - 1.) * pi;
		nhs[i].h = &nh[i];
		ihs[i].h = &ih[i];
		populate(&nh[i], NH, d_cp[i]);
		populate(&ih[i], IH, d_cp[i]);

		oscillate(s, &nhs[i]);
		oscillate(s, &ihs[i]);

		Printf("%d/%d", i+1, N);
	}

	TCanvas* ccc = new TCanvas();
	//ccc->SetLogy();
	TH1* he = new TH1F("he", "", Nbins, 0.6, 8.);
	TH1* hmu = new TH1F("hmu", "", Nbins, 0.6, 8.);
	for (int i=0; i<Nbins; ++i) {
		he->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, nhs[N/2].e[i]);
		hmu->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, ihs[N/2].e[i]);
	}
	hmu->SetMaximum(80);
	hmu->Draw("hist");
	he->Draw("hist same");
	he->SetLineWidth(2);
	hmu->SetLineColor(2);
	
	// Loop over delta_CP values
	for (int i=0; i<N; ++i) {
		// for each d_cp we want to calculate the minimum corresponding delta chi squared
		float min_dc2_n = 1e9;
		float min_dc2_i = 1e9;

		// So we perform another loop over delta_CP and take the minimum chi squared we find
		for (int j=0; j<N; ++j) {
			float delta = ((float)j / (float)(N-1) * 2. - 1.) * pi;

			// When we assume normal hierarchy, the NH spectrum is fixed in i and we
			// search the IH with j
			float dc2_n = mean_dc2(ihs[j].e, nhs[i].e);
			// Vice versa
			float dc2_i = mean_dc2(nhs[j].e, ihs[i].e);
			min_dc2_n = min(dc2_n, min_dc2_n);
			min_dc2_i = min(dc2_i, min_dc2_i);
		}

		//dc2_mh_n[i] = mean_dc2(&nhs[i].e[3], &ihs[i].e[3], 47);
		//dc2_mh_i[i] = mean_dc2(&ihs[i].e[3], &ihs[i].e[3], 47);

		dc2_mh_n[i] = min_dc2_n;
		dc2_mh_i[i] = min_dc2_i;

		dc2_cp[i] = min(mean_dc2(nhs[i].e, nh0s.e), mean_dc2(nhs[i].e, nhpis.e));
		dc2_cp_ih[i] = min(mean_dc2(ihs[i].e, ih0s.e), mean_dc2(ihs[i].e, ihpis.e));
	}

	
	plot_dc2(d_cp, dc2_mh_n, dc2_mh_i, dc2_cp, dc2_cp_ih, N);
}

void plot_dc2(float* d_cp, float* dc2_mh_n, float* dc2_mh_i, float* dc2_cp, float* dc2_cp_ih, int N) {
	float x[N];
	float ymnh[N];
	float ymih[N];
	float ydnh[N];
	float ydih[N];
	for (int i=0; i<N; ++i) {
		x[i] = d_cp[i] / pi;
		ymnh[i] = sqrt(dc2_mh_n[i]);
		ymih[i] = sqrt(dc2_mh_i[i]);
		ydnh[i] = sqrt(dc2_cp[i]);
		ydih[i] = sqrt(dc2_cp_ih[i]);
	}


	TCanvas* c4 = new TCanvas("c4", "", 1000, 400);
	c4->Divide(2, 1);
	TPad* p = (TPad*)c4->GetPad(1);
	p->cd();
	TGraph* gdc2_mnh = new TGraph(N, x, ymnh);
	TGraph* gdc2_mih = new TGraph(N, x, ymih);
	gdc2_mnh->Draw();
	gdc2_mnh->SetTitle("");
	c4->SetTitle("MH sensitivity");
	gdc2_mnh->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mnh->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_mnh->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mnh->SetMinimum(0);
	//gdc2_mnh->SetMaximum(25);
	//gdc2_mnh->SetMaximum(25);
 	p = (TPad*)c4->GetPad(2);	
	p->cd();
	gdc2_mih->SetTitle("");
	gdc2_mih->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mih->SetMinimum(0);
	//gdc2_mih->SetMaximum(25);
	gdc2_mih->Draw();
	

	// Draw delta_CP mean delta chi squared
	TCanvas* c5 = new TCanvas("c5", "", 1000, 400);
	c5->Divide(2, 1);
	c5->cd(1);
	TGraph* gdc2_cp = new TGraph(N, x, ydnh);
	//gdc2_cp->Draw();
	gdc2_cp->Draw("");
	//gdc2_cp->SetMaximum(10);
	//gdc2_cp->SetMinimum(0);
	gdc2_cp->GetXaxis()->SetLimits(-1, 1);
	
	c5->cd(2);
	TGraph* gdc2_cp_ih = new TGraph(N, x, ydih);
	//gdc2_cp->Draw();
	gdc2_cp_ih->Draw("");
	//gdc2_cp->SetMaximum(10);
	//gdc2_cp->SetMinimum(0);
	gdc2_cp_ih->GetXaxis()->SetLimits(-1, 1);
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

