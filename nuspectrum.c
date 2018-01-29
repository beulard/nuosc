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
float mean_dc2(const double* test, const double* tru, int N = 37) {
	float r = 0.;
	for (int i=0; i<N; ++i) {
		r += pow(test[i] - tru[i], 2) / tru[i];
	}
	return r;
}

void nuspectrum(int anti=0) {
	gStyle->SetOptStat(0);

	// initial spectrum in the LBNF's neutrino mode
	initial_spectrum nu;
	// and in the antineutrino mode
	initial_spectrum antinu;

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
	//plot_initial(nu.mu, nu.antimu, nu.e, nu.antie);
	//plot_initial(antinu.mu, antinu.antimu, antinu.e, antinu.antie);
	
	// Do the best-fit oscillation
	//oscillate(&nu, &best_fit_nu);
	

	// TODO
	// Gaussian smearing on the reconstructed energy
	// and try the CDR transition probability

	
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

	 // To plot a test spectrum
	TCanvas* ccc = new TCanvas();
	//ccc->SetLogy();
	TH1* hnh = new TH1F("hnh", "", Nbins, 0.6, 8.);
	TH1* hih = new TH1F("hih", "", Nbins, 0.6, 8.);
	int index = (int)((.0 / 2. + .5) * N);

	Printf("%f", d_cp[index] / pi);
	for (int i=0; i<Nbins; ++i) {
		hnh->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, nhs[index].events[E_SIGNAL][i]);
		hih->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, ihs[index].events[E_SIGNAL][i]);
	}
	hnh->Draw("hist");
	hnh->SetMaximum(110);
	hih->Draw("hist same");
	hnh->SetLineWidth(2);
	hih->SetLineColor(2);
	
	// Loop over delta_CP values
	for (int i=0; i<N; ++i) {
		// for each d_cp we want to calculate the minimum corresponding delta chi squared
		float min_dc2_n = 1e28;
		float min_dc2_i = 1e28;

		// So we perform another loop over delta_CP and take the minimum chi squared we find
		for (int j=0; j<N; ++j) {
			float delta = ((float)j / (float)(N-1) * 2. - 1.) * pi;

			// When we assume normal hierarchy, the NH spectrum is fixed in i and we
			// search the IH with j
			float dc2_n = mean_dc2(ihs[j].events[E_SIGNAL], nhs[i].events[E_SIGNAL]);
			// Vice versa
			float dc2_i = mean_dc2(nhs[j].events[E_SIGNAL], ihs[i].events[E_SIGNAL]);

			min_dc2_n = min(dc2_n, min_dc2_n);
			min_dc2_i = min(dc2_i, min_dc2_i);

		}

		dc2_mh_n[i] = min_dc2_n;
		dc2_mh_i[i] = min_dc2_i;

		dc2_cp[i] = min(mean_dc2(nhs[i].events[E_SIGNAL], nh0s.events[E_SIGNAL]),
			   			mean_dc2(nhs[i].events[E_SIGNAL], nhpis.events[E_SIGNAL]));
		dc2_cp_ih[i] = min(mean_dc2(ihs[i].events[E_SIGNAL], ih0s.events[E_SIGNAL]), 
						   mean_dc2(ihs[i].events[E_SIGNAL], ihpis.events[E_SIGNAL]));
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
	c4->SetFillColor(ci[CI_BACKGROUND]);
	c4->Divide(2, 1);
	
	TPad* p = (TPad*)c4->GetPad(1);
	p->cd();
	TGraph* gdc2_mnh = new TGraph(N, x, ymnh);
	TGraph* gdc2_mih = new TGraph(N, x, ymih);
	gdc2_mnh->Draw();
	gdc2_mnh->SetTitle("True normal hierarchy");
	gdc2_mih->SetTitle("True inverted hierarchy");
	c4->SetTitle("MH sensitivity");
	gdc2_mnh->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mnh->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_mih->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mih->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_mnh->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mnh->SetLineColor(ci[CI_NH]);
	gdc2_mih->SetLineColor(ci[CI_IH]);
	gdc2_mnh->SetMarkerColor(ci[CI_NH]);
	gdc2_mih->SetMarkerColor(ci[CI_IH]);
	gdc2_mnh->SetMinimum(0);
	gdc2_mnh->SetMaximum(25);
	gdc2_mih->SetMaximum(25);
	gdc2_mnh->SetLineWidth(3);
	gdc2_mih->SetLineWidth(3);
 	p = (TPad*)c4->GetPad(2);	
	p->cd();
	gdc2_mih->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mih->SetMinimum(0);
	//gdc2_mih->SetMaximum(25);
	gdc2_mih->Draw();
	

	// Draw delta_CP mean delta chi squared
	TCanvas* c5 = new TCanvas("c5", "", 1000, 400);
	c5->SetFillColor(ci[CI_BACKGROUND]);
	c5->Divide(2, 1);
	c5->cd(1);
	TGraph* gdc2_cp = new TGraph(N, x, ydnh);
	//gdc2_cp->Draw();
	gdc2_cp->Draw();
	gdc2_cp->SetTitle("True normal hierarchy");
	//gdc2_cp->SetMaximum(10);
	//gdc2_cp->SetMinimum(0);
	gdc2_cp->GetXaxis()->SetLimits(-1, 1);
	gdc2_cp->SetLineColor(ci[CI_NH]);
	gdc2_cp->SetLineWidth(3);
	gdc2_cp->SetMarkerColor(ci[CI_NH]);
	gdc2_cp->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_cp->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	
	c5->cd(2);
	TGraph* gdc2_cp_ih = new TGraph(N, x, ydih);
	//gdc2_cp->Draw();
	gdc2_cp_ih->Draw("");
	gdc2_cp_ih->SetTitle("True inverted hierarchy");
	gdc2_cp_ih->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_cp_ih->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	//gdc2_cp->SetMaximum(10);
	//gdc2_cp->SetMaximum(10);
	//gdc2_cp->SetMinimum(0);
	gdc2_cp_ih->SetLineColor(ci[CI_IH]);
	gdc2_cp_ih->SetLineWidth(3);
	gdc2_cp_ih->SetMarkerColor(ci[CI_IH]);
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

