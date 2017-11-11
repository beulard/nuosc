#include "spectrum.h"

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
		c2 += pow(y[i] - l[i], 2) / l[i];
	}
	return c2;
}

// calculates the chi squared given a hierarchy and a measured spectrum 
float spectrum_chisq(hierarchy* h, const initial_spectrum* is, const spectrum* data) {
	spectrum s;
	// assign hierarchy
	s.h = *h;
	
	// propagate spectrum 
	propagate(is, &s);
	
	// now calculate the chi-squared
	return chisq(data->e, s.e, 50);
}

void nuspectrum() {
	// initial spectrum in the LBNF's neutrino mode
	initial_spectrum nu;
	// and in the antineutrino mode
	initial_spectrum antinu;
	// TODO Cleanup antinu things,
	// let's start with simple
	// we'll do only the neutrino mode, but we still need the hierarchy stuff
	// and we can probably salvage some of the deltaCP exploration stuff as well


	gStyle->SetOptStat(0);

	// best fit spectra, with best fit angles, delta, normal hierarchy
	// this will be used as our 'experimental' data when fitting with different parameters
	spectrum best_fit_nu;
	populate(&best_fit_nu.h, NH);

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
	
	//normalize(&nu);

	propagate(&nu, &best_fit_nu);

	// Draw propagated and normalized spectrum
	TCanvas* c2 = new TCanvas();
	c2->SetLogy();
	TH1* hmu = new TH1F("hmup", "", 50, 0., 10.);
	TH1* hamu = new TH1F("hamup", "", 50, 0., 10.);
	TH1* he = new TH1F("hep", "", 50, 0., 10.);
	TH1* hae = new TH1F("haep", "", 50, 0., 10.);
	for (int i=0; i<50; ++i) {
		hmu->Fill(0.2 * i + 1e-2, best_fit_nu.mu[i]);
		hamu->Fill(0.2 * i + 1e-2, best_fit_nu.antimu[i]);
		he->Fill(0.2 * i + 1e-2, best_fit_nu.e[i]);
		hae->Fill(0.2 * i + 1e-2, best_fit_nu.antie[i]);
	}
	hmu->Draw("HIST");
	he->Draw("HIST SAME");
	hae->Draw("HIST SAME");
	hamu->Draw("HIST SAME");
	// For now we need to take one of our propagations as the data set, so we'll pick the one
	// with best fit parameters.
	// Then we'll fit that to a bunch of predictions, changing both the hierarchy and 
	// delta.
	// Then we'll plot the chi-squared difference against delta_cp on two plots.
	// one for each hierarchy.
	//
	// Half of the exposure is done in the neutrino mode, and the other half in
	// antineutrino mode.
	// We want to simultaneously fit the oscillated spectra to the experimental data.
	
	// As discussed in the CDR, we are exploring delta_CP(true) space and calculating
	// a chi squared for each. 
	int N = 100;
	float d_cp[N];

	float c2_nu_nh[N];
	float c2_nu_ih[N];


	// Here we calculate a chi squared between best fit spectrum and 
	// explored spectrum, for a normal and inverted hierarchy propagation
	// but wait we have different normalizations for nh and ih
	for (int i=0; i<N; ++i) {
		d_cp[i] = (float)i / (float)N * 2. * TMath::Pi();
		hierarchy nh;
		populate(&nh, NH, d_cp[i]);
		c2_nu_nh[i] = spectrum_chisq(&nh, &nu, &best_fit_nu);

		hierarchy ih;
		populate(&ih, IH, d_cp[i]);
		c2_nu_ih[i] = spectrum_chisq(&ih, &nu, &best_fit_nu);
	}

	// We also need the delta_cp test values, 0 and pi
	hierarchy nh_0;
	hierarchy nh_pi;
	hierarchy ih_0;
	hierarchy ih_pi;
	populate(&nh_0, NH, 0.);
	populate(&nh_pi, NH, TMath::Pi());
	populate(&ih_0, IH, 0.);
	populate(&ih_pi, IH, TMath::Pi());
	
	float c2_nu_nh_0 = spectrum_chisq(&nh_0, &nu, &best_fit_nu);
	float c2_nu_nh_pi = spectrum_chisq(&nh_pi, &nu, &best_fit_nu);
	float c2_nu_ih_0 = spectrum_chisq(&ih_0, &nu, &best_fit_nu);
	float c2_nu_ih_pi = spectrum_chisq(&ih_pi, &nu, &best_fit_nu);


	float delta_c2_cpv_nh_nu[N];
	float delta_c2_cpv_ih_nu[N];

	float delta_c2_mh_nu[N];

	for (int i=0; i<N; ++i) {
		//Printf("%f %f", c2_nu_test_0 - c2_nu_nh[i], c2_nu_test_pi - c2_nu_nh[i]);
		delta_c2_cpv_nh_nu[i] = (min(c2_nu_nh_0 - c2_nu_nh[i],
				    		    	c2_nu_nh_pi - c2_nu_nh[i]));
		
		delta_c2_cpv_ih_nu[i] = (min(c2_nu_ih_0 - c2_nu_ih[i],
									c2_nu_ih_pi - c2_nu_ih[i]));

		delta_c2_mh_nu[i] = c2_nu_ih[i] - c2_nu_nh[i];
		//Printf("%f", delta_c2_mh_nu[i]);
	}

	// This canvas contains the plots for the neutrino mode: delta_CPV and delta_MH
	TCanvas* c4 = new TCanvas("Neutrino mode", "", 900, 1000);
	c4->Divide(1, 2);
	TPad* pad = (TPad*)c4->GetPad(2);
	pad->cd();
	TH1* hdc2_nh = new TH1F("normal hierarchy", "", N, 0., 2.);
	TH1* hdc2_ih = new TH1F("inverted hierarchy", "", N, 0., 2.);
	for (int i=0; i<N; ++i) {
		hdc2_nh->Fill((float)i / (float)N * 2. + 1e-3, delta_c2_cpv_nh_nu[i]);
		hdc2_ih->Fill((float)i / (float)N * 2. + 1e-3, delta_c2_cpv_ih_nu[i]);
	}
	hdc2_nh->Draw("HIST");
	hdc2_ih->SetLineColor(2);
	hdc2_ih->Draw("HIST SAME");
	// formatting
	pad->BuildLegend();
	pad->SetTicks();
	pad->SetGrid();
	hdc2_nh->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	hdc2_nh->GetYaxis()->SetTitle("| #Delta#chi^{2}_{CPV} |");
	hdc2_nh->SetTitle("#nu mode");
	///

	pad = (TPad*)c4->GetPad(1);
	pad->cd();
	TH1* hdc2_mh = new TH1F("hdc2_mh", "", N, 0., 2.);
	TH1* hdc2_cnh = new TH1F("hdc2_cnh", "", N, 0., 2.);
	TH1* hdc2_cih = new TH1F("hdc2_cih", "", N, 0., 2.);
	for (int i=0; i<N; ++i) {
		hdc2_mh->Fill((float)i / (float)N * 2. + 1e-3, delta_c2_mh_nu[i]);
		hdc2_cnh->Fill((float)i / (float)N * 2. + 1e-3, c2_nu_nh[i]);
		hdc2_cih->Fill((float)i / (float)N * 2. + 1e-3, c2_nu_ih[i]);
	}
	hdc2_mh->Draw("HIST");
	hdc2_cnh->Draw("HIST SAME");
	hdc2_cih->Draw("HIST SAME");
	hdc2_mh->SetMinimum(0);
	pad->BuildLegend();
	hdc2_cnh->SetLineColor(2);
	// formatting
	pad->SetTicks();
	pad->SetGrid();
	hdc2_mh->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	hdc2_mh->GetYaxis()->SetTitle("#Delta#chi^{2}_{MH}");
	///



	// now plot the final flux
	// this plot is now obsolete: just for testing purposes
	/*
	TCanvas* c3 = new TCanvas("c3", "", 800, 800);
	//c3->Divide(1, 2);
	c3->SetLogy();
	c3->SetTicks();


	
	TH1* hmu = new TH1F("Predicted #nu_{#mu} + #bar{#nu}_{#mu}", "", 50, 0., 10.);
	TH1* he = new TH1F("Predicted #nu_{e} + #bar{#nu}_{e}", "", 50, 0., 10.);
	TH1* htau = new TH1F("Predicted #nu_{#tau} + #bar{#nu}_{#tau}", "", 50, 0., 10.);
	for (int i=0; i<50; ++i) {
		hmu->Fill(0.2 * i, best_fit_nu.mu[i]);
		he->Fill(0.2 * i, best_fit_nu.e[i]);
		htau->Fill(0.2 * i, best_fit_nu.tau[i]);
	}
	hmu->SetLineWidth(2);
	he->SetLineWidth(2);
	htau->SetLineWidth(2);
	hmu->SetLineColor(1);
	hmu->SetMaximum(1e4);
	hmu->SetMinimum(8e-2);
	hmu->Draw("HIST");
	he->SetLineColor(2);	
	he->Draw("SAME HIST");
	htau->SetLineColor(8);
	htau->Draw("SAME HIST");

	// plot original spectra as well
	plot_particles(nu.mu, nu.e);
	plot_antiparticles(nu.antimu, nu.antie);

	c3->BuildLegend(0.7, 0.7, 0.9, 0.9);
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

