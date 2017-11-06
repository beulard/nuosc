
// global experiment baseline
double L = 1300.;

struct spectrum {
	hierarchy h;

	// predicted fluxes at the far detector
	float mu_f[50];
	float e_f[50];
	// tau muons can also appear at the FD
	float tau_f[50];
};

struct initial_spectrum {
	float mu[50];
	float antimu[50];
	float e[50];
	float antie[50];
};


// propagation function
// calculates probability of a neutrino flavor a to oscillate to flavor b,
// given an energy E, a hierarchy h and a constant baseline
double P(flavor a, flavor b, float E, hierarchy* h, bool anti) {
	double p = 0.;

	for(int i=0; i<3; ++i) {
		p += pow(abs(h->MNS[a*3 + i] * conj(h->MNS[b*3 + i])), 2);
	}

	if (anti) {
		for (int i=0; i<3; ++i) {
			for (int j=0; j<3; ++j) {
				if (j > i) {
					p += 2. * real( conj(h->MNS[a*3 + i]) * h->MNS[a*3 + j]
						  * h->MNS[b*3 + i] * conj(h->MNS[b*3 + j])
						  * exp(complex<double>(-2.i * 1.2668 * h->dm2_mat[i*3 + j] * L/E)));
				}
			}
		}
	}
	else {
		for(int i=0; i<3; ++i) {
			for(int j=0; j<3; ++j) {
				if(j > i) {
					p += 2. * real(h->MNS[a*3 + i] * std::conj(h->MNS[a*3 + j]) 
							* conj(h->MNS[b*3 + i]) * h->MNS[b*3 + j]
							* exp(complex<double>(-2.i * 1.2668 *  h->dm2_mat[i*3 + j] * L/E)));
				}
			}
		}
	}
	return p;
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
	h1->SetMinimum(3e-4);
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

// Propagate neutrinos for each flavor and for each energy and get the FD flux.
void propagate(const initial_spectrum* is, spectrum* s) {	
	hierarchy* h = &(s->h);
	for (int i=0; i<50; ++i) {
		float E = 0.2 * i;
		// avoid zero energy
		if (i == 0) E = 1e-12;

			
		// eg a muon neutrino will oscillate into electron neutrinos AND tau neutrinos
		// and leave a certain amount of muon neutrinos behind as well.
		//
		// Electron neutrinos will do the same thing so we must add the number of 
		// leftover muon neutrinos to the number of oscillated electron neutrinos
		// to get the final number of muon neutrinos
		

		//TODO fix normalizations
		//		add the antineutrino fluxes back
		//		
		//
		// Find proportion of leftover/oscillated muon neutrinos
		// the proportion of leftover MUON nus is
		double leftover_mu = is->mu[i] * P(f_m, f_m, E, h, false);
		// proportion of oscillated mu -> electron neutrinos
		double mu_e = is->mu[i] * P(f_m, f_e, E, h, false);
		// oscillated mu -> taus
		double mu_tau = is->mu[i] * P(f_m, f_t, E, h, false);

		// same for electron
		double leftover_e = is->e[i] * P(f_e, f_e, E, h, false);
		double e_mu = is->e[i] * P(f_e, f_m, E, h, false);
		double e_tau = is->e[i] * P(f_e, f_t, E, h, false);

		double leftover_antimu = is->antimu[i] * P(f_m, f_m, E, h, true);
		double antimu_antie = is->antimu[i] * P(f_m, f_e, E, h, true);
		double antimu_antitau = is->antimu[i] * P(f_m, f_t, E, h, true);

		double leftover_antie = is->antie[i] * P(f_e, f_e, E, h, true);
		double antie_antimu = is->antie[i] * P(f_e, f_m, E, h, true);
		double antie_antitau = is->antie[i] * P(f_e, f_t, E, h, true);
		
		// we actually want to add the antineutrino fluxes to their respective neutrino
		// fluxes because there is no way to distinguish between them at the detector
		// (no magnet) so
		// store the values
		s->mu_f[i] = leftover_mu + e_mu + leftover_antimu + antie_antimu;
		s->e_f[i] = leftover_e + mu_e + leftover_antie + antimu_antie;
		s->tau_f[i] = mu_tau + e_tau + antimu_antitau + antie_antitau;
	}
}

// Normalize the flux to the CDR predicted event rate (vol.2, p.27).
// Mode is 0 for neutrino mode, 1 for antineutrino mode
void normalize(initial_spectrum* s, bool mode) {
	float flux_mu, flux_e, flux_antimu, flux_antie = 0;

	for (int i=0; i<50; ++i) {
		// 0.2 is dE
		flux_mu += 0.2 * s->mu[i];
		flux_e += 0.2 * s->e[i];
		flux_antimu += 0.2 * s->antimu[i];
		flux_antie += 0.2 * s->antie[i];
	}
	
	// values CDR vol. 2 p. 27 (neutrino mode) for 150 kt MW year exposure
	float mu_signal, antimu_signal, e_signal, antie_signal;
	if (!mode) {
		mu_signal = 10842.;
		antimu_signal = 958.;
		e_signal = 861.;
		antie_signal = 13.;
	} else {
		mu_signal = 2598.;
		antimu_signal = 3754.;
		e_signal = 61.;
		antie_signal = 167.;
	}
	

	for (int i=0; i<50; ++i) {
		s->mu[i] *= 0.2 *  mu_signal / flux_mu;
		s->e[i] *= 0.2 * e_signal / flux_e;
		s->antimu[i] *= 0.2 * antimu_signal / flux_antimu;
		s->antie[i] *= 0.2 * antie_signal / flux_antie;
	}
}

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
	return chisq(data->mu_f, s.mu_f, 50) + chisq(data->e_f, s.e_f, 50);
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
	spectrum best_fit_antinu;
	populate(&best_fit_nu.h, NH);
	populate(&best_fit_antinu.h, NH);

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

	
	// Normalize each spectrum to its integral (total flux)
	normalize(&nu, false);
	normalize(&antinu, true);


	// so now we want to repeat this for the different hierarchies, and different 
	// values of delta CP, and then calculate the chi squared between
	// those and the original propagation with physical delta and normal h
	propagate(&nu, &best_fit_nu);
	propagate(&antinu, &best_fit_antinu);

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
	float c2_antinu_nh[N];
	float c2_nu_ih[N];
	float c2_antinu_ih[N];


	for (int i=0; i<N; ++i) {
		d_cp[i] = (float)i / (float)N * 2. * TMath::Pi();
		hierarchy nh;
		populate(&nh, NH, d_cp[i]);
		c2_nu_nh[i] = spectrum_chisq(&nh, &nu, &best_fit_nu);
		c2_antinu_nh[i] = spectrum_chisq(&nh, &antinu, &best_fit_antinu);

		hierarchy ih;
		populate(&ih, IH, d_cp[i]);
		c2_nu_ih[i] = spectrum_chisq(&ih, &nu, &best_fit_nu);
		c2_antinu_ih[i] = spectrum_chisq(&ih, &antinu, &best_fit_antinu);
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

	float c2_antinu_nh_0 = spectrum_chisq(&nh_0, &antinu, &best_fit_antinu);
	float c2_antinu_nh_pi = spectrum_chisq(&nh_pi, &antinu, &best_fit_antinu);
	float c2_antinu_ih_0 = spectrum_chisq(&ih_0, &antinu, &best_fit_antinu);
	float c2_antinu_ih_pi = spectrum_chisq(&ih_pi, &antinu, &best_fit_antinu);


	float delta_c2_cpv_nh_nu[N];
	float delta_c2_cpv_nh_antinu[N];
	float delta_c2_cpv_ih_nu[N];
	float delta_c2_cpv_ih_antinu[N];

	float delta_c2_mh_nu[N];
	float delta_c2_mh_antinu[N];

	for (int i=0; i<N; ++i) {
		//Printf("%f %f", c2_nu_test_0 - c2_nu_nh[i], c2_nu_test_pi - c2_nu_nh[i]);
		delta_c2_cpv_nh_nu[i] = abs(min(c2_nu_nh_0 - c2_nu_nh[i],
				    		    	c2_nu_nh_pi - c2_nu_nh[i]));
		delta_c2_cpv_ih_nu[i] = abs(min(c2_nu_ih_0 - c2_nu_ih[i],
									c2_nu_ih_pi - c2_nu_ih[i]));
		delta_c2_cpv_nh_antinu[i] = abs(min(c2_antinu_nh_0 - c2_antinu_nh[i],
										    c2_antinu_nh_pi - c2_antinu_nh[i]));
		delta_c2_cpv_ih_antinu[i] = abs(min(c2_antinu_ih_0 - c2_antinu_ih[i],
										    c2_antinu_ih_pi - c2_antinu_ih[i]));

		delta_c2_mh_nu[i] = c2_nu_ih[i] - c2_nu_nh[i];
		delta_c2_mh_antinu[i] = c2_antinu_ih[i] - c2_antinu_nh[i];
		Printf("%f", delta_c2_mh_antinu[i]);
	}

	// This canvas contains the plots for the neutrino mode: delta_CPV and delta_MH
	TCanvas* c4 = new TCanvas("Neutrino mode", "", 1600, 1000);
	c4->Divide(2, 2);
	TPad* pad = (TPad*)c4->GetPad(1);
	pad->cd();
	TH1* hdc2_nh = new TH1F("normal hierarchy", "", N, 0., 2. * TMath::Pi());
	TH1* hdc2_ih = new TH1F("inverted hierarchy", "", N, 0., 2. * TMath::Pi());
	for (int i=0; i<N; ++i) {
		hdc2_nh->Fill((float)i / (float)N * 2. * TMath::Pi() + 1e-2, delta_c2_cpv_nh_nu[i]);
		hdc2_ih->Fill((float)i / (float)N * 2. * TMath::Pi() + 1e-2, delta_c2_cpv_ih_nu[i]);
	}
	hdc2_nh->Draw("HIST");
	hdc2_ih->SetLineColor(2);
	hdc2_ih->Draw("HIST SAME");
	// formatting
	pad->BuildLegend();
	pad->SetTicks();
	pad->SetGrid();
	hdc2_nh->GetXaxis()->SetTitle("#delta_{CP}");
	hdc2_nh->GetYaxis()->SetTitle("| #Delta#chi^{2}_{CPV} |");
	hdc2_nh->SetTitle("#nu mode");
	///

	pad = (TPad*)c4->GetPad(3);
	pad->cd();
	TH1* hdc2_mh = new TH1F("hdc2_mh", "", N, 0., 2. * TMath::Pi());
	for (int i=0; i<N; ++i) {
		hdc2_mh->Fill((float)i / (float)N * 2. * TMath::Pi() + 1e-2, delta_c2_mh_nu[i]);
	}
	hdc2_mh->Draw("HIST");
	// formatting
	pad->SetTicks();
	pad->SetGrid();
	hdc2_mh->GetXaxis()->SetTitle("#delta_{CP}");
	hdc2_mh->GetYaxis()->SetTitle("#Delta#chi^{2}_{MH}");
	///

	pad = (TPad*)c4->GetPad(2);
	pad->cd();
	TH1* hdc2_nh_antinu = new TH1F("normal hierarchy ", "", N, 0., 2. * TMath::Pi());
	TH1* hdc2_ih_antinu = new TH1F("inverted hierarchy ", "", N, 0., 2. * TMath::Pi());
	for (int i=0; i<N; ++i) {
		hdc2_nh_antinu->Fill((float)i / (float)N * 2. * TMath::Pi() + 1e-2,
			   			  	delta_c2_cpv_nh_antinu[i]);
		hdc2_ih_antinu->Fill((float)i / (float)N * 2. * TMath::Pi() + 1e-2,
			   				delta_c2_cpv_ih_antinu[i]);
	}
	hdc2_nh_antinu->Draw("HIST");
	hdc2_ih_antinu->SetLineColor(2);
	hdc2_ih_antinu->Draw("HIST SAME");
	// formatting
	pad->BuildLegend();
	pad->SetTicks();
	pad->SetGrid();
	hdc2_nh_antinu->GetXaxis()->SetTitle("#delta_{CP}");
	hdc2_nh_antinu->GetYaxis()->SetTitle("| #Delta#chi^{2}_{CPV} |");
	hdc2_nh_antinu->SetTitle("#bar{#nu} mode");
	///


	pad = (TPad*)c4->GetPad(4);
	pad->cd();
	TH1* hdc2_mh_antinu = new TH1F("hdc2_mh_anti", "", N, 0., 2. * TMath::Pi());
	for (int i=0; i<N; ++i) {
		hdc2_mh_antinu->Fill((float)i / (float)N * 2. * TMath::Pi() + 1e-2, 
							delta_c2_mh_antinu[i]);
	}
	hdc2_mh_antinu->Draw("HIST");
	// formatting
	pad->SetTicks();
	pad->SetGrid();
	hdc2_mh_antinu->GetXaxis()->SetTitle("#delta_{CP}");
	hdc2_mh_antinu->GetYaxis()->SetTitle("#Delta#chi^{2}_{MH}");
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
		hmu->Fill(0.2 * i, best_fit_nu.mu_f[i]);
		he->Fill(0.2 * i, best_fit_nu.e_f[i]);
		htau->Fill(0.2 * i, best_fit_nu.tau_f[i]);
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
