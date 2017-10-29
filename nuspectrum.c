
// global experiment baseline
double L = 1300.;

// propagation function
// calculates probability of a neutrino flavor a to oscillate to flavor b,
// given an energy E, a hierarchy h and a constant baseline
double P(flavor a, flavor b, float E, hierarchy* h) {
	double p = 0.;

	for(int i=0; i<3; ++i) {
		p += pow(abs(h->MNS[a*3 + i] * std::conj(h->MNS[b*3 + i])), 2);
	}

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

void plot_initial(float* mu, float* antimu, float* e, float* antie) {
	TCanvas* c1 = new TCanvas("c1", "", 700, 600);
	// set log y axis and no statistics legend
	c1->SetLogy();	
	c1->SetTicks();
	gStyle->SetOptStat(0);

	TH1* h1 = new TH1F("h1", "", 50, 0., 10.);
	TH1* h2 = new TH1F("h2", "", 50, 0., 10.);
	TH1* h3 = new TH1F("h3", "", 50, 0., 10.);
	TH1* h4 = new TH1F("h4", "", 50, 0., 10.);
	
	// fill histograms manually with values from the spectrum
	for(int i=0; i<50; ++i) {
		h1->Fill(i * 0.2, mu[i]);
		h2->Fill(i * 0.2, antimu[i]);
		h3->Fill(i * 0.2, e[i]);
		h4->Fill(i * 0.2, antie[i]);
	}

	
	// plot the initial spectrum plots
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

void plot_normalized_particles(float* mu, float* e) {
	TH1* h1 = new TH1F("hp1", "", 50, 0., 10.);
	TH1* h2 = new TH1F("hp2", "", 50, 0., 10.);
	
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
	h1->SetMaximum(1.0);
	h1->SetMinimum(3e-4);
	h1->Draw("HIST");
	h2->SetLineColor(2);	
	h2->Draw("SAME HIST");
	
}

void read_spectrum(float* mu_vals, float* antimu_vals, float* e_vals, float* antie_vals) {
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
		mu_vals[i] = *mu;
		antimu_vals[i] = *antimu;
		e_vals[i] = *e;
		antie_vals[i] = *antie;	

		++i;
	}
}

void nuspectrum() {

	// fluxes at the accelerator
	float mu_i[50];
	float antimu_i[50];
	float e_i[50];
	float antie_i[50];

	// predicted fluxes at the far detector
	float mu_f[50];
	float antimu_f[50];
	float e_f[50];
	float antie_f[50];
	// tau muons can also appear at the FD
	float tau_f[50];
	float antitau_f[50];

	// read data from root file
	read_spectrum(mu_i, antimu_i, e_i, antie_i);

	gStyle->SetOptStat(0);
	// plot the initial spectrum
	//plot_initial(mu_i, antimu_i, e_i, antie_i);

	// we would like to normalize the fluxes so let's find the
	// maximum value and divide everything by that
	// the max value is located in the mu_i array (known from spectrum)
	float max_flux = 0.;
	for(int i=0; i<50; ++i) {
		if (mu_i[i] > max_flux)
			max_flux = mu_i[i];
	}
	Printf("maximum flux (norm): %e", max_flux);
	// now do the normalization
	for(int i=0; i<50; ++i) {
		mu_i[i] /= max_flux;
		antimu_i[i] /= max_flux;
		e_i[i] /= max_flux;
		antie_i[i] /= max_flux;
	}
	// plot the normalized spectrum
	//TCanvas* c = new TCanvas("c", "", 800, 600);
	//c->SetLogy();
	//c->SetTicks();
	//plot_normalized(mu_i, antimu_i, e_i, antie_i);
	

	// setup hierarchy
	hierarchy nh;
	populate(&nh, NH);

	// propagate neutrinos for each flavor and for each energy and get the
	// flux at the far detector
	for (int i=0; i<50; ++i) {
		float E = 0.2 * i;
		// avoid zero energy
		if (i == 0) E = 1e-12;

			
		// eg a muon neutrino will oscillate into electron neutrinos (neglecting taus)
		// and leave a certain amount of muon neutrinos as well.
		// electron neutrinos will do the same thing so we must add the number of 
		// leftover muon neutrinos to the number of oscillated electron
		// neutrinos, and then multiply this by the normalized flux to get
		// the flux at the FD.
		
		// Find proportion of leftover/oscillated muon neutrinos
		// the proportion of leftover MUON nus is
		double leftover_mu = P(f_m, f_m, E, &nh);
		// proportion of oscillated ELECTRON neutrinos
		double oscillated_mu = P(f_m, f_e, E, &nh);

		// same for electron
		double leftover_e = P(f_e, f_e, E, &nh);
		double oscillated_e = P(f_e, f_m, E, &nh);

		// now we need to multiply by the normalized flux and store the values
		mu_f[i] = mu_i[i] * (leftover_mu + oscillated_e);
		e_f[i] = e_i[i] * (leftover_e + oscillated_mu);

		// but mu and e can also oscillate into tau, even if there are no taus
		// at the beam.
		// the flux is the sum of probabilities weighted by the respective 
		// initial fluxes
		tau_f[i] = mu_i[i] * P(f_m, f_t, E, &nh) + e_i[i] * P(f_e, f_t, E, &nh);
	}
	// now plot the final flux
	TCanvas* c3 = new TCanvas("c3", "", 1100, 800);
	//c3->Divide(1, 2);
	c3->SetLogy();
	c3->SetTicks();

	//TPad* p1 = (TPad*)c3->GetPad(1);
	//p1->cd();
	//p1->SetLogy();
	//p1->SetTicks();
	plot_normalized_particles(mu_i, e_i);

	//TPad* p2 = (TPad*)c3->GetPad(2);
	//p2->cd();
	//p2->SetLogy();
	//p2->SetTicks();
	
	TH1* hmu = new TH1F("hmu", "", 50, 0., 10.);
	TH1* he = new TH1F("he", "", 50, 0., 10.);
	TH1* htau = new TH1F("htau", "", 50, 0., 10.);
	for (int i=0; i<50; ++i) {
		hmu->Fill(0.2 * i, mu_f[i]);
		he->Fill(0.2 * i, e_f[i]);
		htau->Fill(0.2 * i, tau_f[i]);
	}
	hmu->SetLineWidth(2);
	he->SetLineWidth(2);
	htau->SetLineWidth(2);

	hmu->SetLineColor(1);
	hmu->SetMaximum(1.0);
	hmu->SetMinimum(3e-4);
	hmu->Draw("HIST SAME");
	he->SetLineColor(2);	
	he->Draw("SAME HIST");
	htau->SetLineColor(8);
	htau->Draw("SAME HIST");

}
