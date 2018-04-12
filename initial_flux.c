struct initial {
	double mu[50];	
	double e[50];	
	double antimu[50];	
	double antie[50];	
};

void read_spectrum(initial* s) {
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

void plot_initial(initial* s) {
	TCanvas* c1 = new TCanvas("c1", "", 700, 600);
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
		h1->Fill(i * 0.2, s->mu[i]);
		h2->Fill(i * 0.2, s->antimu[i]);
		h3->Fill(i * 0.2, s->e[i]);
		h4->Fill(i * 0.2, s->antie[i]);
	}

	
	// plot the initial spectrum plots
	h1->SetLineWidth(3);
	h2->SetLineWidth(3);
	h3->SetLineWidth(3);
	h4->SetLineWidth(3);

	h1->GetXaxis()->SetTitleSize(0.045);
	h1->GetXaxis()->SetLabelSize(0.039);
	h1->GetYaxis()->SetTitleSize(0.045);
	h1->GetYaxis()->SetLabelSize(0.039);

	h1->GetXaxis()->SetTitle("Energy / GeV");
	h1->GetYaxis()->SetTitle("\#nu flux/m^{2}/GeV/10^{20} POT at 1300 km");

	h1->SetLineColor(ci[CI_2]);
	h1->SetMinimum(1e6);
	h1->SetMaximum(1e10);
	h1->Draw("HIST");
	h2->SetLineColor(ci[CI_1]);	
	h2->Draw("SAME HIST");
	h3->SetLineColor(ci[CI_3]);
	h3->Draw("SAME HIST");
	h4->SetLineColor(ci[CI_4]);
	h4->Draw("SAME HIST");

	gPad->SetGrid();
	gPad->BuildLegend();

}

void initial_flux() {
	initial dune;
	read_spectrum(&dune);

	plot_initial(&dune);
}
