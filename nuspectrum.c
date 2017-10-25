

void nuspectrum() {
	TFile* f = TFile::Open("data/spectrum.root");

	TH1* h1 = new TH1F("h1", "", 50, 0., 10.);
	TH1* h2 = new TH1F("h2", "", 50, 0., 10.);
	TH1* h3 = new TH1F("h3", "", 50, 0., 10.);
	TH1* h4 = new TH1F("h4", "", 50, 0., 10.);
	
	float mu_vals[50];
	float antimu_vals[50];
	float e_vals[50];
	float antie_vals[50];

	TTree* t = (TTree*)f->Get("spectrum");
	t->Print();

	TTreeReader r("spectrum", f);
	TTreeReaderValue<float> mu(r, "mu");
	TTreeReaderValue<float> antimu(r, "antimu");
	TTreeReaderValue<float> e(r, "e");
	TTreeReaderValue<float> antie(r, "antie");

	int i = 0;
	while(r.Next()) {
		h1->Fill(i * 0.2, *mu);
		h2->Fill(i * 0.2, *antimu);
		h3->Fill(i * 0.2, *e);
		h4->Fill(i * 0.2, *antie);

		mu_vals[i] = *mu;
		antimu_vals[i] = *antimu;
		e_vals[i] = *e;
		antie_vals[i] = *antie;	

		++i;
	}

	gStyle->SetOptStat(0);
	
	//	initial spectrum plots
	TCanvas* c1 = new TCanvas("c1", "", 700, 600);
	c1->SetLogy();	
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
