

void nuspectrum() {
	TFile* f = TFile::Open("data/mu.root");
	if(f->IsZombie()) {
		Printf("caca");
	}
	TNtuple* mutree = (TNtuple*)f->Get("spectrum");	

	TH1* h1 = new TH1F("h1", "", 50, 0., 10.);
	

	TTreeReader r("spectrum", f);
	TTreeReaderValue<float> mu(r, "f");
	int i = 0;
	while(r.Next()) {
		h1->Fill(i * 0.2, *mu);
		++i;
	}

	gStyle->SetOptStat(0);
	
	TCanvas* c1 = new TCanvas("c1", "", 700, 600);
	c1->SetLogy();	
	h1->SetLineColor(1);
	h1->SetMinimum(1e6);
	h1->SetMaximum(1e10);
	h1->Draw("HIST");
	
	//	TODO maybe arrange data in one root file if possible
	//	draw other curves
	//	start propagating neutrinos
	
}
