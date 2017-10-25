
void csv_to_root() {
	TTree* spectrum = new TTree("spectrum", "");
	spectrum->ReadFile("data/spectrum.csv", "mu/F:antimu:e:antie", ',');
	
	TFile* f = new TFile("data/spectrum.root", "RECREATE");
	spectrum->Write();
	f->Close();
}
