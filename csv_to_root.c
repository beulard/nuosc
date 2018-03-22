
// Converts a CSV file to ROOT format
void csv_to_root(const char* file) {
	TTree* spectrum = new TTree("spectrum", "");
	spectrum->ReadFile("data/spectrum.csv", "mu/D:antimu:e:antie", ',');
	
	TFile* f = new TFile("data/spectrum.root", "RECREATE");
	spectrum->Write();
	f->Close();
	delete spectrum;
}
