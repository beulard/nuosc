
void csv_to_root() {
	TNtuple* mu_data = new TNtuple("spectrum", "", "f");
	TNtuple* antimu_data = new TNtuple("spectrum", "", "f");
	TNtuple* e_data = new TNtuple("spectrum", "", "f");
	TNtuple* antie_data = new TNtuple("spectrum", "", "f");
	
	mu_data->ReadFile("data/mu.csv");
	antimu_data->ReadFile("data/antimu.csv");
	e_data->ReadFile("data/e.csv");
	antie_data->ReadFile("data/antie.csv");
	
	TFile* f = new TFile("data/mu.root", "RECREATE");
	mu_data->Write();
	f->Close();
	f = new TFile("data/antimu.root", "RECREATE");
	antimu_data->Write();
	f->Close();
	f = new TFile("data/e.root", "RECREATE");
	e_data->Write();
	f->Close();
	f = new TFile("data/antie.root", "RECREATE");
	antie_data->Write();

	f->Close();
}
