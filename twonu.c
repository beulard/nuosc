// delta m32^2 (eV^2)
const double dm2 = 2.50e-3;
// mixing = sin^2(2 * theta23) (delta m^2 < 0)
const double mixing = 0.90;
 
double sinsq(double x) {
	return TMath::Power(TMath::Sin(x), 2);
}

// transition probability
double Pt(double x, double mix, double dm2) {
	return mix * sinsq(1.27 * dm2 * x);
}



void twonu() {
	TCanvas* canvas = new TCanvas("c1", "NuOsc2");
	TF1* f1 = new TF1("f1", "Pt(x, [0], [1])", 0., 4000.);
	TF1* f2 = new TF1("f2", "1.0 - Pt(x, [0], [1])", 0., 4000.);

	f1->SetParameter(0, mixing);
	f2->SetParameter(0, mixing);
	f1->SetParameter(1, dm2);
	f2->SetParameter(1, dm2);

	f2->SetTitle("");
	f2->SetMaximum(1.0);
	f2->SetMinimum(0.0);
	f2->SetLineColor(4);
	f2->GetXaxis()->SetTitle("L/E (km/GeV)");
	f2->GetYaxis()->SetTitle("P");
	f2->Draw();
	f1->Draw("same");
}


