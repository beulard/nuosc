#include "spectrum.h"

void plot_dc2(int N, double* d_cp, double* dc2_mh_n, double* dc2_mh_i, 
					 double* dc2_cp_n, double* dc2_cp_i);


double mean_dc2(const double* test, const double* tru, int N) {
	double r=0;

	for (int i=0; i<N; ++i) {
		if (i>4)
		r += pow(test[i] - tru[i], 2) / test[i];
	}

	return r;
}


void hyperk() {

	// HyperK Baseline
	const double L = 295;

	initial_spectrum is;
	is.read("data/hyperk.csv");

	//is.plot();
	
	// For HyperK, the electron neutrino reconstructed energy falls between 0.1 and 1.3 GeV
	// Thus we must restrict our attention to this range when oscillating the spectrum and
	// calculating dc2.

	
	// Number of d_cp samples
	const int N = 151;
	// Values of d_cp
	double d_cp[N];

	// Delta chi squared for the mass hierarchy assuming true normal hierarchy
	double dc2_mh_n[N] = {0};
	double dc2_mh_i[N] = {0};
	// Delta chi squared for delta_CP
	double dc2_cp_n[N] = {0};
	double dc2_cp_i[N] = {0};

	// Parameter combinations for normalization
	parameters nh0, nhpi;
	nh0.populate(NH, 0.);
	nhpi.populate(NH, pi);

	spectrum nh0s, nhpis;

	nh0s.p = &nh0;
	nhpis.p = &nhpi;

	double E1 = 0.1;
	double E2 = 1.2;
	int N_bins = round((E2 - E1) / 0.05);
	Printf("%f", (E2 - E1) / 0.05);
	Printf("%d", N_bins);

	// Oscillate these spectra with specific parameters
	nh0s.oscillate(&is, L, E1, E2, N_bins);
	nhpis.oscillate(&is, L, E1, E2, N_bins);

	// Predicted event rate, from HyperK design report p.181
	const double e_nh = 2300;
	double nh0_int = nh0s.get_integral();
	Printf("%f", nh0_int);
	//Printf("%f", e_nh);

	// Integral for d_cp = 0 that we use to normalize spectra under NH assumption
	nh0s.normalize(nh0_int, e_nh);
	nhpis.normalize(nh0_int, e_nh);

	// IH Parameters combinations for normalization and d_cp sensitivity
	parameters ih0, ihpi;
	ih0.populate(IH, 0.);
	ihpi.populate(IH, pi);

	spectrum ih0s, ihpis;
	ih0s.p = &ih0;
	ihpis.p = &ihpi;
	ih0s.oscillate(&is, L, E1, E2, N_bins);
	ihpis.oscillate(&is, L, E1, E2, N_bins);

	double ih0_int = ih0s.get_integral();
	// For HyperK, we don't have predicted event rate for IH. We must estimate it ourselves.
	// It can be done by taking the ratio of probabilities of oscillation.
	// We want to evaluate the integral over the energy range of the probability under
	// NH and IH, then get the ratio IH/NH and multiply by the NH event rate.
	const int samples = 200;
	double x[samples];
	double int_nh = 0;
	double int_ih = 0;

	for (int i=0; i<samples; ++i) {
		// The energy range of interest is [0.1, 1.3] GeV for HyperK
		double delta_E = (E2 - E1) / samples;
		x[i] = i * delta_E + E1;

		int_nh += P_me(f_m, f_e, x[i], L, &nh0, false) * delta_E;
		int_ih += P_me(f_m, f_e, x[i], L, &ih0, false) * delta_E;
	}
	const double e_ih = e_nh * int_ih / int_nh;
	//const double e_ih = 2300;
	Printf("Norms: %f %f", e_nh, e_ih);
	

	ih0s.normalize(ih0_int, e_ih);
	ihpis.normalize(ih0_int, e_ih);


	parameters nh[N], ih[N];

	spectrum osc_nhs[N];
	spectrum osc_ihs[N];

	for (int i=0; i<N; ++i) {
		Printf("%d/%d", i+1, N);
		d_cp[i] = ((double)i / (double)(N-1) * 2. - 1.) * pi;

		nh[i].populate(NH, d_cp[i]);
		ih[i].populate(IH, d_cp[i]);

		osc_nhs[i].p = &nh[i];
		osc_ihs[i].p = &ih[i];
		osc_nhs[i].oscillate(&is, L, E1, E2, N_bins);
		osc_ihs[i].oscillate(&is, L, E1, E2, N_bins);

		osc_nhs[i].normalize(nh0_int, e_nh);
		osc_ihs[i].normalize(ih0_int, e_ih);
	}

	// Calculate delta chi squared
	for (int i=0; i<N; ++i) {
		
		// For each d_cp we want to calculate the minimum corresponding delta chi squared
		double min_dc2_n = 1e28;
		double min_dc2_i = 1e28;

		
		// So we perform another loop over delta_CP and take the minimum chi squared we find
		for (int j=0; j<N; ++j) {

			// When we assume normal hierarchy, the NH spectrum is fixed in i and we
			// search the IH with j
			double dc2_n = mean_dc2(osc_ihs[j].e,
									osc_nhs[i].e, N_bins);

			// Vice versa
			double dc2_i = mean_dc2(osc_nhs[j].e,
									osc_ihs[i].e, N_bins);

			

			min_dc2_n = min(dc2_n, min_dc2_n);
			min_dc2_i = min(dc2_i, min_dc2_i);

		}

		dc2_mh_n[i] = sqrt(min_dc2_n);
		dc2_mh_i[i] = sqrt(min_dc2_i);

		// For the CP sensitivity we want the j=0 elements (the ones with
		// the correct hierarchy)
		dc2_cp_n[i] = sqrt(min(mean_dc2(osc_nhs[i].e,
										osc_nhs[N/2].e, N_bins),
							   mean_dc2(osc_nhs[i].e,
								   	    osc_nhs[0].e, N_bins)));

		dc2_cp_i[i] = sqrt(min(mean_dc2(osc_ihs[i].e,
										osc_ihs[N/2].e, N_bins),
							   mean_dc2(osc_ihs[i].e,
								   		osc_ihs[0].e, N_bins)));
	}

	// To plot a test spectrum
	
	TCanvas* ccc = new TCanvas();
	//ccc->SetLogy();
	TH1* hnh = new TH1F("hnh", "", N_bins, E1, E2);
	TH1* hih = new TH1F("hih", "", N_bins, E1, E2);
	int index = (int)((-0. / 2. + .5) * N);

	//Printf("%f", d_cp[index] / pi);
	for (int i=0; i<N_bins; ++i) {
		//hnh->Fill(0.1 + (1.3 - 0.1) * (float)i / 25 + 1e-3, osc_nhs[index].e[i]);
		hnh->Fill(E1 + (E2 - E1) * (float)i / N_bins + 1e-3, nh0s.e[i]);
		//hnh->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, ih0s.events[E_SIGNAL][i]);
		hih->Fill(E1 + (E2 - E1) * (float)i / N_bins + 1e-3, osc_ihs[index].e[i]);
		//hih->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, ih0s.events[E_SIGNAL][i]);
	}

	hnh->Draw("hist");
	//hnh->SetMaximum(110);
	hih->Draw("hist same");
	hnh->SetLineWidth(2);
	hih->SetLineWidth(2);
	hnh->SetLineColor(ci[CI_NH]);
	hih->SetLineColor(ci[CI_IH]);



	plot_dc2(N, d_cp, dc2_mh_n, dc2_mh_i, dc2_cp_n, dc2_cp_i);
	
}


void plot_dc2(int N, double* d_cp, double* dc2_mh_n, double* dc2_mh_i, 
					 double* dc2_cp_n, double* dc2_cp_i) {
	
	double x[N];
	double ymnh[N];
	double ymih[N];
	double ydnh[N];
	double ydih[N];
	for (int i=0; i<N; ++i) {
		x[i] = d_cp[i] / pi;
		ymnh[i] = dc2_mh_n[i];
		ymih[i] = dc2_mh_i[i];
		ydnh[i] = dc2_cp_n[i];
		ydih[i] = dc2_cp_i[i];
	}


	TCanvas* c4 = new TCanvas("c4", "", 1000, 400);
	c4->SetFillColor(ci[CI_BACKGROUND]);
	c4->Divide(2, 1);
	
	TPad* p = (TPad*)c4->GetPad(1);
	p->cd();
	TGraph* gdc2_mnh = new TGraph(N, x, ymnh);
	//TGraph* gdc2_mnh_e = new TGraphErrors(N, x, ymnh, NULL, sd_mh_n);
	//gdc2_mnh_e->SetFillColor(ci[CI_NH]);
	//gdc2_mnh_e->Draw("a3");
	gdc2_mnh->Draw();
	gdc2_mnh->SetTitle("True normal hierarchy");
	gdc2_mnh->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mnh->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_mnh->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mnh->SetMinimum(0);
	//gdc2_mnh->SetLineColor(4);
	gdc2_mnh->SetLineColor(ci[CI_NH]);
	gdc2_mnh->SetMarkerColor(4);
	gdc2_mnh->SetLineWidth(3);

	TGraph* gdc2_mih = new TGraph(N, x, ymih);
	//TGraph* gdc2_mih_e = new TGraphErrors(N, x, ymih, NULL, sd_mh_i);
	gdc2_mih->SetTitle("True inverted hierarchy");
	c4->SetTitle("MH sensitivity");
	gdc2_mih->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mih->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	//gdc2_mih->SetLineColor(kRed+1);
	gdc2_mih->SetLineColor(ci[CI_IH]);
	gdc2_mih->SetMarkerColor(2);
	//gdc2_mih_e->SetFillColor(ci[CI_IH]);
	//gdc2_mnh_e->SetMaximum(25);
	//gdc2_mih_e->SetMaximum(25);
	gdc2_mih->SetLineWidth(3);
 	p = (TPad*)c4->GetPad(2);	
	p->cd();
	gdc2_mih->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mih->SetMinimum(0);
	//gdc2_mih->SetMaximum(25);
	//gdc2_mih_e->Draw("a3");
	gdc2_mih->Draw();
	

	// Draw delta_CP mean delta chi squared
	TCanvas* c5 = new TCanvas("c5", "", 1000, 400);
	c5->SetFillColor(ci[CI_BACKGROUND]);
	c5->Divide(2, 1);
	c5->cd(1);
	TGraph* gdc2_cp = new TGraph(N, x, ydnh);
	//gdc2_cp->Draw();
	gdc2_cp->Draw();
	gdc2_cp->SetTitle("True normal hierarchy");
	//gdc2_cp->SetMaximum(10);
	gdc2_cp->SetMinimum(0);
	gdc2_cp->GetXaxis()->SetLimits(-1, 1);
	gdc2_cp->SetLineColor(ci[CI_NH]);
	gdc2_cp->SetLineWidth(3);
	gdc2_cp->SetMarkerColor(ci[CI_NH]);
	gdc2_cp->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_cp->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	
	c5->cd(2);
	TGraph* gdc2_cp_ih = new TGraph(N, x, ydih);
	//gdc2_cp->Draw();
	gdc2_cp_ih->Draw("");
	gdc2_cp_ih->SetTitle("True inverted hierarchy");
	gdc2_cp_ih->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_cp_ih->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	//gdc2_cp->SetMaximum(10);
	//gdc2_cp->SetMaximum(10);
	gdc2_cp_ih->SetMinimum(0);
	gdc2_cp_ih->SetLineColor(ci[CI_IH]);
	gdc2_cp_ih->SetLineWidth(3);
	gdc2_cp_ih->SetMarkerColor(ci[CI_IH]);
	gdc2_cp_ih->GetXaxis()->SetLimits(-1, 1);


}
