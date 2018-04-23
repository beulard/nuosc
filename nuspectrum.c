#include "spectrum.h"
#include "helper.h"

// plotting and helper functions
void plot_initial(float* mu, float* antimu, float* e, float* antie);
void plot_normalized(float* mu, float* antimu, float* e, float* antie);
void plot_particles(float* mu, float* e);
void plot_antiparticles(float* antimu, float* antie);
void read_spectrum(initial_spectrum* s);
void plot_dc2(double* d_cp, double* dc2_mh_n, double* sd_mh_n, double* dc2_mh_i, double* sd_mh_i, double* dc2_cp, double* sd_cp_n, double* dc2_cp_ih, double* sd_cp_i, int N);
void do_dc2(const initial_spectrum* s, double baseline);
void do_dc2_theta(const initial_spectrum* s, double baseline);


// Mean delta chi-squared between a test spectrum and a "true" spectrum and a set of 
// standard deviations,
// as given in 1210.3651 p.9.
double mean_dc2(const double* test, const double* tru, const double* sd = NULL, int N = 37) {
	double r = 0.;
	//if (sd == NULL) {
		double avg_rate = 0.;
		for (int i=0; i<N; ++i) {
			avg_rate += test[i] / N;
		}
		for (int i=0; i<N; ++i) {
			//r += pow(test[i] - tru[i], 2) / /*tru[i]*/ avg_rate;
			if (i>6 && tru[i] > 1)
			r += pow(test[i] - tru[i], 2) / test[i];//2./avg_rate;// / pow(test[i] / 5., 2);
		}
	/*} else {
		for (int i=0; i<N; ++i) {
			if (i > 4)
			r += pow(test[i] - tru[i], 2) / pow(sd[i], 2);
		}
	}*/
	return r;
}



// theta=0: plot sensitivity with error band
// theta=1: plot sensitivity for best fit and 1sigma theta23 as band
void nuspectrum(int theta=0) {
	gStyle->SetOptStat(0);


	initial_spectrum is_dune;
	is_dune.read("data/dune_mu.csv");

	//is_dune.plot();

	initial_spectrum is_hyperk;
	is_hyperk.read("data/hyperk.csv");

	//is_hyperk.plot();

	//initial_spectrum is_spl;
	//is_spl.read("data/southpole_mu.csv");

	//is_spl.plot();



	// switch particle/antiparticle spectrum around for the antinu mode
	/*for (int i=0; i<50; ++i) {
		antinu.mu[i] = nu.antimu[i];
		antinu.antimu[i] = nu.mu[i];
		antinu.e[i] = nu.antie[i];
		antinu.antie[i] = nu.e[i];
	}*/
	
	if(theta == 0)
		do_dc2(&is_dune, 1300);
	else if(theta == 1)
		do_dc2_theta(&is_dune, 1300);
}



// Calculates the mean oscillation spectrum over a given set of N parameter combinations
// Also provides the standard deviation through 'out_sd'
// Output arrays should have size double[Nbins]
/*void get_MC_spectrum(parameters* p, int N, double* out_spectrum, double* out_sd) {
	double spectra[Nbins][N];
	for (int i=0; i<N; ++i) {
		for (int j=0; j<Nbins; ++j) {
			spectra[j][i] = P_me(f_m, f_e, (j + firstbin) * 0.2, &p[i], false);
		}
	}
	for (int i=0; i<Nbins; ++i) {
		out_spectrum[i] = mean(spectra[i], N);
		out_sd[i] = std_dev(spectra[i], out_spectrum[i], N);
	}
}*/

/*void plot_MC_spectrum(double* spectrum, double* sd) {

	double x[Nbins];
	for (int i=0; i<Nbins; ++i) {
		x[i] = (i + firstbin) * 0.2;
	}

	TCanvas* c = new TCanvas();
	TGraph* gs = new TGraph(Nbins, x, spectrum);
	TGraph* gsd = new TGraphErrors(Nbins, x, spectrum, NULL, sd);
	gsd->SetFillColor(3);
	gsd->Draw("a3");
	gs->Draw("same");

}*/


void do_dc2(const initial_spectrum* is, double L) {
	// Number of times we will evaluate the sensitivity to compose the final plot
	const int N_rec = 1;
	// Number of d_cp samples
	const int N = 31;
	// Values of d_cp
	double d_cp[N];

	// "Debug" variable: do we reconstruct spectra?
	const bool reconstruct = true;

	// Delta chi squared for the mass hierarchy assuming true normal hierarchy
	double dc2_mh_n[N][N_rec] = {0};
	double dc2_mh_i[N][N_rec] = {0};
	// Delta chi squared for delta_CP
	double dc2_cp_n[N][N_rec] = {0};
	double dc2_cp_i[N][N_rec] = {0};


	// Parameter combinations for normalization
	parameters nh0, nhpi;
	// Hierarchy that has NH parameters (theta23) but inverted mass hierarchy
	parameters nhi0;
	nh0.populate(NH, 0.);
	nhpi.populate(NH, pi);
	nhi0.populate(NH, 0.);
	nhi0.flip_hierarchy();

	spectrum nh0s, nhpis;
	spectrum nhi0s;

	nh0s.p = &nh0;
	nhpis.p = &nhpi;
	nhi0s.p = &nhi0;

	// Oscillate these spectra with specific parameters
	nh0s.oscillate(is, L);
	nhpis.oscillate(is, L);
	nhi0s.oscillate(is, L);

	const double dune_e_nh = 861;
	const double dune_e_ih = 495;

	// Integral for d_cp = 0 that we use to normalize spectra under NH assumption
	double nh0_int = nh0s.get_integral();
	nh0s.normalize(nh0_int, dune_e_nh);
	nhpis.normalize(nh0_int, dune_e_nh);

	double nhi0_int = nhi0s.get_integral();
	// We normalize nhi0s to the IH event rates because it is oscillated under IH
	nhi0s.normalize(nhi0_int, IH);

	if (reconstruct) { 
	nh0s.reconstruct(500);
	nhpis.reconstruct(500);
	//nhi0s.reconstruct(100);
	}


	parameters ih0, ihpi;
	// Hierarchy that has IH parameters (theta23) but NH mass hierarchy
	parameters ihn0;
	ih0.populate(IH, 0.);
	ihpi.populate(IH, pi);
	ihn0.populate(IH, 0.);
	ihn0.flip_hierarchy();

	spectrum ih0s, ihpis;
	spectrum ihn0s;
	ih0s.p = &ih0;
	ihpis.p = &ihpi;
	ihn0s.p = &ihn0;
	ih0s.oscillate(is, L);
	ihpis.oscillate(is, L);
	ihn0s.oscillate(is, L);

	double ih0_int = ih0s.get_integral();
	ih0s.normalize(ih0_int, dune_e_ih);
	ihpis.normalize(ih0_int, dune_e_ih);

	double ihn0_int = ihn0s.get_integral();
	ihn0s.normalize(ihn0_int, NH);

	if (reconstruct) {
	ih0s.reconstruct(500);
	ihpis.reconstruct(500);
	//ihn0s.reconstruct(100);
	}


	// N * 2 set of parameters
	// The nh array contains parameters that will be used in calculating
	// sensitivity under the assumption that the mass hierarchy is normally ordered.
	// Technically, the second element (nh[i][1]) has inverted ordering.
	parameters nh[N][2], ih[N][2];

	// Oscillated and reconstructed (multiple times) spectra
	spectrum osc_nhs[N][2];
	spectrum osc_ihs[N][2];

	// Spectra that will be reconstructed only once
	spectrum *nhs, *ihs;
	nhs = (spectrum*)malloc(N * N_rec * sizeof(spectrum) * 2);
	ihs = (spectrum*)malloc(N * N_rec * sizeof(spectrum) * 2);
	

	// We store the event rate uncertainties here
	double osc_nhs_sd[N][2];
	double osc_ihs_sd[N][2];

	// Number of Monte Carlo iterations for each parameter
	const int N_MC = 20;


	const int Nbins=37;
	// TODO REMOVE
	// Average event rates.
	// The last word (eg t23) denotes the parameter we will vary in a Monte Carlo
	// simulation.
	double mean_rate_nh_t23[N][Nbins];
	double mean_rate_ih_t23[N][Nbins];
	double mean_rate_t13[N][Nbins];
	memset(mean_rate_nh_t23, 0, N * Nbins * sizeof(double));
	memset(mean_rate_ih_t23, 0, N * Nbins * sizeof(double));

	// Standard deviations around the means
	double sd_nh_t23[N][Nbins];
	double sd_ih_t23[N][Nbins];
	memset(sd_nh_t23, 0, N * Nbins * sizeof(double));
	memset(sd_ih_t23, 0, N * Nbins * sizeof(double));


	// Initialize the hierarchies
	for (int i=0; i<N; ++i) {
		//Printf("%d/%d", i+1, N);
		d_cp[i] = ((float)i / (float)(N-1) * 2. - 1.) * pi;

		// We can either populate with NH best fit values and then flip the sign of dm2_31
		// or we can populate with respective best fit values (latter makes more sense)

		/* //UNDER TESTING
		for (int j=0; j<2; ++j) {
			nh[i][j].populate(NH, d_cp[i]);
			ih[i][j].populate(IH, d_cp[i]);
		}
		nh[i][1].flip_hierarchy();
		ih[i][1].flip_hierarchy();
		*/
		nh[i][0].populate(NH, d_cp[i]);
		nh[i][1].populate(IH, d_cp[i]);
		ih[i][0].populate(IH, d_cp[i]);
		ih[i][1].populate(NH, d_cp[i]);



		for (int j=0; j<2; ++j) {
			osc_nhs[i][j].p = &nh[i][j];
			osc_ihs[i][j].p = &ih[i][j];
			osc_nhs[i][j].oscillate(is, L);
			osc_ihs[i][j].oscillate(is, L);

			// Normalize NH spectra to NH0 and IH spectra to IH0
			if (j==0) {
				osc_nhs[i][j].normalize(nh0_int, dune_e_nh);
				osc_ihs[i][j].normalize(ih0_int, dune_e_ih);
			}
			/*else {
				osc_nhs[i][j].normalize(nhi0_int, IH);
				osc_ihs[i][j].normalize(ihn0_int, NH);
			}*/
			else {
				osc_nhs[i][j].normalize(ih0_int, IH);
				osc_ihs[i][j].normalize(nh0_int, NH);
			}

			if (reconstruct) {
			osc_nhs[i][j].reconstruct(500);
			osc_ihs[i][j].reconstruct(500);
			}

			

			// Oscillate/reconstruct the spectra
			// It is important that we get ALL the spectra before we go
			// into the next loop because it will go through them all to calculate
			// sensitivity.
			for (int n=0; n<N_rec; ++n) {
				nhs[n*N*2 + i*2 + j].p = &nh[i][j];
				ihs[n*N*2 + i*2 + j].p = &ih[i][j];
				//osc_nhs[i][j].reconstruct(&nhs[n*N*2 + i*2 + j], 1);
				//osc_ihs[i][j].reconstruct(&ihs[n*N*2 + i*2 + j], 1);
			}
		}

		
	}



	Printf("N_rec = %d", N_rec);
	for (int n=0; n<N_rec; ++n) {
		// Loop over delta_CP values
		for (int i=0; i<N; ++i) {

			
			// for each d_cp we want to calculate the minimum corresponding delta chi squared
			double min_dc2_n = 1e28;
			double min_dc2_i = 1e28;

			// So we perform another loop over delta_CP and take the minimum chi squared we find
			for (int j=0; j<N; ++j) {

				// When we assume normal hierarchy, the NH spectrum is fixed in i and we
				// search the IH with j
				double dc2_n = mean_dc2(//nhs[n*N*2 + j*2 + 1].events[E_SIGNAL], 
										osc_ihs[j][0].e,
										//nhs[n*N*2 + i*2 + 0].e);
										osc_nhs[i][0].e,
										sd_nh_t23[i]);

				// Vice versa
				double dc2_i = mean_dc2(//ihs[n*N*2 + j*2 + 1].e, 
										osc_nhs[j][0].e,
										//ihs[n*N*2 + i*2 + 0].e);
										osc_ihs[i][0].e, sd_ih_t23[i]);

				

				min_dc2_n = min(dc2_n, min_dc2_n);
				min_dc2_i = min(dc2_i, min_dc2_i);

			}

			dc2_mh_n[i][n] = sqrt(min_dc2_n);
			dc2_mh_i[i][n] = sqrt(min_dc2_i);

			// For the CP sensitivity we want the j=0 elements (the ones with
			// the correct hierarchy)
			dc2_cp_n[i][n] = sqrt(min(min(mean_dc2(osc_nhs[i][0].e,
										  osc_nhs[N/2][0].e),
				   	 	         mean_dc2(osc_nhs[i][0].e,
										  osc_nhs[0][0].e)),
								 mean_dc2(osc_nhs[i][0].e,
									 	  osc_nhs[N-1][0].e)));

			//dc2_cp_n[i][n] = sqrt(min(mean_dc2(osc_nhs[i][0].e,
			//							  osc_nhs[N/2][0].e),
			//	   	 	         mean_dc2(osc_nhs[i][0].e,
			//							  osc_nhs[N-1][0].e)));

			dc2_cp_i[i][n] = sqrt(min(min(mean_dc2(//ihs[n*N*2 + i*2 + 0].e, 
											   osc_ihs[i][0].e,
									      //ih0s.e),
										  osc_ihs[N/2][0].e),
							     mean_dc2(//ihs[n*N*2 + i*2 + 0].e, 
									      osc_ihs[i][0].e,
								   	      //ihpis.e)));
										  osc_ihs[0][0].e)),
							    	mean_dc2(osc_ihs[i][0].e,
											 osc_ihs[N-1][0].e)));

			//dc2_cp_i[i][n] = sqrt(min(mean_dc2(//ihs[n*N*2 + i*2 + 0].e, 
			//								   osc_ihs[i][0].e,
			//						      //ih0s.e),
			//							  osc_ihs[N/2][0].e),
			//				     mean_dc2(//ihs[n*N*2 + i*2 + 0].e, 
			//						      osc_ihs[i][0].e,
			//					   	      //ihpis.e)));
			//							  osc_ihs[N-1][0].e)));
		}

	}


	// Mean sensitivities
	double mean_dc2_mh_n[N];
	double mean_dc2_mh_i[N];
	double mean_dc2_cp_n[N];
	double mean_dc2_cp_i[N];

	// Standard deviations from the mean
	double sd_mh_n[N];
	double sd_mh_i[N];
	double sd_cp_n[N];
	double sd_cp_i[N];

	// After calculating all sensitivities, we want to create an average and a
	// standard deviation for each sample d_cp point
	for (int i=0; i<N; ++i) {
		double dc2_max = 0.;
		double dc2_min = 10000.;
		for (int j=0; j<N_rec; ++j) {
			dc2_min = min(dc2_min, dc2_mh_n[i][j]);
			dc2_max = max(dc2_max, dc2_mh_n[i][j]);
		}
		//mean_dc2_mh_n[i] = mean(dc2_mh_n[i], N_rec);
		mean_dc2_mh_n[i] = (dc2_max + dc2_min) / 2.;
		mean_dc2_mh_i[i] = mean(dc2_mh_i[i], N_rec);
		mean_dc2_cp_n[i] = mean(dc2_cp_n[i], N_rec);
		mean_dc2_cp_i[i] = mean(dc2_cp_i[i], N_rec);

		sd_mh_n[i] = 0;//(dc2_max - dc2_min) / 2.;
		//sd_mh_n[i] = sqrt(std_dev(dc2_mh_n[i], mean_dc2_mh_n[i], N_rec));
		sd_mh_i[i] = 0;//std_dev(dc2_mh_i[i], mean_dc2_mh_i[i], N_rec);
		sd_cp_n[i] = 0;//std_dev(dc2_cp_n[i], mean_dc2_cp_n[i], N_rec);
		sd_cp_i[i] = 0;//std_dev(dc2_cp_i[i], mean_dc2_cp_i[i], N_rec);
	}

	// To plot a test spectrum
	
	TCanvas* ccc = new TCanvas();
	//ccc->SetLogy();
	TH1* hnh = new TH1F("NH", "", Nbins, 0.6, 8.);
	TH1* hih = new TH1F("IH", "", Nbins, 0.6, 8.);
	int index = (int)((+.5 / 2. + .5) * N);
	Printf("d_cp = %f", d_cp[index]);
	Printf("d_cp(N/2) = %f", d_cp[N/2]);
	Printf("d_cp(N-1) = %f", d_cp[N-1]);

	//Printf("%f", d_cp[index] / pi);
	for (int i=0; i<Nbins; ++i) {
		hnh->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, osc_nhs[index][0].e[i]);
		//hnh->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, nh0s.e[i]);
		hih->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, osc_ihs[index][0].e[i]);
		//hih->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, ih0s.events[E_SIGNAL][i]);
	}

	hnh->Draw("hist");
	hnh->SetMaximum(90); hih->Draw("hist same"); hnh->SetLineWidth(3);
	hih->SetLineWidth(3);
	hnh->SetLineColor(ci[CI_NH]);
	hih->SetLineColor(ci[CI_IH]);
	hnh->GetYaxis()->SetLabelSize(0.045);
	hnh->GetXaxis()->SetLabelSize(0.045);

	hnh->GetYaxis()->SetTitle("Event rate");
	hnh->GetYaxis()->SetTitleSize(0.049);
	hnh->GetYaxis()->SetTitleOffset(0.9);
	hnh->GetXaxis()->SetTitle("E (GeV)");
	hnh->GetXaxis()->SetTitleSize(0.049);
	
	ccc->BuildLegend();

	

	
	plot_dc2(d_cp, mean_dc2_mh_n, sd_mh_n, mean_dc2_mh_i, sd_mh_i, 
			       mean_dc2_cp_n, sd_cp_n,  mean_dc2_cp_i, sd_cp_i, N);
	

	free(nhs);
	free(ihs);
}

void plot_dc2(double* d_cp, double* dc2_mh_n, double* sd_mh_n, double* dc2_mh_i, double* sd_mh_i, double* dc2_cp, double* sd_cp_n, double* dc2_cp_ih, double* sd_cp_i, int N) {
	double x[N];
	double ymnh[N];
	double ymih[N];
	double ydnh[N];
	double ydih[N];
	for (int i=0; i<N; ++i) {
		x[i] = d_cp[i] / pi;
		ymnh[i] = dc2_mh_n[i];
		ymih[i] = dc2_mh_i[i];
		ydnh[i] = dc2_cp[i];
		ydih[i] = dc2_cp_ih[i];
	}

	//srand(time(NULL));
	//int r = rand();
	//char title1[64];
	//char title2[64];
	//snprintf(title1, 64, "sens_mh %d", r);
	//snprintf(title2, 64, "sens_cp %d", r);

	//TCanvas* c4 = new TCanvas(title1, "", 1000, 400);
	//c4->SetFillColor(ci[CI_BACKGROUND]);
	//c4->Divide(2, 1);
	
	//TPad* p = (TPad*)c4->GetPad(1);
	//p->cd();
	TCanvas* c4 = new TCanvas();
	TGraph* gdc2_mnh = new TGraph(N, x, ymnh);
	//TGraph* gdc2_mnh_e = new TGraphErrors(N, x, ymnh, NULL, sd_mh_n);

	gdc2_mnh->GetXaxis()->SetTitleSize(0.046);
	gdc2_mnh->GetXaxis()->SetLabelSize(0.041);
	gdc2_mnh->GetYaxis()->SetTitleSize(0.045);
	gdc2_mnh->GetYaxis()->SetLabelSize(0.041);

	gdc2_mnh->SetFillColor(ci[CI_NH]);
	//gdc2_mnh_e->Draw("a3");
	gdc2_mnh->SetTitle("True normal hierarchy");
	gdc2_mnh->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mnh->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_mnh->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mnh->SetMinimum(0);
	//gdc2_mnh->SetLineColor(4);
	gdc2_mnh->SetLineColor(ci[CI_NH]);
	gdc2_mnh->SetMarkerColor(4);
	gdc2_mnh->SetLineWidth(3);
	gdc2_mnh->Draw();

	TCanvas* c5 = new TCanvas();
	TGraph* gdc2_mih = new TGraph(N, x, ymih);
	//TGraph* gdc2_mih_e = new TGraphErrors(N, x, ymih, NULL, sd_mh_i);

	gdc2_mih->GetXaxis()->SetTitleSize(0.046);
	gdc2_mih->GetXaxis()->SetLabelSize(0.041);
	gdc2_mih->GetYaxis()->SetTitleSize(0.049);
	gdc2_mih->GetYaxis()->SetLabelSize(0.041);

	gdc2_mih->SetTitle("True inverted hierarchy");
	//c4->SetTitle("MH sensitivity");
	gdc2_mih->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mih->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	//gdc2_mih->SetLineColor(kRed+1);
	gdc2_mih->SetLineColor(ci[CI_IH]);
	gdc2_mih->SetMarkerColor(2);
	gdc2_mih->SetFillColor(ci[CI_IH]);
	//gdc2_mnh_e->SetMaximum(25);
	//gdc2_mih_e->SetMaximum(25);
	gdc2_mih->SetLineWidth(3);
 	//p = (TPad*)c4->GetPad(2);	
	//p->cd();
	gdc2_mih->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mih->SetMinimum(0);
	//gdc2_mih->SetMaximum(25);
	//gdc2_mih_e->Draw("a3");
	gdc2_mih->Draw();
	

	// Draw delta_CP mean delta chi squared
	//TCanvas* c5 = new TCanvas(title2, "", 1000, 400);
	//c5->SetFillColor(ci[CI_BACKGROUND]);
	//c5->Divide(2, 1);
	//c5->cd(1);
	TCanvas* c6 = new TCanvas();
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
	
	gdc2_cp->GetXaxis()->SetTitleSize(0.046);
	gdc2_cp->GetXaxis()->SetLabelSize(0.041);
	gdc2_cp->GetYaxis()->SetTitleSize(0.049);
	gdc2_cp->GetYaxis()->SetLabelSize(0.041);
	gdc2_cp->SetMaximum(9.);
	
	//c5->cd(2);
	TCanvas* c7 = new TCanvas();
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

	gdc2_cp_ih->SetMaximum(9.);
	gdc2_cp_ih->GetXaxis()->SetTitleSize(0.046);
	gdc2_cp_ih->GetXaxis()->SetLabelSize(0.041);
	gdc2_cp_ih->GetYaxis()->SetTitleSize(0.049);
	gdc2_cp_ih->GetYaxis()->SetLabelSize(0.041);
}

// Calculates MH mean delta chi squared and standard deviation given a hierarchy and theta_23
void calculate_mh_sens(const initial_spectrum* is, double L, h_type H, double** dc2_out, double theta=0., int N=21) {
	
	// Values of d_cp
	double d_cp[N];
	

	double* dc2 = (double*)malloc(N * sizeof(double));

	// We need these for normalizations
	parameters nh0, ih0;
	nh0.populate(NH, 0.);
	ih0.populate(NH, 0.);

	ih0.flip_hierarchy();

	spectrum nh0s, ih0s;
	nh0s.p = &nh0;
	ih0s.p = &ih0;
	nh0s.oscillate(is, L);
	ih0s.oscillate(is, L);

	// Normalization for the NH spectra
	double nh0_int = nh0s.get_integral();
	nh0s.normalize(nh0_int, NH);

	double ih0_int = ih0s.get_integral();
	ih0s.normalize(ih0_int, IH);



	parameters nh[N], ih[N];

	spectrum osc_nhs[N];
	spectrum osc_ihs[N];


	// Initialize the hierarchies
	for (int i=0; i<N; ++i) {
		Printf("%d/%d", i+1, N);
		d_cp[i] = ((float)i / (float)(N-1) * 2. - 1.) * pi;

		nh[i].populate(H, d_cp[i]);
		ih[i].populate(H, d_cp[i]);

		// Change theta23
		nh[i].t23 = theta;
		ih[i].t23 = theta;
		nh[i].populate_common();
		ih[i].populate_common();

		// Flip MH of ih hierarchies
		ih[i].flip_hierarchy();

		osc_nhs[i].p = &nh[i];
		osc_ihs[i].p = &ih[i];
		osc_nhs[i].oscillate(is, L);
		osc_ihs[i].oscillate(is, L);


		if (H == NH) {
			osc_nhs[i].normalize(nh0_int, NH);
			osc_ihs[i].normalize(ih0_int, IH);
		}
		else {
			osc_nhs[i].normalize(ih0_int, IH);
			osc_ihs[i].normalize(nh0_int, NH);
		}
			
		//osc_nhs[i].reconstruct(300);
		//osc_ihs[i].reconstruct(300);

	}
	
	// TODO
	// if we reconstruct all spectra, we get pretty nice results, but the sensitivity is
	// very low (and the code takes long to run).
	// maybe try to play around with energy reconstruction 

	// Calculate sensitivity
	// Loop over delta_CP values
	for (int i=0; i<N; ++i) {
		
		// for each d_cp we want to calculate the minimum corresponding delta chi squared
		double min_dc2 = 1e28;

		// So we perform another loop over delta_CP and take the minimum chi squared we find
		for (int j=0; j<N; ++j) {

			// When we assume normal hierarchy, the NH spectrum is fixed in i and we
			// search the IH with j
			double dc2_n = mean_dc2(osc_ihs[j].e, 
									osc_nhs[i].e);

			min_dc2 = min(dc2_n, min_dc2);

		}

		dc2[i] = sqrt(min_dc2);
	}


	// Set the output values but check for NULL pointers
	// Caller function should remember to free these memory blocks
	if (dc2_out)
		*dc2_out = dc2;
	else
		free(dc2);

}

void do_dc2_theta(const initial_spectrum* is, double L) {
	const int N=21;

	double theta23 = 0.738;
	//double theta23_hi = 0.738 * (1 + 3 * .059);
	//double theta23_lo = 0.738 * (1 - 3 * .059);
	double theta23_hi = 0.798;
	double theta23_lo = 0.683;
	
	double *dc2_mean, *dc2_hi, *dc2_lo;

	calculate_mh_sens(is, L, NH, &dc2_mean, theta23, N);
	calculate_mh_sens(is, L, NH, &dc2_hi, theta23_hi, N);
	calculate_mh_sens(is, L, NH, &dc2_lo, theta23_lo, N);


	// Plotting
	double x[N];
	double d_cp[N];
	for (int i=0; i<N; ++i) {
		d_cp[i] = ((float)i / (float)(N-1) * 2. - 1.) * pi;
		x[i] = d_cp[i] / pi;
	}
	

	TCanvas* c4 = new TCanvas("c9", "", 600, 500);
	c4->SetTitle("MH sensitivity");
	//c4->SetFillColor(ci[CI_BACKGROUND]);
	
	TGraph* gdc2_mean = new TGraph(N, x, dc2_mean);
	double e_y[N];
	double e_e[N];

	for (int i=0; i<N; ++i) {
		e_y[i] = (dc2_hi[i] + dc2_lo[i]) / 2.;
		e_e[i] = (dc2_hi[i] - dc2_lo[i]) / 2.;
	}

	TGraphErrors* ge = new TGraphErrors(N, x, e_y, NULL, e_e);

	ge->SetTitle("True normal hierarchy");
	ge->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	ge->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	ge->GetXaxis()->SetLimits(-1., 1.);
	ge->SetMinimum(0);
	//gdc2_mean->SetMaximum(20);
	gdc2_mean->SetLineColor(4);
	gdc2_mean->SetMarkerColor(4);
	gdc2_mean->SetLineWidth(2);

	ge->SetFillStyle(1001);
	ge->SetFillColor(ci[CI_NH]);
	ge->Draw("a3");
	gdc2_mean->Draw("same");


	free(dc2_mean);
	free(dc2_hi);
	free(dc2_lo);
}


/*void read_spectrum(initial_spectrum* s) {
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
}*/

/*void plot_initial(float* mu, float* antimu, float* e, float* antie) {
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

}*/

