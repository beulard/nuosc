#include "spectrum.h"
#include "helper.h"

// plotting and helper functions
void plot_initial(float* mu, float* antimu, float* e, float* antie);
void plot_normalized(float* mu, float* antimu, float* e, float* antie);
void plot_particles(float* mu, float* e);
void plot_antiparticles(float* antimu, float* antie);
void read_spectrum(initial_spectrum* s);
void plot_dc2(double* d_cp, double* dc2_mh_n, double* sd_mh_n, double* dc2_mh_i, double* sd_mh_i, double* dc2_cp, double* sd_cp_n, double* dc2_cp_ih, double* sd_cp_i, int N);
void do_dc2(const initial_spectrum* s);
void do_dc2_theta(const initial_spectrum* s);


// Mean delta chi-squared between a test spectrum and a "true" spectrum,
// as given in 1210.3651 p.9
double mean_dc2(const double* test, const double* tru, int N = 37) {
	double r = 0.;
	for (int i=0; i<N; ++i) {
		r += pow(test[i] - tru[i], 2) / tru[i];
	}
	return r;
}

// Statistical mean of a set of N values
double mean(double* vals, int N) {
	double r = 0.;
	for (int i=0; i<N; ++i) {
		r += vals[i] / N;
	}
	return r;
}

// Statistical standard deviation from the mean
double std_dev(double* vals, double mean, int N) {
	double r = 0.;
	for (int i=0; i<N; ++i) {
		r += pow(vals[i] - mean, 2) / (N-1);
	}

	return sqrt(r);
}

// theta=0: plot sensitivity with error band
// theta=1: plot sensitivity for best fit and 1sigma theta23 as band
void nuspectrum(int theta=0) {
	gStyle->SetOptStat(0);

	// initial spectrum in the LBNF's neutrino mode
	initial_spectrum nu;
	// and in the antineutrino mode
	initial_spectrum antinu;

	// read data from root file
	read_spectrum(&nu);
	
	// switch particle/antiparticle spectrum around for the antinu mode
	for (int i=0; i<50; ++i) {
		antinu.mu[i] = nu.antimu[i];
		antinu.antimu[i] = nu.mu[i];
		antinu.e[i] = nu.antie[i];
		antinu.antie[i] = nu.e[i];
	}
	
	// plot the initial spectrum (un-normalized)
	//plot_initial(nu.mu, nu.antimu, nu.e, nu.antie);
	//plot_initial(antinu.mu, antinu.antimu, antinu.e, antinu.antie);
	
	// Do the best-fit oscillation
	//oscillate(&nu, &best_fit_nu);
	

	// Main functions
	/*if(anti)
		do_dc2(&antinu);
	else
		do_dc2(&nu);
		*/
	if(theta == 0)
		do_dc2(&nu);
	else if(theta == 1)
		do_dc2_theta(&nu);
}

double get_integral(const spectrum* s) {
	double r=0.;
	for (int i=0; i<Nbins; ++i) {
		r += s->events[E_SIGNAL][i];
	}
	return r;
}

// Calculates MH mean delta chi squared and standard deviation given a hierarchy and theta_23
void calculate_mh_sens(const initial_spectrum* s, h_type H, double** dc2_out, double** sd_out, double theta=0., int N=61, int N_rec=100) {
	
	// Values of d_cp
	double d_cp[N];
	

	double dc2[N][N_rec];
	memset(&dc2, 0, N * N_rec * sizeof(double));

	// We need these for normalizations
	hierarchy nh0, ih0;
	populate(&nh0, NH, 0.);
	populate(&ih0, NH, 0.);
	nh0.t23 = theta;
	ih0.t23 = theta;
	populate_common(&nh0);
	populate_common(&ih0);

	flip_hierarchy(&ih0);

	spectrum nh0s, ih0s;
	nh0s.h = &nh0;
	ih0s.h = &ih0;
	oscillate(s, &nh0s);
	oscillate(s, &ih0s);

	hierarchy nh[N], ih[N];

	// Oscillated, but not reconstructed spectra
	spectrum osc_nhs[N];
	spectrum osc_ihs[N];
	memset(&osc_nhs, 0, N * sizeof(spectrum));
	memset(&osc_ihs, 0, N * sizeof(spectrum));

	// Normalization for the NH spectra
	double nh0_int = get_integral(&nh0s);
	normalize(&nh0s, nh0_int, NH);

	double ih0_int = get_integral(&ih0s);
	normalize(&ih0s, ih0_int, IH);


	// Reconstructed spectra, build from non-reconstructed ones
	spectrum *nhs, *ihs;
	nhs = (spectrum*)malloc(N * N_rec * sizeof(spectrum));
	ihs = (spectrum*)malloc(N * N_rec * sizeof(spectrum));

	// Initialize the hierarchies
	for (int i=0; i<N; ++i) {
		Printf("%d/%d", i+1, N);
		d_cp[i] = ((float)i / (float)(N-1) * 2. - 1.) * pi;

		populate(&nh[i], H, d_cp[i]);
		populate(&ih[i], H, d_cp[i]);

		// Change theta23
		nh[i].t23 = theta;
		ih[i].t23 = theta;
		populate_common(&nh[i]);
		populate_common(&ih[i]);

		// Flip MH of ih hierarchies
		flip_hierarchy(&ih[i]);

		osc_nhs[i].h = &nh[i];
		osc_ihs[i].h = &ih[i];
		oscillate(s, &osc_nhs[i]);
		oscillate(s, &osc_ihs[i]);

		normalize(&osc_nhs[i], nh0_int, NH);
		normalize(&osc_ihs[i], ih0_int, IH);

		// Reconstruct the spectra
		// It is important that we get ALL the spectra before we go
		// into the next loop because it will go through them all to calculate
		// sensitivity.
		for (int n=0; n<N_rec; ++n) {
			nhs[i + n*N].h = &nh[i];
			ihs[i + n*N].h = &ih[i];
			reconstruct(&osc_nhs[i], &nhs[i + n*N]);
			reconstruct(&osc_ihs[i], &ihs[i + n*N]);
		}
	}
	
	// TODO
	// if we reconstruct all spectra, we get pretty nice results, but the sensitivity is
	// very low (and the code takes long to run).
	// maybe try to play around with energy reconstruction 

	// Calculate sensitivity
	for (int n=0; n<N_rec; ++n) {
		// Loop over delta_CP values
		for (int i=0; i<N; ++i) {
			
			// for each d_cp we want to calculate the minimum corresponding delta chi squared
			double min_dc2 = 1e28;

			// So we perform another loop over delta_CP and take the minimum chi squared we find
			for (int j=0; j<N; ++j) {

				// When we assume normal hierarchy, the NH spectrum is fixed in i and we
				// search the IH with j
				double dc2_n = mean_dc2(ihs[n*N + j].events[E_SIGNAL], 
										nhs[n*N + i].events[E_SIGNAL]);

				min_dc2 = min(dc2_n, min_dc2);

			}

			dc2[i][n] = min_dc2;
		}
	}


	// Mean sensitivity
	double* mean_dc2_mh = (double*)malloc(N * sizeof(double));

	// Standard deviations from the mean
	double* sd_mh = (double*)malloc(N * sizeof(double));

	// After calculating all sensitivities, we want to create an average and a
	// standard deviation for each sample d_cp point
	for (int i=0; i<N; ++i) {
		mean_dc2_mh[i] = sqrt(mean(dc2[i], N_rec));

		sd_mh[i] = sqrt(std_dev(dc2[i], mean_dc2_mh[i], N_rec));
	}

	// Set the output values but check for NULL pointers
	// Caller function should remember to free these memory blocks
	if (dc2_out)
		*dc2_out = mean_dc2_mh;
	if (sd_out)
		*sd_out = sd_mh;

	free(nhs);
	free(ihs);
}

void do_dc2_theta(const initial_spectrum* s) {
	const int N=61;
	const int N_rec=300;

	double theta23 = 0.738;
	//double theta23_hi = 0.738 * (1 + 3 * .059);
	//double theta23_lo = 0.738 * (1 - 3 * .059);
	double theta23_hi = 0.927;
	double theta23_lo = 0.66;
	
	double *dc2_mean, *dc2_hi, *dc2_lo;

	calculate_mh_sens(s, NH, &dc2_mean, NULL, theta23, N, N_rec);
	calculate_mh_sens(s, NH, &dc2_hi, NULL, theta23_hi, N, N_rec);
	calculate_mh_sens(s, NH, &dc2_lo, NULL, theta23_lo, N, N_rec);


	// Plotting
	double x[N];
	double d_cp[N];
	for (int i=0; i<N; ++i) {
		d_cp[i] = ((float)i / (float)(N-1) * 2. - 1.) * pi;
		x[i] = d_cp[i] / pi;
	}
	

	TCanvas* c4 = new TCanvas("c9", "", 600, 500);
	c4->SetFillColor(ci[CI_BACKGROUND]);
	
	TGraph* gdc2_mean = new TGraph(N, x, dc2_mean);
	TGraph* gdc2_hi = new TGraph(N, x, dc2_hi);
	TGraph* gdc2_lo = new TGraph(N, x, dc2_lo);
	gdc2_mean->SetFillColor(ci[CI_NH]);
	gdc2_mean->Draw();
	gdc2_mean->SetTitle("True normal hierarchy");
	gdc2_mean->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mean->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_mean->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mean->SetMinimum(0);
	//gdc2_mean->SetMaximum(20);
	gdc2_mean->SetLineColor(4);
	gdc2_mean->SetMarkerColor(4);
	gdc2_mean->SetLineWidth(2);
	gdc2_hi->SetLineWidth(2);
	gdc2_hi->SetLineColor(3);
	gdc2_lo->SetLineWidth(2);
	//gdc2_lo->SetLineColor(5);
	gdc2_hi->Draw("same");
	gdc2_lo->Draw("same");

	c4->SetTitle("MH sensitivity");
}


void do_dc2(const initial_spectrum* s) {
	// Number of times we will evaluate the sensitivity to compose the final plot
	const int N_rec = 300;
	// Number of d_cp samples
	const int N = 61;
	// Values of d_cp
	double d_cp[N];

	// Delta chi squared for the mass hierarchy assuming true normal hierarchy
	double dc2_mh_n[N][N_rec] = {0};
	double dc2_mh_i[N][N_rec] = {0};
	// Delta chi squared for delta_CP
	double dc2_cp_n[N][N_rec] = {0};
	double dc2_cp_i[N][N_rec] = {0};


	// Calculate mean delta chi squared at each d_cp
	hierarchy nh0, nhpi;
	// Hierarchy that has NH parameters (theta23) but inverted mass hierarchy
	hierarchy nhi0;
	populate(&nh0, NH, 0.);
	populate(&nhpi, NH, pi);
	populate(&nhi0, NH, 0.);
	flip_hierarchy(&nhi0);

	spectrum nh0s, nhpis;
	spectrum nhi0s;
	nh0s.h = &nh0;
	nhpis.h = &nhpi;
	nhi0s.h = &nhi0;
	oscillate(s, &nh0s);
	oscillate(s, &nhpis);
	oscillate(s, &nhi0s);

	// Integral for d_cp = 0 that we use to normalize spectra under NH assumption
	double nh0_int = get_integral(&nh0s);
	normalize(&nh0s, nh0_int, NH);
	normalize(&nhpis, nh0_int, NH);

	double nhi0_int = get_integral(&nhi0s);
	// We normalize nhi0s to the IH event rates because it is oscillated under IH
	normalize(&nhi0s, nhi0_int, IH);

	reconstruct(&nh0s, &nh0s, 100);
	reconstruct(&nhpis, &nhpis, 100);
	reconstruct(&nhi0s, &nhi0s, 100);


	hierarchy ih0, ihpi;
	// Hierarchy that has IH parameters (theta23) but NH mass hierarchy
	hierarchy ihn0;
	populate(&ih0, IH, 0.);
	populate(&ihpi, IH, pi);
	populate(&ihn0, IH, 0.);
	flip_hierarchy(&ihn0);

	spectrum ih0s, ihpis;
	spectrum ihn0s;
	ih0s.h = &ih0;
	ihpis.h = &ihpi;
	ihn0s.h = &ihn0;
	oscillate(s, &ih0s);
	oscillate(s, &ihpis);
	oscillate(s, &ihn0s);

	double ih0_int = get_integral(&ih0s);
	normalize(&ih0s, ih0_int, IH);
	normalize(&ihpis, ih0_int, IH);

	double ihn0_int = get_integral(&ihn0s);
	normalize(&ihn0s, ihn0_int, NH);

	reconstruct(&ih0s, &ih0s, 100);
	reconstruct(&ihpis, &ihpis, 100);
	reconstruct(&ihn0s, &ihn0s, 100);

	Printf("%e %e %e %e", nh0_int, nhi0_int, ih0_int, ihn0_int);

	// N * 2 hierarchies
	// The nh array contains hierarchies that will be used in calculating
	// sensitivity under the assumption that the mass hierarchy is normally ordered.
	// Technically, the second element (nh[i][1]) has inverted ordering.
	hierarchy nh[N][2], ih[N][2];

	// Oscillated, but not reconstructed spectra
	spectrum osc_nhs[N][2] = {0};
	spectrum osc_ihs[N][2] = {0};

	// Reconstructed spectra, build from non-reconstructed ones
	spectrum *nhs, *ihs;
	nhs = (spectrum*)malloc(N * N_rec * sizeof(spectrum) * 2);
	ihs = (spectrum*)malloc(N * N_rec * sizeof(spectrum) * 2);
	

	// Initialize the hierarchies
	for (int i=0; i<N; ++i) {
		Printf("%d/%d", i+1, N);
		d_cp[i] = ((float)i / (float)(N-1) * 2. - 1.) * pi;

		for (int j=0; j<2; ++j) {
			populate(&nh[i][j], NH, d_cp[i]);
			populate(&ih[i][j], IH, d_cp[i]);
		}
		flip_hierarchy(&nh[i][1]);
		flip_hierarchy(&ih[i][1]);


		for (int j=0; j<2; ++j) {
			osc_nhs[i][j].h = &nh[i][j];
			osc_ihs[i][j].h = &ih[i][j];
			oscillate(s, &osc_nhs[i][j]);
			oscillate(s, &osc_ihs[i][j]);

			// Normalize NH spectra to NH0 and IH spectra to IH0
			if (j==0) {
				normalize(&osc_nhs[i][j], nh0_int, NH);
				normalize(&osc_ihs[i][j], ih0_int, IH);
			}
			else {
				normalize(&osc_nhs[i][j], nhi0_int, IH);
				normalize(&osc_ihs[i][j], ihn0_int, NH);
			}

			// Oscillate/reconstruct the spectra
			// It is important that we get ALL the spectra before we go
			// into the next loop because it will go through them all to calculate
			// sensitivity.
			for (int n=0; n<N_rec; ++n) {
				nhs[n*N*2 + i*2 + j].h = &nh[i][j];
				ihs[n*N*2 + i*2 + j].h = &ih[i][j];
				reconstruct(&osc_nhs[i][j], &nhs[n*N*2 + i*2 + j], 3);
				reconstruct(&osc_ihs[i][j], &ihs[n*N*2 + i*2 + j], 3);
			}
		}
	}
	Printf("%f %f", get_integral(&osc_ihs[N/2][0]), get_integral(&osc_ihs[N/2][1]));


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
										osc_nhs[j][1].events[E_SIGNAL],
										nhs[n*N*2 + i*2 + 0].events[E_SIGNAL]);
				// Vice versa
				double dc2_i = mean_dc2(//ihs[n*N*2 + j*2 + 1].events[E_SIGNAL], 
										osc_ihs[j][1].events[E_SIGNAL],
										ihs[n*N*2 + i*2 + 0].events[E_SIGNAL]);

				

				min_dc2_n = min(dc2_n, min_dc2_n);
				min_dc2_i = min(dc2_i, min_dc2_i);

			}

			dc2_mh_n[i][n] = min_dc2_n;
			dc2_mh_i[i][n] = min_dc2_i;

			// For the CP sensitivity we want the j=0 elements (the ones with
			// the correct hierarchy)
			dc2_cp_n[i][n] = min(mean_dc2(nhs[n*N*2 + i*2 + 0].events[E_SIGNAL],
					   				      nh0s.events[E_SIGNAL]),
				   	 	         mean_dc2(nhs[n*N*2 + i*2 + 0].events[E_SIGNAL], 
								          nhpis.events[E_SIGNAL]));

			dc2_cp_i[i][n] = min(mean_dc2(ihs[n*N*2 + i*2 + 0].events[E_SIGNAL], 
									      ih0s.events[E_SIGNAL]),
							     mean_dc2(ihs[n*N*2 + i*2 + 0].events[E_SIGNAL], 
								   	      ihpis.events[E_SIGNAL]));
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
		mean_dc2_mh_n[i] = mean(dc2_mh_n[i], N_rec);
		mean_dc2_mh_i[i] = mean(dc2_mh_i[i], N_rec);
		mean_dc2_cp_n[i] = mean(dc2_cp_n[i], N_rec);
		mean_dc2_cp_i[i] = mean(dc2_cp_i[i], N_rec);

		sd_mh_n[i] = sqrt(std_dev(dc2_mh_n[i], mean_dc2_mh_n[i], N_rec));
		sd_mh_i[i] = sqrt(std_dev(dc2_mh_i[i], mean_dc2_mh_i[i], N_rec));
		sd_cp_n[i] = sqrt(std_dev(dc2_cp_n[i], mean_dc2_cp_n[i], N_rec));
		sd_cp_i[i] = sqrt(std_dev(dc2_cp_i[i], mean_dc2_cp_i[i], N_rec));
	}	

	// To plot a test spectrum
	
	TCanvas* ccc = new TCanvas();
	//ccc->SetLogy();
	TH1* hnh = new TH1F("hnh", "", Nbins, 0.6, 8.);
	TH1* hih = new TH1F("hih", "", Nbins, 0.6, 8.);
	int index = (int)((0.4 / 2. + .5) * N);

	//Printf("%f", d_cp[index] / pi);
	for (int i=0; i<Nbins; ++i) {
		hnh->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, nhs[index*2].events[E_SIGNAL][i]);
		//hnh->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, ih0s.events[E_SIGNAL][i]);
		hih->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, nh0s.events[E_SIGNAL][i]);
		//hih->Fill(0.6 + (8. - 0.6) * (float)i / Nbins + 1e-3, ih0s.events[E_SIGNAL][i]);
	}

	hnh->Draw("hist");
	hnh->SetMaximum(110);
	hih->Draw("hist same");
	hnh->SetLineWidth(2);
	hih->SetLineColor(2);
	

	
	plot_dc2(d_cp, mean_dc2_mh_n, sd_mh_n, mean_dc2_mh_i, sd_mh_i, mean_dc2_cp_n, sd_cp_n,  mean_dc2_cp_i, sd_cp_i, N);
	

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
		ymnh[i] = sqrt(dc2_mh_n[i]);
		ymih[i] = sqrt(dc2_mh_i[i]);
		ydnh[i] = sqrt(dc2_cp[i]);
		ydih[i] = sqrt(dc2_cp_ih[i]);
	}


	TCanvas* c4 = new TCanvas("c4", "", 1000, 400);
	c4->SetFillColor(ci[CI_BACKGROUND]);
	c4->Divide(2, 1);
	
	TPad* p = (TPad*)c4->GetPad(1);
	p->cd();
	TGraph* gdc2_mnh = new TGraph(N, x, ymnh);
	TGraph* gdc2_mnh_e = new TGraphErrors(N, x, ymnh, NULL, sd_mh_n);
	gdc2_mnh_e->SetFillColor(ci[CI_NH]);
	gdc2_mnh_e->Draw("a4");
	gdc2_mnh->Draw("same");
	gdc2_mnh_e->SetTitle("True normal hierarchy");
	gdc2_mnh_e->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mnh_e->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_mnh_e->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mnh_e->SetMinimum(0);
	gdc2_mnh->SetLineColor(4);
	gdc2_mnh->SetMarkerColor(4);
	gdc2_mnh->SetLineWidth(3);

	TGraph* gdc2_mih = new TGraph(N, x, ymih);
	TGraph* gdc2_mih_e = new TGraphErrors(N, x, ymih, NULL, sd_mh_i);
	gdc2_mih_e->SetTitle("True inverted hierarchy");
	c4->SetTitle("MH sensitivity");
	gdc2_mih_e->GetYaxis()->SetTitle("#sqrt{#bar{#Delta #chi^{2}}}");
	gdc2_mih_e->GetXaxis()->SetTitle("#delta_{CP} / #pi");
	gdc2_mih->SetLineColor(kRed+1);
	gdc2_mih->SetMarkerColor(2);
	gdc2_mih_e->SetFillColor(ci[CI_IH]);
	//gdc2_mnh->SetMaximum(25);
	//gdc2_mih->SetMaximum(25);
	gdc2_mih->SetLineWidth(3);
 	p = (TPad*)c4->GetPad(2);	
	p->cd();
	gdc2_mih_e->GetXaxis()->SetLimits(-1., 1.);
	gdc2_mih_e->SetMinimum(0);
	//gdc2_mih->SetMaximum(25);
	gdc2_mih_e->Draw("a4");
	gdc2_mih->Draw("same");
	

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


void read_spectrum(initial_spectrum* s) {
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

void plot_initial(float* mu, float* antimu, float* e, float* antie) {
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

}

