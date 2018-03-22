#include "spectrum.h"

double mean_dc2(const double* test, const double* tru, int N=37) {
	double r=0;

	for (int i=0; i<N; ++i) {
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
	const int N = 21;
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

	// Oscillate these spectra with specific parameters
	nh0s.oscillate(&is, L);
	nhpis.oscillate(&is, L);

	// Predicted event rate, from HyperK design report p.181
	const double e_nh = 2300;

	// Integral for d_cp = 0 that we use to normalize spectra under NH assumption
	double nh0_int = nh0s.get_integral();
	nh0s.normalize(nh0_int, e_nh);
	nhpis.normalize(nh0_int, e_nh);

	// IH Parameters combinations for normalization and d_cp sensitivity
	parameters ih0, ihpi;
	ih0.populate(IH, 0.);
	ihpi.populate(IH, pi);

	spectrum ih0s, ihpis;
	ih0s.p = &ih0;
	ihpis.p = &ihpi;
	ih0s.oscillate(&is, L);
	ihpis.oscillate(&is, L);

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
		double delta_E = (1.3 - 0.1) / samples;
		x[i] = i * delta_E + 0.1;

		int_nh += P_me(f_m, f_e, x[i], L, &nh0, false);
		int_ih += P_me(f_m, f_e, x[i], L, &ih0, false);
	}
	const double e_ih = e_nh * int_ih / int_nh;
	//Printf("Norms: %f %f", e_nh, e_ih);
	

	ih0s.normalize(ih0_int, e_ih);
	ihpis.normalize(ih0_int, e_ih);


	parameters nh[N], ih[N];

	spectrum osc_nhs[N] = {0};
	spectrum osc_ihs[N] = {0};

	for (int i=0; i<N; ++i) {
		Printf("%d/%d", i+1, N);
		d_cp[i] = ((double)i / (double)(N-1) * 2. - 1.) * pi;

		nh[i].populate(NH, d_cp[i]);
		ih[i].populate(IH, d_cp[i]);

		osc_nhs[i].p = &nh[i];
		osc_ihs[i].p = &ih[i];
		osc_nhs[i].oscillate(&is, L);
		osc_ihs[i].oscillate(&is, L);

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
			double dc2_n = mean_dc2(osc_ihs[j].events[E_SIGNAL],
									osc_nhs[i].events[E_SIGNAL]);

			// Vice versa
			double dc2_i = mean_dc2(osc_nhs[j].events[E_SIGNAL],
									osc_ihs[i].events[E_SIGNAL]);

			

			min_dc2_n = min(dc2_n, min_dc2_n);
			min_dc2_i = min(dc2_i, min_dc2_i);

		}

		dc2_mh_n[i] = sqrt(min_dc2_n);
		dc2_mh_i[i] = sqrt(min_dc2_i);

		// For the CP sensitivity we want the j=0 elements (the ones with
		// the correct hierarchy)
		dc2_cp_n[i] = sqrt(min(min(mean_dc2(osc_nhs[i].events[E_SIGNAL],
									  osc_nhs[N/2].events[E_SIGNAL]),
			   	 	         mean_dc2(osc_nhs[i].events[E_SIGNAL],
									  osc_nhs[0].events[E_SIGNAL])),
							 mean_dc2(osc_nhs[i].events[E_SIGNAL],
								 	  osc_nhs[N-1].events[E_SIGNAL])));

		dc2_cp_i[i] = sqrt(min(min(mean_dc2(//ihs[n*N*2 + i*2 + 0].events[E_SIGNAL], 
										   osc_ihs[i].events[E_SIGNAL],
								      //ih0s.events[E_SIGNAL]),
									  osc_ihs[N/2].events[E_SIGNAL]),
						     mean_dc2(//ihs[n*N*2 + i*2 + 0].events[E_SIGNAL], 
								      osc_ihs[i].events[E_SIGNAL],
							   	      //ihpis.events[E_SIGNAL])));
									  osc_ihs[0].events[E_SIGNAL])),
						    	mean_dc2(osc_ihs[i].events[E_SIGNAL],
										 osc_ihs[N-1].events[E_SIGNAL])));
	}

	
}
