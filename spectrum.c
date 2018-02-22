#include "spectrum.h"

// calculates probability of a neutrino flavor a to oscillate to flavor b,
// given an energy E, a hierarchy h and a constant baseline
double P(flavor a, flavor b, double E, hierarchy* h, bool anti) {
	double p = 0.;

	for(int i=0; i<3; ++i) {
		p += pow(abs(h->MNS[a*3 + i] * conj(h->MNS[b*3 + i])), 2);
	}

	if (anti) {
		for (int i=0; i<3; ++i) {
			for (int j=0; j<3; ++j) {
				if (j > i) {
					p += 2. * real( conj(h->MNS[a*3 + i]) * h->MNS[a*3 + j]
						  * h->MNS[b*3 + i] * conj(h->MNS[b*3 + j])
						  * exp(complex<double>(-2.i * 1.2668 * h->dm2_mat[i*3 + j] * L/E)));
				}
			}
		}
	}
	else {
		for(int i=0; i<3; ++i) {
			for(int j=0; j<3; ++j) {
				if(j > i) {
					p += 2. * real(h->MNS[a*3 + i] * conj(h->MNS[a*3 + j]) 
							* conj(h->MNS[b*3 + i]) * h->MNS[b*3 + j]
							* exp(complex<double>(-2.i * 1.2668 *  h->dm2_mat[i*3 + j] * L/E)));
				}
			}
		}
	}
	return p;
}


// Electron density in mantle in units cm^-3 / (Avogadro's number)
const double Ne = 2.2;
// Oscillation probability for nu_mu -> nu_e or antinu_mu -> antinu_e
// including matter effects (Boris Kayser Neutrino Oscillation phenomenology p. 13)
double P_me(flavor fa, flavor fb, double E, hierarchy* h, bool anti) {

	double x = 1.53e-4 * Ne * E / h->dm2_31;
	double D = 1.269 * h->dm2_31 * L / E;
	double a = h->dm2_21 / h->dm2_31;
	x = anti ? -x : x;
	double d_cp = anti ? -h->d_cp : h->d_cp;

	double T1 = sinsq(h->t23) * sinsq((1. - x) * D) / pow(1. - x, 2);
	double T2 = sin(d_cp) * sin(2 * h->t12) * sin(2 * h->t23) * sin(D)
				* sin(x * D) / x * sin((1 - x) * D) / (1 - x);
	double T3 = cos(d_cp) * sin(2 * h->t12) * sin(2 * h->t23) * cos(D)
				* sin(x * D) / x * sin((1 - x) * D) / (1 - x);
	double T4 = cossq(h->t23) * sinsq(2 * h->t12) * sinsq(x * D) / pow(x, 2);

	double p = sinsq(2 * h->t13) * T1 - a * sin(2 * h->t13) * T2 + a * sinsq(2 * h->t13) * T3
			+ pow(a, 2) * T4;

	return p;
}

// CDR transition probability
double P_me2(flavor fa, flavor fb, double E, hierarchy* h, bool anti) {
	double aL = 1.93e-4 * Ne * L;
	aL = anti ? -aL : aL;

	double d_cp = anti ? -h->d_cp : h->d_cp;
	double D31 = 1.269 * h->dm2_31 * L / E;
	double D21 = 1.269 * h->dm2_21 * L / E;

	double p = sinsq(h->t23) * sinsq(2*h->t13) * sinsq(D31 - aL) / pow(D31 - aL, 2) * pow(D31, 2)

		+ sin(2*h->t23) * sin(2 * h->t13) * sin(2*h->t12) * sin(D31 - aL) / (D31 - aL) * D31 * sin(aL) / aL * D21 * cos(D31 + d_cp)

		+ cossq(h->t23) * sinsq(2*h->t12) * sinsq(aL) / pow(aL, 2) * pow(D21, 2);
	return p;
}

// Store the integrals of the fluxes in the 'integrals' array (size N_RATES)
void get_integrals(spectrum* s, double* integrals) {
	for (int i=0; i<N_RATES; ++i) {
		for (int j=0; j<Nbins; ++j) {
			integrals[i] += s->events[i][j];
		}
	}
}

// Oscillate neutrinos for each flavor and for each energy and get the FD flux.
void oscillate(const initial_spectrum* is, spectrum* os, bool normal, bool antimode) {	
	

	initial_spectrum s;
	memcpy(&s, is, sizeof(initial_spectrum));


	for (int i=0; i<Nbins; ++i) {
		// Instead of using the probability value at the discrete point given by i * 0.2,
		// we average the probabilities over the range (i*0.2 + eps, (i+1) * 0.2)
		int samples = 50;

		double avg_P[N_RATES] = {0};
		for (int j=0; j<samples; ++j) {
			double E = 0.2 * (firstbin + i + (double)j / (double)samples);
			
			
			// Now we change the effective mass splitting and mixing angles and 
			// recalculate everything. We don't want to modify the original hierarchy though.
			hierarchy h;
			memcpy(&h, os->h, sizeof(hierarchy));
			
			avg_P[E_SIGNAL] += P_me(f_m, f_e, E, &h, false);
			avg_P[ANTIE_SIGNAL] += P_me(f_m, f_e, E, &h, true);

			avg_P[MU_SIGNAL] += P(f_m, f_m, E, &h, false);

		}
		for (int j=0; j<N_RATES; ++j) {
			avg_P[j] /= (double)samples;
		}

			

		//	Number of events at FD if we could detect every single neutrino.
		os->events[E_SIGNAL][i] = s.mu[i + firstbin] * avg_P[E_SIGNAL];
		os->events[ANTIE_SIGNAL][i] = s.antimu[i + firstbin] * avg_P[ANTIE_SIGNAL];
		os->events[MU_SIGNAL][i] = s.mu[i + firstbin] * avg_P[MU_SIGNAL];

	}


	// We want to store a d_cp=0 oscillated spectrum so we can normalize others to it
	/*static int call = 0;
	static hierarchy h_0[2];
	static spectrum s_0[2];
	static double integrals_0[2][N_RATES];
	if (call == 0) {
		call = 1;
		memset(integrals_0, 0, 2 * N_RATES * sizeof(double));
		for (int k=0; k<2; ++k) {
			populate(&h_0[k], (h_type)k, 0.);

			s_0[k].h = &h_0[k];
			oscillate(is, &s_0[k], false);


			get_integrals(&s_0[k], integrals_0[k]);
		}
	}*/

	// /!\ CHANGE: normalization should not be performed in this function,
	
	// If the normal flag is set, we normalize the spectrum to the d_cp=0 one
	//if (normal) {
	//	normalize(os, integrals_0[os->h->type]);
	//}
}

// Energy reconstruction function with smearing of the reconstructed energy to
// simulate experimental uncertainties
void reconstruct(const spectrum* s, spectrum* os, int smooth) {

	// Random number generator
	static int seed = 0;
	TRandom ra(seed);
	seed++;


	// TODO instead of smoothing a spectrum
	// repeat the sensitivity plot for every reconstructed spectrum
	// and take the mean as the center and evaluate the standard deviation
	// and plot that as a band around the mean
	//
	// TODO try to plot sensitivity for 'extreme' values for the mixing
	// angle (which? check out CDR) and plot a band around the sensitivity
	// at the best fit

	const int N_rec = smooth;

	// Reconstructed spectrum that will overwrite the original one
	spectrum s_f = {0};

	
	for (int n=0; n<N_rec; ++n) {
		// Temporary reconstructed spectrum
		spectrum s_rec = {0};

		// For each event, we displace the energy by a random number from a gaussian 
		// distribution
		for (int i=0; i<Nbins; ++i) {
			// True value of the energy (center of each bin)
			double E = 0.2 * (firstbin + i) + 0.1;
			//Printf("%f", E);
			for (int j=0; j<ceil(s->events[E_SIGNAL][i]); ++j) {
				// Can play around with the sigma on the gaussian for different results
				double E_rec = ra.Gaus(E, 0.35);
				//double E_rec = ra.Uniform(0.5 * E, E + 1.5); 
				// This actually gives nice results...
				//double E_rec = ra.Uniform(0.6, 8);
				//double E_rec = E;
				
				// After picking E_rec from a uniform distribution,
				// we take it slightly closer to the true value
				// as indicated by the DUNE simulations
				// This increases the overall sensitivity 
				// by a small amount the more we move E_rec but it also
				// makes the spectrum more similar to those of DUNE.
				//E_rec += (E - E_rec) / 1.5;
		
				// Find the bin corresponding to the reconstructed energy
				int bin = round((E_rec - 0.7) / 7.4 * Nbins + 1e-3);

				// If our reconstructed energy was out of bounds, 
				// we just assign the same bin.
				// This gives the true energy a slightly better chance at being found
				if (bin >= Nbins)
					bin = i;
				if (bin < 0)
					bin = i;
		
				// The value of each event is 1, except for the last one which might be
				// worth less because of our floating point precision
				double value = 1.;
				if (j == ceil(s->events[E_SIGNAL][i]) - 1){
					value = s->events[E_SIGNAL][i] - ceil(s->events[E_SIGNAL][i]) + 1.;
				}
				// We want to move this event to the bin we picked
				s_rec.events[E_SIGNAL][bin] += value;
			}
			// Avoid zero events for stability (crashes in delta chi^2 calculation)
			if (s_rec.events[E_SIGNAL][i] < .1)
				s_rec.events[E_SIGNAL][i] = 1.;
		}
		
		for (int i=0; i<Nbins; ++i) {
			s_f.events[E_SIGNAL][i] += s_rec.events[E_SIGNAL][i] / N_rec;
		}

	}
	// After obtaining final spectrum, copy it over to output
	memcpy(os->events, s_f.events, Nbins * N_RATES * sizeof(double));
}


void normalize(spectrum* s, double integral, h_type h) {
	for (int i=0; i<Nbins; ++i) {
		s->events[E_SIGNAL][i] *= norms[h][E_SIGNAL] / integral;
	}
}
