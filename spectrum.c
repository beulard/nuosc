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

// Store the integrals of the fluxes in the 'integrals' array (size N_RATES)
void get_integrals(spectrum* s, double* integrals) {
	for (int i=0; i<N_RATES; ++i) {
		for (int j=0; j<Nbins; ++j) {
			integrals[i] += s->events[i][j];
		}
	}
}

// Normalize a spectrum to the CDR rates given the integrals of the d_cp=0 spectra
void normalize(spectrum* s, const double* integrals) {
	for (int i=0; i<N_RATES; ++i) {
		for (int j=0; j<Nbins; ++j) {
			s->events[i][j] *= norms[s->h->type][i] / integrals[i];
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
			
			// If we want to implement matter effects, we need to recalculate
			// the MNS matrix for every energy since the C factor depends on E
			
			// TODO for antineutrinos, this is negative
			double A = 2.52e-28 * Ne * E / (os->h->dm2_31);
			double C = sqrt(pow(cos(2. * os->h->t13) - A, 2) + sinsq(2. * os->h->t12));

			// Now we change the effective mass splitting and mixing angles and 
			// recalculate everything. We don't want to modify the original hierarchy though.
			hierarchy h;
			memcpy(&h, os->h, sizeof(hierarchy));
			//populate_common(&h);
			//h.dm2_31 *= C;
			//h.t12 = 0.5 * asin(sin(2. * h.t12) / C);
			//populate_common(&h);
			//Printf("%f", h.t12 / TMath::Pi() * 180.);
			
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
		//	To get a more realistic picture, we normalize these to predicted event rates.


	}


	// We want to store a d_cp=0 oscillated spectrum so we can normalize others to it
	static int call = 0;
	static hierarchy h_0[2];
	static spectrum s_0[2];
	static double integrals_0[2][N_RATES];
	if (call == 0) {
		call = 1;
		for (int k=0; k<2; ++k) {
			populate(&h_0[k], (h_type)k, 0.);

			s_0[k].h = &h_0[k];
			oscillate(is, &s_0[k], false);

			get_integrals(&s_0[k], integrals_0[k]);
		}
	}

	if (normal) {
		normalize(os, integrals_0[os->h->type]);
	}
}
