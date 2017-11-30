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

/*
double P(flavor a, flavor b, double E, hierarchy* h, bool anti) {
	double x = L/E;
	double p = (a == b ? 1.0 : 0.0);

	for (int i=0; i<3; ++i) {
		for (int j=0; j<3; ++j) {
			if (i > j) {
				complex<double> K = h->MNS[a*3 + i] * conj(h->MNS[b*3 + i]) 
					  * conj(h->MNS[a*3 + j]) * h->MNS[b*3 + j];
					
				p -= 4.0 * real(K) * sinsq(1.2668 * h->dm2_mat[i*3 + j] * x);
				p += 4.0 * imag(K) * TMath::Sin(1.2668 * h->dm2_mat[i*3 + j] * x)
					     * TMath::Cos(1.2668 * h->dm2_mat[i*3 + j] * x);
			}
		}
	}
	return p;
}*/

// Electron density in mantle in units cm^-3 / (Avogadro's number)
const double Ne = 2.2;
// Oscillation probability for nu_mu -> nu_e or antinu_mu -> antinu_e
// including matter effects (Boris Kayser Neutrino Oscillation phenomenology p. 13)
double P_me(flavor fa, flavor fb, double E, hierarchy* h, bool anti) {

	double x = 1.53e-4 * Ne * E / h->dm2_31;
	double D = 1.269 * h->dm2_31 * L / E;
	double a = h->dm2_21 / h->dm2_31;
	a = 1. / 30.;
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

// Avogadro's number
//const double N_A = 6.02e23;
// Electron density in mantle
//const double Ne = 2.2 * N_A;

double P_me2(double E, hierarchy* h, bool anti) {
	double x = 2.52e-28 * Ne * E / (h->dm2_31);
	double D = 1.269 * h->dm2_31 * L / E;
	double A = h->dm2_21 / h->dm2_31;

	double T1 = sinsq(h->t23) * sinsq((1. - x) * D);
	double T2 = sin(h->d_cp) * sin(2. * h->t12) * sin(2. * h->t23)
				* sin(D) * sin(x * D) / x * sin((1. - x) * D) / (1. - x);
	double T3 = cos(h->d_cp) * sin(2 * h->t12) * sin(2. * h->t23)
				* cos(D) * sin(x * D)/x * sin((1. - x) * D)/(1. - x);
	double T4 = cossq(h->t23) * sinsq(2. * h->t12) * sinsq(x * D)/pow(x, 2);

	double p = sinsq(2. * h->t13) * T1 - A * sin(2. * h->t13) * T2
				+ A * sin(2. * h->t13) * T3 + pow(A, 2) * T4;

	return p;
}


// Oscillate neutrinos for each flavor and for each energy and get the FD flux.
void oscillate(const initial_spectrum* is, spectrum* os) {	
	
	initial_spectrum s;
	memcpy(&s, is, sizeof(initial_spectrum));

	// Propagated fluxes	
	float fluxes[N_FLUXES][Nbins] = {0};
	const float* norms;
	if (os->h->type == NH) {
		norms = nh_norm;
	} else if (os->h->type == IH) {
		norms = nh_norm;
	}

	// TODO check if this takes the correct values from the initial spectrum
	// when doing the oscillation, i.e. if it does start at 0.6 GeV and not at 0.

	//Printf("\n%f", os->h->t12 / TMath::Pi() * 180.);
	for (int i=0; i<Nbins; ++i) {
		// Instead of using the probability value at the discrete point given by i * 0.2,
		// we average the probabilities over the range (i*0.2 + eps, (i+1) * 0.2)
		//int samples = (int)(50. / (i * 0.2 + .3));
		int samples = 40;
		//Printf("\t%d", samples);
		double avg_P[N_FLUXES] = {0};
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

			avg_P[MU_SURVIVAL] += P(f_m, f_m, E, &h, false);
			avg_P[MU_E] += P_me(f_m, f_e, E, &h, false);
			//avg_P[MU_SURVIVAL] += 1. - P_me(E, &h, false);
			//avg_P[MU_E] += P_me(E, &h, false);
			//if(i==0) {
			//	Printf("%f", P_me2(E, &h, false));
			//}
			avg_P[MU_TAU] += P(f_m, f_t, E, &h, false);
			avg_P[E_SURVIVAL] += P(f_e, f_e, E, &h, false);
			avg_P[E_MU] += P(f_e, f_m, E, &h, false);
			avg_P[E_TAU] += P(f_e, f_t, E, &h, false);
			avg_P[ANTIMU_SURVIVAL] += P(f_m, f_m, E, &h, true);
			avg_P[ANTIMU_ANTIE] += P(f_m, f_e, E, &h, true);
			//avg_P[ANTIMU_SURVIVAL] += 1. - P_me(E, &h, true);
			//avg_P[ANTIMU_ANTIE] += P_me(E, &h, true);
			avg_P[ANTIMU_ANTITAU] += P(f_m, f_t, E, &h, true);
			avg_P[ANTIE_SURVIVAL] += P(f_e, f_e, E, &h, true);
			avg_P[ANTIE_ANTIMU] += P(f_e, f_m, E, &h, true);
			avg_P[ANTIE_ANTITAU] += P(f_e, f_t, E, &h, true);
		}
		for (int j=0; j<N_FLUXES; ++j) {
			avg_P[j] /= (double)samples;
		}

			
		// eg a muon neutrino will oscillate into electron neutrinos AND tau neutrinos
		// and leave a certain amount of muon neutrinos behind as well.
		//
		// Electron neutrinos will do the same thing so we must add the number of 
		// survival muon neutrinos to the number of oscillated electron neutrinos
		// to get the final number of muon neutrinos

		fluxes[MU_SURVIVAL][i] = s.mu[i] * avg_P[MU_SURVIVAL];
		fluxes[MU_E][i] = s.mu[i] * avg_P[MU_E];
		fluxes[MU_TAU][i] = s.mu[i] * avg_P[MU_TAU];
		fluxes[E_SURVIVAL][i] = s.e[i] * avg_P[E_SURVIVAL];
		fluxes[E_MU][i] = s.e[i] * avg_P[E_MU];
		fluxes[E_TAU][i] = s.e[i] * avg_P[E_TAU];
		fluxes[ANTIMU_SURVIVAL][i] = s.antimu[i] * avg_P[ANTIMU_SURVIVAL];
		fluxes[ANTIMU_ANTIE][i] = s.antimu[i] * avg_P[ANTIMU_ANTIE];
		fluxes[ANTIMU_ANTITAU][i] = s.antimu[i] * avg_P[ANTIMU_ANTITAU];
		fluxes[ANTIE_SURVIVAL][i] = s.antie[i] * avg_P[ANTIE_SURVIVAL];
		fluxes[ANTIE_ANTIMU][i] = s.antie[i] * avg_P[ANTIE_ANTIMU];
		fluxes[ANTIE_ANTITAU][i] = s.antie[i] * avg_P[ANTIE_ANTITAU];

	}

	// Now we normalize the fluxes to their respective FD predicted values
	float integrals[N_FLUXES] = {0};
	for (int i=0; i<N_FLUXES; ++i) {
		// For the nu e antinu e background, aka the survival of nu e's and antinu e's,
		// we only have one normalization, thus they need to be normalized together
		if (i == E_SURVIVAL) {
			// normalize e survival and anti e survival
			// add the fluxes and normalize the sum
			float e_antie_survival[Nbins] = {0};
			// array of fractions of the fluxes so we can retrieve them at the end
			float f[Nbins];
			for (int j=0; j<Nbins; ++j) {
				e_antie_survival[j] = fluxes[E_SURVIVAL][j] + fluxes[ANTIE_SURVIVAL][j];
				f[j] = fluxes[E_SURVIVAL][j] / fluxes[ANTIE_SURVIVAL][j];
			}
			float integral = 0.;
			for (int j=0; j<Nbins; ++j) {
				integral += e_antie_survival[j];
			}
			for (int j=0; j<Nbins; ++j) {
				e_antie_survival[j] /= integral;
			}
			// Now we have the normalized sum of fluxes but we can retrieve the individual
			// fluxes by using the following
			for (int j=0; j<Nbins; ++j) {
				fluxes[E_SURVIVAL][j] = e_antie_survival[j] * f[j] / (f[j] + 1.);
				fluxes[ANTIE_SURVIVAL][j] = e_antie_survival[j] / f[j] / (1./f[j] + 1.);

				// And normalize to data
				fluxes[E_SURVIVAL][j] *= norms[E_ANTIE_BACKGROUND];
				fluxes[ANTIE_SURVIVAL][j] *= norms[E_ANTIE_BACKGROUND];
			}
		} else if (i == ANTIE_SURVIVAL) {;} else {	
			for (int j=0; j<Nbins; ++j) {
				integrals[i] += fluxes[i][j];
			}
			// normalize to 1
			for (int j=0; j<Nbins; ++j) {
				fluxes[i][j] /= integrals[i];
			}
		}
	}
	// Then we multiply the fluxes by their respective normalizations
	for (int i=0; i<Nbins; ++i) {
		fluxes[MU_E][i] *= norms[E_SIGNAL];
		fluxes[ANTIMU_ANTIE][i] *= norms[ANTIE_SIGNAL];
		fluxes[MU_SURVIVAL][i] *= norms[MU_SIGNAL];
		fluxes[ANTIMU_SURVIVAL][i] *= norms[ANTIMU_BACKGROUND];

		// Omit fluxes that weren't normalized because they are negligible
		os->mu[i] = /*fluxes[E_MU][i]*/ + fluxes[MU_SURVIVAL][i];
		os->e[i] = fluxes[MU_E][i]/* + fluxes[E_SURVIVAL][i]*/;
		os->antimu[i] = /*fluxes[ANTIE_ANTIMU][i] +*/ fluxes[ANTIMU_SURVIVAL][i];
		os->antie[i] = fluxes[ANTIMU_ANTIE][i]/* + fluxes[ANTIE_SURVIVAL][i]*/;
	}
}
