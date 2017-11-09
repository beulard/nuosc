#include "spectrum.h"

// calculates probability of a neutrino flavor a to oscillate to flavor b,
// given an energy E, a hierarchy h and a constant baseline
double P(flavor a, flavor b, float E, hierarchy* h, bool anti) {
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
					p += 2. * real(h->MNS[a*3 + i] * std::conj(h->MNS[a*3 + j]) 
							* conj(h->MNS[b*3 + i]) * h->MNS[b*3 + j]
							* exp(complex<double>(-2.i * 1.2668 *  h->dm2_mat[i*3 + j] * L/E)));
				}
			}
		}
	}
	return p;
}


enum spectrum_flux {
	MU_SURVIVAL,	/* mu_signal */
	MU_E,			/* e_signal */
	MU_TAU,			// no tau for now
	E_SURVIVAL,		/* part of e_antie_background */
	E_MU,			// e->mu probability is too small
	E_TAU,			// no tau for now
	ANTIMU_SURVIVAL,/* antimu_background */
	ANTIMU_ANTIE,	/* antie_signal */
	ANTIMU_ANTITAU,	// no tau for now
	ANTIE_SURVIVAL,	/* part of e_antie_background */
	ANTIE_ANTIMU,	// e->mu probability is too small
	ANTIE_ANTITAU,	// no tau for now
	N_FLUXES
};

// Propagate neutrinos for each flavor and for each energy and get the FD flux.
void propagate(const initial_spectrum* is, spectrum* os) {	
	hierarchy* h = &(os->h);
	
	initial_spectrum s;
	memcpy(&s, is, sizeof(initial_spectrum));

	// Propagated fluxes	

	float fluxes[N_FLUXES][50];

	for (int i=0; i<50; ++i) {
		float E = 0.2 * i;
		// avoid zero energy
		if (i == 0) E = 1e-12;

			
		// eg a muon neutrino will oscillate into electron neutrinos AND tau neutrinos
		// and leave a certain amount of muon neutrinos behind as well.
		//
		// Electron neutrinos will do the same thing so we must add the number of 
		// survival muon neutrinos to the number of oscillated electron neutrinos
		// to get the final number of muon neutrinos


		// Find proportion of survival/oscillated muon neutrinos
		// the proportion of survival MUON nus is
		fluxes[MU_SURVIVAL][i] = s.mu[i] * P(f_m, f_m, E, h, false);
		// proportion of oscillated mu -> electron neutrinos
		fluxes[MU_E][i] = s.mu[i] * P(f_m, f_e, E, h, false);
		// oscillated mu -> taus
		fluxes[MU_TAU][i] = s.mu[i] * P(f_m, f_t, E, h, false);

		// same for electron
		fluxes[E_SURVIVAL][i] = s.e[i] * P(f_e, f_e, E, h, false);
		fluxes[E_MU][i] = s.e[i] * P(f_e, f_m, E, h, false);
		fluxes[E_TAU][i] = s.e[i] * P(f_e, f_t, E, h, false);

		// and same for antiparticles
		fluxes[ANTIMU_SURVIVAL][i] = s.antimu[i] * P(f_m, f_m, E, h, true);
		fluxes[ANTIMU_ANTIE][i] = s.antimu[i] * P(f_m, f_e, E, h, true);
		fluxes[ANTIMU_ANTITAU][i] = s.antimu[i] * P(f_m, f_t, E, h, true);

		fluxes[ANTIE_SURVIVAL][i] = s.antie[i] * P(f_e, f_e, E, h, true);
		fluxes[ANTIE_ANTIMU][i] = s.antie[i] * P(f_e, f_m, E, h, true);
		fluxes[ANTIE_ANTITAU][i] = s.antie[i] * P(f_e, f_t, E, h, true);
	}

	// Now we normalize the fluxes to their respective FD predicted values
	// First we calculate integrals
	float integrals[N_FLUXES] = {0};
	for (int i=0; i<N_FLUXES; ++i) {
		// For the nu e antinu e background, aka the survival of nu e's and antinu e's,
		// we only have one normalization, thus they need to be normalized together
		if (i == E_SURVIVAL) {
			// normalize e survival and anti e survival
			// add the fluxes and normalize the sum
			float e_antie_survival[50] = {0};
			// array of fractions of the fluxes so we can retrieve them at the end
			float f[50];
			for (int j=0; j<50; ++j) {
				e_antie_survival[j] = fluxes[E_SURVIVAL][j] + fluxes[ANTIE_SURVIVAL][j];
				f[j] = fluxes[E_SURVIVAL][j] / fluxes[ANTIE_SURVIVAL][j];
			}
			float integral = 0.;
			for (int j=0; j<50; ++j) {
				integral += 0.2 * e_antie_survival[j];
			}
			for (int j=0; j<50; ++j) {
				e_antie_survival[j] /= integral;
			}
			// Now we have the normalized sum of fluxes but we can retrieve the individual
			// fluxes by using the following
			for (int j=0; j<50; ++j) {
				fluxes[E_SURVIVAL][j] = e_antie_survival[j] * f[j] / (f[j] + 1.);
				fluxes[ANTIE_SURVIVAL][j] = e_antie_survival[j] / f[j] / (1./f[j] + 1.);
			}
		} else if (i == ANTIE_SURVIVAL) {/* do nothing */} else {	
			for (int j=0; j<50; ++j) {
				integrals[i] += 0.2 * fluxes[i][j];
			}
			// normalize to 1
			for (int j=0; j<50; ++j) {
				fluxes[i][j] /= integrals[i];
			}
		}
	}
	// Then we multiply the fluxes by their respective normalizations
	for (int i=0; i<50; ++i) {
		fluxes[MU_E][i] *= e_signal;
		fluxes[ANTIMU_ANTIE][i] *= antie_signal;
		fluxes[MU_SURVIVAL][i] *= mu_signal;
		fluxes[ANTIMU_SURVIVAL][i] *= antimu_background;
	}

}

// Normalize the flux to the CDR predicted event rate (vol.2, p.27).
// Mode is 0 for neutrino mode, 1 for antineutrino mode
void normalize(initial_spectrum* s) {
	float flux_mu, flux_e, flux_antimu, flux_antie = 0;

	// Compute integrals
	for (int i=0; i<50; ++i) {
		// 0.2 is dE
		flux_mu += 0.2 * s->mu[i];
		flux_e += 0.2 * s->e[i];
		flux_antimu += 0.2 * s->antimu[i];
		flux_antie += 0.2 * s->antie[i];
	}

	// Now the integrals of all spectra are 1
	for (int i=0; i<50; ++i) {
		s->mu[i] /= flux_mu;
		s->e[i] /= flux_e;
		s->antimu[i] /= flux_antimu;
		s->antie[i] /= flux_antie;
	}
}
