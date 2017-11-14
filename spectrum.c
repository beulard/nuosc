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
					p += 2. * real(h->MNS[a*3 + i] * conj(h->MNS[a*3 + j]) 
							* conj(h->MNS[b*3 + i]) * h->MNS[b*3 + j]
							* exp(complex<double>(-2.i * 1.2668 *  h->dm2_mat[i*3 + j] * L/E)));
				}
			}
		}
	}
	return p;
}

// Oscillation probability for nu_mu -> nu_e or antinu_mu -> antinu_e
// including matter effects (CDR vol. 2 p.19)
float P_me(float E, hierarchy* h, bool anti) {
	float d_31 = 1.2669 * h->dm2_mat[2 * 3 + 0] * L / E;
	float d_21 = 1.2669 * h->dm2_mat[1 * 3 + 0] * L / E; 
	float a = anti ? -4.255e-4 : 4.255e-4;
	float d_cp = anti ? -h->d_cp : h->d_cp;

	float p = pow(sin(h->t23) * sin(2*h->t13) * sin(d_31 - a*L) / (d_31 - a*L) * d_31, 2) 
			 
			 + sin(2*h->t23) * sin(2*h->t13) * sin(2*h->t12)
			 * sin(d_31 - a*L) / (d_31 - a*L) * d_31 
			 * sin(a*L) / (a*L) * d_21 * cos(d_31 + d_cp)

			 + pow(cos(h->t23) * sin(2*h->t12) * sin(a*L) / (a*L) * d_21, 2);

	return p;
}


// Propagate neutrinos for each flavor and for each energy and get the FD flux.
void propagate(const initial_spectrum* is, spectrum* os) {	
	hierarchy* h = os->h;
	
	initial_spectrum s;
	memcpy(&s, is, sizeof(initial_spectrum));

	// Propagated fluxes	
	float fluxes[N_FLUXES][50] = {0};
	const float* norms;
	if (h->type == NH) {
		norms = nh_norm;
	} else if (h->type == IH) {
		norms = nh_norm;
	}

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

		// Test for matter effects
		//float N_mu_e = s.mu[i] * P_me(E, h, false);
		//float N_anti_mu_e = s.antimu[i] * P_me(E, h, true);

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

				// And normalize to data
				fluxes[E_SURVIVAL][j] *= norms[E_ANTIE_BACKGROUND];
				fluxes[ANTIE_SURVIVAL][j] *= norms[E_ANTIE_BACKGROUND];
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
