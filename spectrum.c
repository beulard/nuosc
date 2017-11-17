#include "spectrum.h"

// calculates probability of a neutrino flavor a to oscillate to flavor b,
// given an energy E, a hierarchy h and a constant baseline
/*double P(flavor a, flavor b, float E, hierarchy* h, bool anti) {
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
}*/

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


// Oscillate neutrinos for each flavor and for each energy and get the FD flux.
void oscillate(const initial_spectrum* is, spectrum* os) {	
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
		// Instead of using the probability value at the discrete point given by i * 0.2,
		// we average the probabilities over the range (i*0.2 + eps, (i+1) * 0.2)
		float eps = 1e-6;
		int samples = 10;
		float avg_P[N_FLUXES] = {0};
		for (int j=0; j<samples; ++j) {
			float E = eps + 0.2 * (i + (float)j / (float)samples);

			avg_P[MU_SURVIVAL] += P(f_m, f_m, E, h, false);
			avg_P[MU_E] += P(f_m, f_e, E, h, false);
			avg_P[MU_TAU] += P(f_m, f_t, E, h, false);
			avg_P[E_SURVIVAL] += P(f_e, f_e, E, h, false);
			avg_P[E_MU] += P(f_e, f_m, E, h, false);
			avg_P[E_TAU] += P(f_e, f_t, E, h, false);
			avg_P[ANTIMU_SURVIVAL] += P(f_m, f_m, E, h, true);
			avg_P[ANTIMU_ANTIE] += P(f_m, f_e, E, h, true);
			avg_P[ANTIMU_ANTITAU] += P(f_m, f_t, E, h, true);
			avg_P[ANTIE_SURVIVAL] += P(f_e, f_e, E, h, true);
			avg_P[ANTIE_ANTIMU] += P(f_e, f_m, E, h, true);
			avg_P[ANTIE_ANTITAU] += P(f_e, f_t, E, h, true);
		}
		for (int j=0; j<N_FLUXES; ++j) {
			avg_P[j] /= (float)samples;
		}

		/*if(call == 1) {
			Printf("E = %.1f, P_avg(mu->mu) = %f", i*0.2, avg_P[MU_TAU]);
			Printf("E = %.1f, P(e->mu) = %f", i*0.2, P(f_m, f_t, i*0.2 + 1e-12, h, false));
			if(i < 6)
				Printf("P %f", avg_P[0]);
		}*/
			
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

		/*float p = 0.;
		if (i > 0)
			p = P_me(0.2 * i, h, true);
		fluxes[MU_SURVIVAL][i] = s.mu[i] * (1. - p);
		fluxes[MU_E][i] = s.mu[i] * p;
		if (p > 1.) {
			Printf("%e\t%.1f", p, 0.2 * i);
		}*/

		/*
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
		*/
	}

	// debug
	/*if (call==1) {
		for(int j=0; j<6; ++j) {
			//Printf("P", avg_P[MU_SURVIVAL][j]);
			Printf("%f %e", 0.2 * j, fluxes[MU_SURVIVAL][j]);
			Printf("%e %e", s.mu[j], is->mu[j]); 
		}
	}*/

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
				integral += e_antie_survival[j];
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
		} else if (i == ANTIE_SURVIVAL) {;} else {	
			for (int j=0; j<50; ++j) {
				integrals[i] += fluxes[i][j];
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
