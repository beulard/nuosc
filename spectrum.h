#ifndef SPECTRUM
#define SPECTRUM
#include "osc_defs.h"

// global experiment baseline
double L = 1300.;


// Normalizations indices 
enum {
	// Oscillated nu e's
	E_SIGNAL,
	// Oscillated antinu e's
	ANTIE_SIGNAL,
	// Survival nu e + antinu e 's
	// Since these are grouped together, they need to be normalized together
	E_ANTIE_BACKGROUND,
	// Survival nu mu's
	MU_SIGNAL,
	// Survival antinu mu's
	ANTIMU_BACKGROUND,
	N_NORMS
};

// Flux indices
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

// Normalizations from CDR, in order given by enum above
const float nh_norm[N_NORMS] = { 861., 13., 159., 10842., 958. };
const float ih_norm[N_NORMS] = { 495., 26., 159., 10842., 958. };


struct spectrum {
	hierarchy* h;

	float fluxes[N_FLUXES][50];
	// predicted fluxes at the far detector
	float mu[50];
	float e[50];
	// tau muons can also appear at the FD
	float tau[50];

	// now we want to keep the antiparticles spectra
	float antimu[50];
	float antie[50];
	float antitau[50];
};

struct initial_spectrum {
	float mu[50];
	float antimu[50];
	float e[50];
	float antie[50];
};

// calculates probability of a neutrino flavor a to oscillate to flavor b,
// given an energy E, a hierarchy h and a constant baseline
double P(flavor a, flavor b, double E, hierarchy* h, bool anti);

// Oscillate neutrinos for each flavor and for each energy and get the FD flux.
void oscillate(const initial_spectrum* is, spectrum* s);

// Normalize the flux to the CDR predicted event rate (vol.2, p.27).
// Mode is 0 for neutrino mode, 1 for antineutrino mode
void normalize(initial_spectrum* s);

#endif
