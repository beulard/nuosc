#ifndef SPECTRUM
#define SPECTRUM
#include "osc_defs.h"

// global experiment baseline
double L = 1300.;

// Normalizations
// Oscillated nu e's
float e_signal = 861.;
// Oscillated antinu e's
float antie_signal = 13.;
// Survival nu e + antinu e 's
// Since these are grouped together, they need to be normalized together
float e_antie_background = 159.;

// Survival nu mu's
float mu_signal = 10842.;
// Survival antinu mu's
float antimu_background = 958.;


struct spectrum {
	hierarchy h;

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
double P(flavor a, flavor b, float E, hierarchy* h, bool anti);

// Propagate neutrinos for each flavor and for each energy and get the FD flux.
void propagate(const initial_spectrum* is, spectrum* s);

// Normalize the flux to the CDR predicted event rate (vol.2, p.27).
// Mode is 0 for neutrino mode, 1 for antineutrino mode
void normalize(initial_spectrum* s);

#endif
