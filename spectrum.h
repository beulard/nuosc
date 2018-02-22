#ifndef SPECTRUM
#define SPECTRUM
#include "osc_defs.h"

// global experiment baseline
const double L = 1300.;
const int Nbins = 37;
const int firstbin = 3;


// Event rates indices
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
	N_RATES
};


// Normalizations from CDR, in order given by enum above
// TODO we also need the antineutrino mode normalizations
const double norms[2][N_RATES] = { { 861., 13., 159., 10842., 958. },
								   { 495., 26., 159., 10842., 958. } };


struct spectrum {
	hierarchy* h;

	double events[N_RATES][Nbins];
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
void oscillate(const initial_spectrum* is, spectrum* s, 
				bool normal=true, bool antimode=false);

// Energy reconstruction function with smearing of the reconstructed energy to
// simulate experimental uncertainties
// If the 'smooth' parameter is set, the spectrum is reconstructed multiple times
// and averaged over to make it smooth
void reconstruct(const spectrum* s, spectrum* os, int smooth=1);

// Normalize a spectrum to the CDR predicted event rates, given the d_cp=0
// spectrum integral that we want to normalize to. Usually that will be the
// spectrum of the same mass hierarchy.
void normalize(spectrum* s, double integral, h_type h);


#endif
