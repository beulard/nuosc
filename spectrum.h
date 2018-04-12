#ifndef SPECTRUM
#define SPECTRUM
#include "osc_defs.h"

// TODO change this. We must allow different number of bins and different firstbin
// and different energy range etc.
//const int Nbins = 37;
//const int firstbin = 3;


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



class initial_spectrum {
public:
	//float mu[50];
	//float antimu[50];
	//float e[50];
	//float antie[50];


	int N_entries;
	double* mu_flux;
	double* E;
	ROOT::Math::Interpolator *interp;
		

	// Reads a spectrum from a csv file. Only fills the 'mu' array.
	void read(const char* file);
	void plot();
};


// calculates probability of a neutrino flavor a to oscillate to flavor b,
// given an energy E, a hierarchy h and a constant baseline
double P(flavor a, flavor b, double E, double L, parameters* p, bool anti);

// Same, including matter effects under constant density assumption
double P_me(flavor a, flavor b, double E, double L, parameters* p, bool anti);

class spectrum {
public:
	spectrum();
	~spectrum();

	parameters* p;
	//double events[N_RATES][Nbins];
	// Electron signal array
	int N_bins;
	double* e;


	// Oscillate neutrinos for each flavor and for each energy and get the FD flux.
	//void oscillate(const initial_spectrum* is, double baseline,
	//				bool normal=true, bool antimode=false);
	// CHANGED we need to be able to control the energy range and bin count for oscillations
	// Oscillate neutrinos in range [x1, x2] given baseline and bin_count
	// DROPPED we will not need this if we don't model another experiment (eg icecube)
	void oscillate(const initial_spectrum* is, double baseline);
				   //double x1, double x2, int bin_count);
	
	// Energy reconstruction function with smearing of the reconstructed energy to
	// simulate experimental uncertainties.
	// If the 'smooth' parameter is set, the spectrum is reconstructed multiple times
	// and averaged over to make it smooth.
	// Reconstructs into another spectrum, in case we want to conserve the current one.
	void reconstruct(spectrum* os, int smooth=1);
	// Calls reconstruct(s, s, smooth);
	void reconstruct(int smooth=1);

	// Get total number of events in E_SIGNAL event rate
	double get_integral();
	
	// Normalize a spectrum to the CDR predicted event rates, given the d_cp=0
	// spectrum integral that we want to normalize to. Usually that will be the
	// spectrum of the same mass hierarchy.
	// CHANGED TO FOLLOWING TO GENERALIZE TO HYPERK
	// Normalize to a given predicted event rate given the integral of the d_cp=0 spectrum
	void normalize(double integral, double predicted);

};

#endif
