#include "spectrum.h"

// calculates probability of a neutrino flavor a to oscillate to flavor b,
// given an energy E, a hierarchy h and a constant baseline
double P(flavor a, flavor b, double E, double L, parameters* p, bool anti) {
	double r = 0.;

	for(int i=0; i<3; ++i) {
		r += pow(abs(p->MNS[a*3 + i] * conj(p->MNS[b*3 + i])), 2);
	}

	if (anti) {
		for (int i=0; i<3; ++i) {
			for (int j=0; j<3; ++j) {
				if (j > i) {
					r += 2. * real( conj(p->MNS[a*3 + i]) * p->MNS[a*3 + j]
						  * p->MNS[b*3 + i] * conj(p->MNS[b*3 + j])
						  * exp(complex<double>(-2.i * 1.2668 * p->dm2_mat[i*3 + j] * L/E)));
				}
			}
		}
	}
	else {
		for(int i=0; i<3; ++i) {
			for(int j=0; j<3; ++j) {
				if(j > i) {
					r += 2. * real(p->MNS[a*3 + i] * conj(p->MNS[a*3 + j]) 
							* conj(p->MNS[b*3 + i]) * p->MNS[b*3 + j]
							* exp(complex<double>(-2.i * 1.2668 *  p->dm2_mat[i*3 + j] * L/E)));
				}
			}
		}
	}
	return r;
}


// Electron density in mantle in units cm^-3 / (Avogadro's number)
const double Ne = 2.2;
// Oscillation probability for nu_mu -> nu_e or antinu_mu -> antinu_e
// including matter effects (Boris Kayser Neutrino Oscillation phenomenology p. 13)
double P_me(flavor fa, flavor fb, double E, double L, parameters* p, bool anti) {

	double x = 1.53e-4 * Ne * E / p->dm2_31;
	double D = 1.269 * p->dm2_31 * L / E;
	double a = p->dm2_21 / p->dm2_31;
	x = anti ? -x : x;
	double d_cp = anti ? -p->d_cp : p->d_cp;

	double T1 = sinsq(p->t23) * sinsq((1. - x) * D) / pow(1. - x, 2);
	double T2 = sin(d_cp) * sin(2 * p->t12) * sin(2 * p->t23) * sin(D)
				* sin(x * D) / x * sin((1 - x) * D) / (1 - x);
	double T3 = cos(d_cp) * sin(2 * p->t12) * sin(2 * p->t23) * cos(D)
				* sin(x * D) / x * sin((1 - x) * D) / (1 - x);
	double T4 = cossq(p->t23) * sinsq(2 * p->t12) * sinsq(x * D) / pow(x, 2);

	double r = sinsq(2 * p->t13) * T1 - a * sin(2 * p->t13) * T2 + a * sinsq(2 * p->t13) * T3
			+ pow(a, 2) * T4;

	return r;
}

// CDR transition probability
double P_me2(flavor fa, flavor fb, double E, double L, parameters* p, bool anti) {
	double aL = 1.93e-4 * Ne * L;
	aL = anti ? -aL : aL;

	double d_cp = anti ? -p->d_cp : p->d_cp;
	double D31 = 1.269 * p->dm2_31 * L / E;
	double D21 = 1.269 * p->dm2_21 * L / E;

	double r = sinsq(p->t23) * sinsq(2*p->t13) * sinsq(D31 - aL) / pow(D31 - aL, 2) * pow(D31, 2)

		+ sin(2*p->t23) * sin(2 * p->t13) * sin(2*p->t12) * sin(D31 - aL) / (D31 - aL) * D31 * sin(aL) / aL * D21 * cos(D31 + d_cp)

		+ cossq(p->t23) * sinsq(2*p->t12) * sinsq(aL) / pow(aL, 2) * pow(D21, 2);

	return r;
}


void initial_spectrum::read(const char* file) {
	TTree* spectrum = new TTree("spectrum", "");

	// Read file with the descriptor given in the first line and delimiter ','
	spectrum->ReadFile(file, "", ',');

	TTreeReader r(spectrum);
	N_entries = spectrum->GetEntries();

	mu_flux = new double[N_entries];
	E = new double[N_entries];
	interp = new ROOT::Math::Interpolator(N_entries, ROOT::Math::Interpolation::kCSPLINE);

	TTreeReaderValue<double> r_E(r, "E");
	TTreeReaderValue<double> r_mu(r, "mu");

	int i=0;
	while(r.Next()) {
		E[i] = *r_E;	
		mu_flux[i] = *r_mu;
		i++;
	}

	interp->SetData(N_entries, E, mu_flux);

	spectrum->Delete();
}

void initial_spectrum::plot() {
	TCanvas* c = new TCanvas();
	c->SetLogy();
	
	//c->SetLogx();
	
	
	const int N = 100;
	double x[N];
	double y[N];

	for (int i=0; i<N; ++i) {
		//x[i] = pow(10, i * (log10(E[N_entries-1]) - log10(E[0])) / N + log10(E[0]));
		x[i] = i * (E[N_entries-1] - E[0]) / N + E[0];
		y[i] = interp->Eval(x[i]);
	}

	TGraph* g = new TGraph(N, x, y);
	//g->SetMaximum(1e3);
	//g->SetMinimum(1e0);

	g->Draw();
}

double spectrum::get_integral() {
	double r=0.;
	for (int i=0; i<Nbins; ++i) {
		r += events[E_SIGNAL][i];
	}
	return r;
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
//void spectrum::oscillate(const initial_spectrum* is, double L, bool normal, bool antimode) {	
void spectrum::oscillate(const initial_spectrum* is, double L,
					     double x1, double x2, int bin_count) {

	e = new double[bin_count];
	for (int i=0; i<bin_count; ++i) {
		// Instead of using the probability value at the discrete point given by i * 0.2,
		// we average the probabilities over the range (i*0.2 + eps, (i+1) * 0.2)
		int samples = sqrt(bin_count);

		double avg_P[N_RATES] = {0};
		double avg_flux = 0;
		for (int j=0; j<samples; ++j) {
			double E = 0.2 * (firstbin + i + (double)j / (double)samples);
			
			avg_flux += is->interp->Eval(E) / samples;
			
			
			avg_P[E_SIGNAL] += P_me(f_m, f_e, E, L, p, false) / samples;
			avg_P[ANTIE_SIGNAL] += P_me(f_m, f_e, E, L, p, true) / samples;

			avg_P[MU_SIGNAL] += P(f_m, f_m, E, L, p, false) / samples;

		}
		//for (int j=0; j<N_RATES; ++j) {
	//		avg_P[j] /= (double)samples;
	//	}

			

		//	Number of events at FD if we could detect every single neutrino.
		//events[E_SIGNAL][i] = is->mu[i + firstbin] * avg_P[E_SIGNAL];
		//events[ANTIE_SIGNAL][i] = is->antimu[i + firstbin] * avg_P[ANTIE_SIGNAL];
		//events[MU_SIGNAL][i] = is->mu[i + firstbin] * avg_P[MU_SIGNAL];
		
		events[E_SIGNAL][i] = avg_flux * avg_P[E_SIGNAL];

	}

	// TODO
	// We are missing events as we are reducing the range of energies without taking
	// into account the neutrinos produced outside of that range that may have been 
	// reconstructed into the range.


}

// Energy reconstruction function with smearing of the reconstructed energy to
// simulate experimental uncertainties
void spectrum::reconstruct(spectrum* os, int smooth) {

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
			for (int j=0; j<ceil(events[E_SIGNAL][i]); ++j) {
				// Can play around with the sigma on the gaussian for different results
				double E_rec = -1.;
				while (E_rec > E + 2. || E_rec < 0.1 * E || E_rec < 0.6) {
					E_rec = ra.Gaus(E, 0.1 * E);
				}
				//E_rec = E;
				

				// Find the bin corresponding to the reconstructed energy
				int bin = round((E_rec - 0.7) / 7.4 * Nbins + 1e-3);

		
				// The value of each event is 1, except for the last one which might be
				// worth less because of our floating point precision
				double value = 1.;
				if (j == ceil(events[E_SIGNAL][i]) - 1){
					value = events[E_SIGNAL][i] - floor(events[E_SIGNAL][i]);
				}
				// We want to move this event to the bin we picked
				s_rec.events[E_SIGNAL][bin] += value;
			}
			// Avoid zero events for stability (crashes in delta chi^2 calculation)
			if (s_rec.events[E_SIGNAL][i] < 1)
				s_rec.events[E_SIGNAL][i] = 1.;
		}
		
		for (int i=0; i<Nbins; ++i) {
			s_f.events[E_SIGNAL][i] += s_rec.events[E_SIGNAL][i] / N_rec;
		}

	}
	// After obtaining final spectrum, copy it over to output
	memcpy(os->events, s_f.events, Nbins * N_RATES * sizeof(double));
}

void spectrum::reconstruct(int smooth) {
	reconstruct(this, smooth);
}

void spectrum::normalize(double integral, double predicted) {
	for (int i=0; i<Nbins; ++i) {
		//events[E_SIGNAL][i] *= norms[h][E_SIGNAL] / integral;
		events[E_SIGNAL][i] *= predicted / integral;
	}
}
