#ifndef HELPER
#define HELPER
#include <complex>
#include "TMath.h"
#include "color_index.h"

const double pi = TMath::Pi();


// Print complex number
void complex_print(std::complex<double> z);

// Statistical mean of a set of N values
double mean(double* vals, int N);

// Statistical standard deviation from the mean
double std_dev(double* vals, double mean, int N);

#endif
