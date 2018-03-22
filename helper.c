#include "TColor.h"
#include "helper.h"



void complex_print(std::complex<double> z) {
	Printf("%f + %fi", std::real(z), std::imag(z));
}

double mean(double* vals, int N) {
	double r = 0.;
	for (int i=0; i<N; ++i) {
		r += vals[i] / N;
	}
	return r;
}

double std_dev(double* vals, double mean, int N) {
	double r = 0.;
	for (int i=0; i<N; ++i) {
		r += pow(vals[i] - mean, 2) / (N-1);
	}

	return sqrt(r);
}

