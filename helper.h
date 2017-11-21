#ifndef HELPER
#define HELPER
#include <complex>

const double pi = TMath::Pi();

static void complex_print(std::complex<double> z) {
	Printf("%f + %fi", std::real(z), std::imag(z));
}

#endif
