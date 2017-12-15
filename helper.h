#ifndef HELPER
#define HELPER
#include <complex>
#include "TMath.h"

const double pi = TMath::Pi();

enum color_index {
	CI_BACKGROUND,
	CI_E,
	CI_MU,
	CI_TAU,
	CI_NH,
	CI_IH,
	CI_1,
	CI_2,
	CI_3,
	CI_N
};

int ci[CI_N];

static void complex_print(std::complex<double> z) {
	Printf("%f + %fi", std::real(z), std::imag(z));
}

#endif
