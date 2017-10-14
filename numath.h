#ifndef NUMATH
#define NUMATH

static double sinsq(double x) {
	return TMath::Power(TMath::Sin(x), 2);
}

static double cossq(double x) {
	return TMath::Power(TMath::Cos(x), 2);
}

//	matrix * vector multiplication
static void dot(const double* mat, const double* vec, double* out, int dim = 3) {
	vec_zero(out);

	for(int i=0; i < dim; ++i) {
		for(int j=0; j < dim; ++j) {
			out[i] += mat[j + i * dim] * vec[j];
		}
	}
}

#endif
