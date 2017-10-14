#include "numath.h"
#include "vector.h"
#include "TMath.h"

double sinsq(double x) {
	return TMath::Power(TMath::Sin(x), 2);
}

double cossq(double x) {
	return TMath::Power(TMath::Cos(x), 2);
}

//	matrix * vector multiplication
void dot(const double* mat, const double* vec, double* out, int dim) {
	vec_zero(out);

	for(int i=0; i < dim; ++i) {
		for(int j=0; j < dim; ++j) {
			out[i] += mat[j + i * dim] * vec[j];
		}
	}
}
