#include "numath.h"
#include "vector.h"
#include "TMath.h"
#include <complex>

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

void mat_mult(const complex<double>* m1, const complex<double>* m2, complex<double>* out, int dim=3) {
	for(int i=0; i<dim; ++i) {
		for(int j=0; j<dim; ++j) {
			out[i*3 + j] = 0.;
			for(int k=0; k<dim; ++k) {
				out[i*3 + j] += m1[i*3 + k] * m2[j + k*3];
			}
		}
	}
}		
