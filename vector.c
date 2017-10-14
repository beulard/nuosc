#include "vector.h"
#include "TString.h"

//	print a vector's elements
void vec_print(const double* vec, int dim) {
	for(int i=0; i<dim; ++i) {
		Printf("%f", vec[i]);
	}
}

//	make a vector zero
void vec_zero(double* vec, int dim) {
	for(int i=0; i < dim; ++i) {
		vec[i] = 0.0;
	}
}
