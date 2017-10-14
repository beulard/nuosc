#ifndef VECTOR
#define VECTOR

//	print a vector's elements
static void vec_print(const double* vec, int dim = 3) {
	for(int i=0; i<dim; ++i) {
		Printf("%f", vec[i]);
	}
}

//	make a vector zero
static void vec_zero(double* vec, int dim = 3) {
	for(int i=0; i < dim; ++i) {
		vec[i] = 0.0;
	}
}

#endif
