#ifndef NUMATH
#define NUMATH
#include <complex>

double sinsq(double x);

double cossq(double x);

void dot(const double* mat, const double* vec, double* out, int dim = 3);

void mat_mult(const complex<double>* m1, const complex<double>* m2, complex<double>* out, int dim=3);

#endif
