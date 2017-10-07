using namespace TMath;

//	mass differences, in eV^2
const double dm2_12 = 7.37e-5;
const double dm2_23 = 2.50e-3;
const double dm2_13 = 2.50e-3;
// mass difference matrix, for indexing
const double dm2[] = { 0.0,    dm2_12, dm2_13,
					   dm2_12, 0.0,    dm2_23,
					   dm2_13, dm2_23, 0.0  };

//	mixing angles
const double t12 = ASin(Sqrt(0.297));
const double t23 = ASin(Sqrt(0.569));
const double t13 = ASin(Sqrt(0.0218));

// trigonometric functions
const double s12 = Sin(t12);
const double s23 = Sin(t23);
const double s13 = Sin(t13);

const double c12 = Cos(t12);
const double c23 = Cos(t23);
const double c13 = Cos(t13);

//	make a vector zero
void zero(double* vec, int dim = 3) {
	for(int i=0; i < dim; ++i) {
		vec[i] = 0.0;
	}
}

double sinsq(double x) {
	return Power(Sin(x), 2);
}

//	matrix * vector multiplication
void dot(const double* mat, const double* vec, double* out, int dim = 3) {
	zero(out);

	for(int i=0; i < dim; ++i) {
		for(int j=0; j < dim; ++j) {
			out[i] += mat[j + i * dim] * vec[j];
		}
	}
}

//	print a vector's elements
void print_vec(const double* vec, int dim = 3) {
	for(int i=0; i<dim; ++i) {
		Printf("%f", vec[i]);
	}
}

//	MNS matrix
const double MNS[] = {  c12 * c13, 					  s12 * c13, 				   s13,
					   -s12 * c23 - c12 * s23 * s13,  c12 * c23 - s12 * s23 * s13, s23 * c13,
					    s12 * s23 - c12 * s23 * s13, -c12 * s23 - s12 * c23 * s13, c23 * c13 };

enum flavor {
	f_e = 0,
	f_m = 1,
	f_t = 2
};

//	probability of a neutrino initially in flavor a to transition to flavor b, as a function
//	of L/E (km/GeV)
double P(flavor a, flavor b, double x) {
	double delta = a == b ? 1.0 : 0.0;
	double p = delta;

	for(int i=0; i < 3; ++i) {
		for(int j=0; j < 3; ++j) {
			if(i > j) {
				Printf("%d %d: %f", i, j, dm2[i*3+j]);
				//Printf("%e", 4.0 * MNS[a*3 + i] * MNS[b*3 + i] * MNS[a*3 + j] * MNS[b*3 + j]);
				p -= 4.0 * MNS[a*3 + i] * MNS[b*3 + i] * MNS[a*3 + j] * MNS[b*3 + j] 
						 * sinsq(1.27 * dm2[i*3 + j] * x);
			}
		}
	}
	
	return p;	
}

double plotP(double* x, double* par) {
	return P((flavor)par[0], (flavor)par[1], x[0]);
}


// main
void nu() {

	double flavor[] = { 0., 0., 1. };
	double mass[3];

	// turn the flavor eigenstate into mass eigenstate
	dot(MNS, flavor, mass);
	print_vec(mass);
	
	TF1* e_e = new TF1("e_e", plotP, 0., 50000., 2);
	e_e->SetParameters(f_e, f_e);
	TF1* e_m = new TF1("e_m", plotP, 0., 50000., 2);
	e_m->SetParameters(f_e, f_m);
	TF1* e_t = new TF1("e_t", plotP, 0., 50000, 2);
	e_t->SetParameters(f_e, f_t);
	e_e->SetNpx(1000);
	e_e->SetMinimum(0.0);
	e_e->SetMaximum(1.0);
	e_t->SetNpx(10000);
	e_m->SetNpx(10000);

	//e_t->SetLineColor(4);
	//e_m->SetLineColor(7);

	P(f_e, f_e, 50.);
	e_e->Draw();
	//e_m->Draw("same");
	//e_t->Draw("same");

	//TCanvas* c2 = new TCanvas();
	//TF1* t_t = new TF1("t_t", plotP, 0., 50000., 2);
	//t_t->SetParameters(f_t, f_t);
	//t_t->SetNpx(10000);
	//t_t->SetMinimum(-0.5);
	//t_t->SetMaximum(1.0);
	//t_t->Draw();
}



