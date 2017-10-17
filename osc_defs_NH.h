//	neutrino oscillation parameters in the NORMAL HIERARCHY
#ifndef OSC_DEFS_NH
#define OSC_DEFS_NH

//	mass differences, in eV^2
const double dm2_21 = 7.37e-5;
//	absolute value of dm^2 as given in pdg neutrino mixing paper
const double dm2 = 2.50e-3;
const double dm2_31 = dm2 + dm2_21 / 2.;
const double dm2_32 = dm2 - dm2_21 / 2.;

// mass difference matrix, for indexing
const double dm2_mat[] = { 0.0,    dm2_21, dm2_31,
					   dm2_21, 0.0,    dm2_32,
					   dm2_31, dm2_32, 0.0  };

//	mixing angles
const double t12 = TMath::ASin(TMath::Sqrt(0.297));
const double t23 = TMath::ASin(TMath::Sqrt(0.437));
const double t13 = TMath::ASin(TMath::Sqrt(0.0214));

// trigonometric functions
const double s12 = TMath::Sin(t12);
const double s23 = TMath::Sin(t23);
const double s13 = TMath::Sin(t13);

const double c12 = TMath::Cos(t12);
const double c23 = TMath::Cos(t23);
const double c13 = TMath::Cos(t13);

const double d_cp = TMath::Pi() * 1.35;
//const double d_cp = 0.;

enum flavor {
	f_e = 0,
	f_m = 1,
	f_t = 2
};

#endif
