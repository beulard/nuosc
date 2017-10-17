//	neutrino oscillation parameters in the NORMAL HIERARCHY
#ifndef OSC_DEFS_NH
#define OSC_DEFS_NH

//	mass differences, in eV^2
const double dm2_21 = 7.37e-5;
//	absolute value of dm^2 as given in pdg neutrino mixing paper
const double dm2 = 2.50e-3;
const double dm2_31 = dm2 + dm2_21 / 2.;
const double dm2_32 = dm2 - dm2_21 / 2.;

//	mixing angles
const double t12 = TMath::ASin(TMath::Sqrt(0.297));
const double t23 = TMath::ASin(TMath::Sqrt(0.437));
const double t13 = TMath::ASin(TMath::Sqrt(0.0214));

//const double d_cp = TMath::Pi() * 1.35;
const double d_cp = 0.;

#endif
