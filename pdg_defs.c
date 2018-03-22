#include "osc_defs.h"
#include "helper.h"

// Initializes pdg values
void pdg_defs() {
	pdg[NH].dm2_21 = 7.37e-5;
	pdg[NH].dm2_31 = 2.56e-3;
	pdg[NH].dm2_32 = pdg[NH].dm2_31 - pdg[NH].dm2_21;

	pdg[NH].t12 = asin(sqrt(0.297));
	pdg[NH].t23 = asin(sqrt(0.425));
	pdg[NH].t13 = asin(sqrt(0.0215));

	pdg[NH].d_cp = pi * 1.38;

	pdg[NH].populate_common();


	pdg[IH].dm2_21 = 7.37e-5;
	pdg[IH].dm2_32 = -2.54e-3;
	pdg[IH].dm2_31 = pdg[IH].dm2_21 + pdg[IH].dm2_32;

	pdg[IH].t12 = asin(sqrt(0.297));
	pdg[IH].t23 = asin(sqrt(0.589));
	pdg[IH].t13 = asin(sqrt(0.0216));

	pdg[IH].d_cp = pi * 1.31;

	pdg[IH].populate_common();


	// Values of the 1-sigma error from the PDG review
	pdg_sd[NH].dm2_21 = 0.20e-5;
	pdg_sd[NH].dm2_31 = 0.04e-3;
	pdg_sd[NH].t12 = 0.021;
	pdg_sd[NH].t23 = 0.064;
	pdg_sd[NH].t13 = 0.00275;

	pdg_sd[IH].dm2_21 = 0.20e-5;
	pdg_sd[IH].dm2_31 = 0.04e-3;
	pdg_sd[IH].t12 = 0.021;
	pdg_sd[IH].t23 = 0.069;
	pdg_sd[IH].t13 = 0.00275;
}

