
{
	gROOT->ProcessLine(".L helper.c");
	gROOT->ProcessLine(".L vector.c+");
	gROOT->ProcessLine(".L numath.c+");
	gROOT->ProcessLine(".L osc_defs.c");
	//gROOT->ProcessLine(".L spectrum.c");
	gROOT->ProcessLine(".x color_index.c");
	gROOT->ProcessLine(".x pdg_defs.c");
}
