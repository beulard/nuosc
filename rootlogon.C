
{
	gROOT->ProcessLine(".x helper.c");
	gROOT->ProcessLine(".L vector.c+");
	gROOT->ProcessLine(".L numath.c+");
	gROOT->ProcessLine(".L osc_defs.c+");
	gROOT->ProcessLine(".x spectrum.c");
}
