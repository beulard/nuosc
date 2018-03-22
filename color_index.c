#include "color_index.h"

// Initializes the global color index array
void color_index() {
	for (int i=0; i<CI_N; ++i) {
		ci[i] = 1789 + i;
	}
	new TColor(ci[CI_BACKGROUND], 250. / 255., 250. / 255, 250. / 255.);
	new TColor(ci[CI_E], 96. / 255., 149. / 255., 201. / 255);
	new TColor(ci[CI_MU], 205. / 255., 102. / 255., 95. / 255.);
	new TColor(ci[CI_TAU], 170. / 255., 196. / 255., 108. / 255);
	new TColor(ci[CI_NH], 50. / 255., 142. / 255., 237. / 255.);
	new TColor(ci[CI_IH], 239. / 255., 83. / 255., 138. / 255.);
	new TColor(ci[CI_1], 239. / 255., 84. / 255., 138. / 255.);
	new TColor(ci[CI_2], 50. / 255., 142. / 255., 237. / 255.);
	new TColor(ci[CI_3], 134. / 255., 211. / 255., 27. / 255.);
	new TColor(ci[CI_4], 100. / 255., 70. / 255., 173. / 255.);
}
