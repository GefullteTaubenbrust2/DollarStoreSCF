#include "Element.hpp"

namespace flo {
	const char* element_symbols[118] = {
		"H",                                                                                 "He",
		"Li","Be",                                                   "B", "C", "N", "O", "F","Ne",
		"Na","Mg",                                                  "Al","Si", "P", "S","Cl","Ar",
		 "K","Ca","Sc","Ti", "V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
		"Rb","Sr", "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te", "I","Xe",
		"Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
				  "Lu","Hf","Ta", "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
		"Fr","Ra","Ac","Th","Pa", "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
				  "Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
	};

	const std::string element_names[118] = {
		"hydrogen",
		"helium",
		"lithium",
		"beryllium",
		"boron",
		"carbon",
		"nitrogen",
		"oxygen",
		"fluorine",
		"neon",
		"sodium",
		"magnesium",
		"aluminium",
		"silicon",
		"phosphorus",
		"sulfur",
		"chlorine",
		"argon",
		"potassium",
		"calcium",
		"scandium",
		"titanium",
		"vanadium",
		"chromium",
		"manganese",
		"iron",
		"cobalt",
		"nickel",
		"copper",
		"zinc",
		"gallium",
		"germanium",
		"arsenic",
		"selenium",
		"bromine",
		"krypton",
		"rubidium",
		"strontium",
		"yttrium",
		"zirconium",
		"niobium",
		"molybdenum",
		"technetium",
		"ruthenium",
		"rhodium",
		"palladium",
		"silver",
		"cadmium",
		"indium",
		"tin",
		"antimony",
		"tellurium",
		"iodine",
		"xenon",
		"caesium",
		"barium",
		"lanthanum",
		"cerium",
		"praeseodymium",
		"neodymium",
		"promethium",
		"samarium",
		"europium",
		"gadolinum",
		"terbium",
		"dysprosium",
		"holmium",
		"erbium",
		"thulium",
		"ytterbium",
		"lutetium",
		"hafnium",
		"tantalum",
		"tungsten",
		"rhenium",
		"osmium",
		"irdium",
		"platinum",
		"gold",
		"mercury",
		"thallium",
		"lead",
		"bismuth",
		"polonium",
		"astatine",
		"radon",
		"francium",
		"radium",
		"actinium",
		"thorium",
		"protactinium",
		"uranium",
		"neptunium",
		"plutonium",
		"americium",
		"curium",
		"berkelium",
		"californium",
		"einsteinium",
		"fermium",
		"mendelevium",
		"nobellium",
		"lawrencium",
		"rutherfordium",
		"dubnium",
		"seaborgium",
		"bohrium",
		"hassium",
		"meitnerium",
		"darmstadtium",
		"roentgenium",
		"copernicium",
		"nihonium",
		"flerovium",
		"moscovium",
		"livermorium",
		"tennessine",
		"oganesson"
	};

	// [1] J. S. Coursey, D. J. Schwab, J. J. Tsai, and R. A. Dragoset
	// NIST Physical Measurement Laboratory
	// 
	// From the website:
	// 
	// "
	// The atomic weights are available for elements 1 through 118 and isotopic compositions or abundances are given when appropriate. 
	// The atomic weights data were published by 
	// [2] J. Meija et al in Atomic Weights of the Elements 2013, 
	// and the isotopic compositions data were published by 
	// [3] M. Berglund and M.E. Wieser in Isotopic Compositions of the Elements 2009. 
	// The relative atomic masses of the isotopes data were published by 
	// [4] M. Wang, G. Audi, A.H. Wapstra, F.G. Kondev, M. MacCormick, X. Xu1, and B. Pfeiffer in The AME2012 Atomic Mass Evaluation
	// "
	// 
	// Mean values are used. For elements beyond Pu, the most stable known nuclei are assumed and non-exact masses used.
	const double common_atomic_masses[118] = {
		  1.008,                                                                                                                                                   4.002,
		  6.968,  9.012,                                                                                             10.814,  12.011,  14.007,  15.999,  18.998,  20.180,
		 22.990,  24.306,                                                                                            26.982,  28.085,  30.974,  32.069,  35.452,  39.948,
		 39.098,  40.078,  49.956,  47.867,  50.942,  51.996,  54.938,  55.845,  58.933,  58.693,  65.546,  65.382,  69.723,  72.631,  74.922,  78.972,  79.904,  83.798,
		 85.468,  87.621,  88.906,  91.224,  92.906,  98.951,  98.000, 101.072, 102.906, 106.421, 107.868, 112.414, 114.818, 118.711, 121.760, 127.603, 126.904, 131.294,
		132.905, 137.328, 138.905, 140.116, 140.908, 144.242, 145.000, 150.362, 151.964, 157.253, 158.925, 162.500, 164.930, 167.259, 168.934, 173.054,
						  174.967, 178.492, 180.948, 183.341, 186.207, 190.233, 192.217, 195.085, 196.967, 200.592, 204.384, 207.210, 208.980, 209.000, 210.000, 222.000,
		233.000, 226.000, 227.000, 232.038, 231.036, 238.029, 237.000, 244.000, 243.000, 247.000, 247.000, 251.000, 252.000, 257.000, 258.000, 259.000,
						  266.000, 267.000, 268.000, 269.000, 278.000, 278.000, 282.000, 282.000, 286.000, 286.000, 286.000, 290.000, 290.000, 293.000, 294.000, 294.000
	};

	// [5] V. Pekka Pyykko, M. Atsumi, Chem. Eur. J. 15, 2008, 186-197.
	const double covalent_radii_A[118] = {
		0.32,                                                                                                 0.46,
		1.33, 1.02,                                                             0.85, 0.75, 0.71, 0.63, 0.64, 0.67,
		1.55, 1.39,                                                             1.26, 1.16, 1.11, 1.03, 0.99, 0.96,
		1.96, 1.71, 1.48, 1.36, 1.34, 1.22, 1.19, 1.16, 1.11, 1.10, 1.12, 1.18, 1.24, 1.21, 1.21, 1.16, 1.14, 1.17,
		2.10, 1.85, 1.63, 1.54, 1.47, 1.38, 1.28, 1.25, 1.25, 1.20, 1.28, 1.36, 1.42, 1.40, 1.40, 1.36, 1.33, 1.31,
		2.32, 1.96, 1.80, 1.63, 1.76, 1.74, 1.73, 1.72, 1.68, 1.69, 1.68, 1.67, 1.66, 1.65, 1.64, 1.70,
					1.62, 1.52, 1.46, 1.37, 1.31, 1.29, 1.22, 1.23, 1.24, 1.33, 1.44, 1.44, 1.51, 1.45, 1.47, 1.42,
		2.23, 2.01, 1.86, 1.75, 1.69, 1.70, 1.71, 1.72, 1.66, 1.66, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76,
					1.61, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28, 1.21, 1.22, 1.36, 1.43, 1.62, 1.75, 1.65, 1.57
	};

	int getElement(const std::string& symbol) {
		if (symbol.size() && symbol.size() < 3) {
			for (int i = 0; i < 118; ++i) {
				if (element_symbols[i] == symbol) return i + 1;
			}
			return 0;
		}
		else return 0;
	}
}