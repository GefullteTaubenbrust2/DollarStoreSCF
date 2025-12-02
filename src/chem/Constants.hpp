#pragma once

namespace flo {
	#ifndef PI
	#define PI 3.14159265358979
	#endif

	// From CODATA Internationally recommended 2022 values of the fundamental physical constants
	// [1] P. J. Mohr, D. B. Newell, B. N. Taylor, CODATA recommended values of the fundamental physical constants: 2022, Rev. Mod. Phys. 97, 2025, 025002-1-62.
	const double avogadro = 6.02214076e+23;
	const double Ha_to_eV = 27.211386245989;
	const double Ha_to_J = 4.3597447222060e-18;
	const double Ha_to_kJmol = Ha_to_J * avogadro / 1000.0;
	const double a0_to_A = 0.529177210545;
	const double A_to_a0 = 1.889726126;
	const double a0_to_m = 5.29177210545e-11;
	const double Da_to_me = 1822.888487;
	const double au_to_s = 2.41888432658e-17;
	const double c0_metric = 299792458;
	const double kB_in_Ha_K = 8.617333262e-5 / Ha_to_eV;

	const double Pa_to_au = (a0_to_m * a0_to_m * a0_to_m) / Ha_to_J;
	const double au_to_per_cm = 0.005 / (PI * au_to_s * c0_metric);
}