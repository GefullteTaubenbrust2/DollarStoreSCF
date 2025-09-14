#include "Damping.hpp"

#include "SCFCommon.hpp"
#include "Population.hpp"

using namespace flo;

namespace scf {
	double damping_quantity = 0.0;

	// M. C. Zerner, M. Hehenberger, A dynamical damping scheme for converging molecular scf calculations, Chem. Phys. Lett. 62 (3), 550-554 (1979).
	double getZernerHehenbergerDamping() {
		calculateMullikenPopulations();

		uint atom_index = 0;
		uint weight = 0;
		double gross_population = 0.0;
		double quantity = 0.0;
		for (int i = 0;; ++i) {
			bool append = i == matrix_size;
			if (!append) append |= atom_index != basis_atoms[i];
			if (append) {
				quantity += gross_population * (double)weight;

				weight = 0;
				gross_population = 0.0;

				if (i == matrix_size) break;

				atom_index = basis_atoms[i];
			}

			gross_population += mulliken_populations[0].at(i, i);
			if (spin_treatment == SpinTreatment::unrestricted) gross_population += mulliken_populations[1].at(i, i);
			++weight;
		}

		quantity /= (double)matrix_size;

		double damping = quantity / (quantity - damping_quantity);
		if (damping < 0.0 || damping >= 1.0 || isnan(damping)) damping = 0.0;

		damping_quantity = damping;

		return damping;
	}
}