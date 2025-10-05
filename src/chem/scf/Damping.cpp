#include "Damping.hpp"

#include "SCFCommon.hpp"
#include "Population.hpp"

#include "../Molecule.hpp"

using namespace flo;

namespace scf {
	std::vector<double> populations;
	extern Molecule molecule;

	// M. C. Zerner, M. Hehenberger, A dynamical damping scheme for converging molecular scf calculations, Chem. Phys. Lett. 62 (3), 550-554 (1979).
	// The method has been somewhat modified because I am lazy.
	double getZernerHehenbergerDamping() {
		populations.resize(molecule.size());

		calculateMullikenPopulations();

		uint atom_index = 0;
		uint weight = 0;
		double damping = 0.0;
		double atom_population = 0.0;
		for (int i = 0;; ++i) {
			bool append = i == matrix_size;
			if (!append) append |= atom_index != basis_atoms[i];
			if (append) {
				double alpha = atom_population / (atom_population - populations[atom_index]);
				if (alpha < 0.0 || alpha >= 1.0 || isnan(alpha)) alpha = 0.0;

				populations[atom_index] = alpha * populations[atom_index] + (1.0 - alpha) * atom_population;

				damping += alpha * (double)weight;

				weight = 0;
				atom_population = 0.0;

				if (i == matrix_size) break;

				atom_index = basis_atoms[i];
			}

			double out_population = mulliken_populations[0].at(i, i);
			if (spin_treatment == SpinTreatment::unrestricted) out_population += mulliken_populations[1].at(i, i);

			atom_population += out_population;

			++weight;
		}

		damping /= (double)matrix_size;

		return damping;
	}
}