#include "Energy.hpp"

#include "../../lalib/Lalib.hpp"

#include "../Molecule.hpp"

#include "SCFCommon.hpp"

using namespace flo;

namespace scf {
	extern Molecule molecule;

	double kinetic_energy = 0.0;
	double nuclear_attraction_energy = 0.0;
	double electron_repulsion_energy = 0.0;
	double exchange_energy = 0.0;
	double correlation_energy = 0.0;
	double nuclear_energy = 0.0;

	void calculateEnergies() {
		kinetic_energy = trace(total_density_matrix * kinetic_energy_matrix);

		nuclear_attraction_energy = trace(total_density_matrix * nuclear_attraction_matrix);

		computeCoulombTerms(Spin::alpha);
		if (spin_treatment == SpinTreatment::unrestricted) computeCoulombTerms(Spin::beta);
		electron_repulsion_energy = 0.5 * trace(density_matrix[0] * two_electron_hamiltonian[0]);
		if (spin_treatment == SpinTreatment::unrestricted) electron_repulsion_energy += 0.5 * trace(density_matrix[1] * two_electron_hamiltonian[1]);
		else electron_repulsion_energy *= 2.0;

		two_electron_hamiltonian[0] = 0.0;
		computeExchangeTerms(Spin::alpha);
		if (spin_treatment == SpinTreatment::unrestricted) {
			two_electron_hamiltonian[1] = 0.0;
			computeExchangeTerms(Spin::beta);
		}
		exchange_energy = 0.5 * trace(density_matrix[0] * two_electron_hamiltonian[0]);
		if (spin_treatment == SpinTreatment::unrestricted) exchange_energy += 0.5 * trace(density_matrix[1] * two_electron_hamiltonian[1]);
		else exchange_energy *= 2.0;

		correlation_energy = 0.0;

		nuclear_energy = 0.0;
		for (int i = 0; i < molecule.size(); ++i) {
			for (int j = i + 1; j < molecule.size(); ++j) {
				nuclear_energy += (int)molecule[i].element * (int)molecule[j].element / length(molecule[i].position - molecule[j].position);
			}
		}
	}

	double getKineticEnergy() {
		return kinetic_energy;
	}

	double getNuclearAttractionEnergy() {
		return nuclear_attraction_energy;
	}

	double getElectronRepulsionEnergy() {
		return electron_repulsion_energy;
	}

	double getExchangeEnergy() {
		return exchange_energy;
	}

	double getCorrelationEnergy() {
		return correlation_energy;
	}

	double getNuclearRepulsionEnergy() {
		return nuclear_energy;
	}

	double getTotalEnergy() {
		return kinetic_energy + nuclear_attraction_energy + electron_repulsion_energy + exchange_energy + correlation_energy + nuclear_energy;
	}
}