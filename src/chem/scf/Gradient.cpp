#include "Gradient.hpp"

#include "SCFCommon.hpp"
#include "ExactCoulomb.hpp"

#include "../Molecule.hpp"
#include "../GaussianIntegral.hpp"

using namespace flo;

namespace scf {
	extern Molecule molecule;

	void getEnergyGradient(flo::VectorNd& result) {
		result.resize(molecule.size() * 3);

		result *= 0.0;

		SymmetricMatrixNd& weighted_density_matrix = buffer[0];

		if (spin_treatment == SpinTreatment::unrestricted) {
			for (int mu = 0; mu < basis.size(); ++mu) {
				for (int nu = mu; nu < basis.size(); ++nu) {
					weighted_density_matrix.at(mu, nu) = coefficient_matrix[0].at(mu, 0) * coefficient_matrix[0].at(nu, 0) * mo_levels[0][0];
					if (electron_count[1]) weighted_density_matrix.at(mu, nu) += coefficient_matrix[1].at(mu, 0) * coefficient_matrix[1].at(nu, 0) * mo_levels[1][0];
				}
			}
			for (int i = 1; i < electron_count[0]; ++i) {
				for (int mu = 0; mu < basis.size(); ++mu) {
					for (int nu = mu; nu < basis.size(); ++nu) {
						weighted_density_matrix.at(mu, nu) += coefficient_matrix[0].at(mu, i) * coefficient_matrix[0].at(nu, i) * mo_levels[0][i];
					}
				}
			}
			for (int i = 1; i < electron_count[1]; ++i) {
				for (int mu = 0; mu < basis.size(); ++mu) {
					for (int nu = mu; nu < basis.size(); ++nu) {
						weighted_density_matrix.at(mu, nu) += coefficient_matrix[1].at(mu, i) * coefficient_matrix[1].at(nu, i) * mo_levels[1][i];
					}
				}
			}
		}
		else {
			for (int mu = 0; mu < basis.size(); ++mu) {
				for (int nu = mu; nu < basis.size(); ++nu) {
					weighted_density_matrix.at(mu, nu) = 2.0 * coefficient_matrix[0].at(mu, 0) * coefficient_matrix[0].at(nu, 0) * mo_levels[0][0];
				}
			}
			for (int i = 1; i < electron_count[0]; ++i) {
				for (int mu = 0; mu < basis.size(); ++mu) {
					for (int nu = mu; nu < basis.size(); ++nu) {
						weighted_density_matrix.at(mu, nu) += 2.0 * coefficient_matrix[0].at(mu, i) * coefficient_matrix[0].at(nu, i) * mo_levels[0][i];
					}
				}
			}
		}

		for (int i = 0; i < basis.size(); ++i) {
			uint atom_i = basis_atoms[i];
			ContractedGaussian& chi_mu = basis[i];
			for (int j = i; j < basis.size(); ++j) {
				uint atom_j = basis_atoms[j];
				ContractedGaussian& chi_nu = basis[j];

				double density = total_density_matrix.at(i, j);
				if (i != j) density *= 2.0;

				if (atom_i != atom_j) {
					vec3 kinetic_gradient = density * kineticEnergyGradient(chi_mu, chi_nu);

					result[atom_i * 3]     += kinetic_gradient.x;
					result[atom_i * 3 + 1] += kinetic_gradient.y;
					result[atom_i * 3 + 2] += kinetic_gradient.z;

					result[atom_j * 3]     -= kinetic_gradient.x;
					result[atom_j * 3 + 1] -= kinetic_gradient.y;
					result[atom_j * 3 + 2] -= kinetic_gradient.z;

					double weighted_density = weighted_density_matrix.at(i, j);
					if (i != j) weighted_density *= 2.0;

					vec3 overlap_gradient = weighted_density * overlapGradient(chi_mu, chi_nu);

					result[atom_i * 3]     -= overlap_gradient.x;
					result[atom_i * 3 + 1] -= overlap_gradient.y;
					result[atom_i * 3 + 2] -= overlap_gradient.z;

					result[atom_j * 3]     += overlap_gradient.x;
					result[atom_j * 3 + 1] += overlap_gradient.y;
					result[atom_j * 3 + 2] += overlap_gradient.z;
				}

				for (int c = 0; c < molecule.size(); ++c) {
					if (atom_i == atom_j && c == atom_i) continue;

					auto nuclear_gradient = nuclearPotentialGradient(chi_mu, chi_nu, molecule[c].position);

					vec3 gradient_mu = -density * molecule[c].charge * nuclear_gradient[0];
					vec3 gradient_nu = -density * molecule[c].charge * nuclear_gradient[1];

					result[c * 3]     -= gradient_mu.x + gradient_nu.x;
					result[c * 3 + 1] -= gradient_mu.y + gradient_nu.y;
					result[c * 3 + 2] -= gradient_mu.z + gradient_nu.z;

					result[atom_i * 3]     += gradient_mu.x;
					result[atom_i * 3 + 1] += gradient_mu.y;
					result[atom_i * 3 + 2] += gradient_mu.z;

					result[atom_j * 3]     += gradient_nu.x;
					result[atom_j * 3 + 1] += gradient_nu.y;
					result[atom_j * 3 + 2] += gradient_nu.z;
				}
			}
		}

		addExactTwoElectronGradient(result);

		for (int i = 0; i < molecule.size(); ++i) {
			for (int j = i + 1; j < molecule.size(); ++j) {
				if (i == j) continue;

				vec3 gradient = (int)molecule[j].charge * molecule[i].getRepulsionGradient(molecule[j].position);

				result[i * 3]     -= gradient.x;
				result[i * 3 + 1] -= gradient.y;
				result[i * 3 + 2] -= gradient.z;

				result[j * 3]     += gradient.x;
				result[j * 3 + 1] += gradient.y;
				result[j * 3 + 2] += gradient.z;
			}
		}
	}
}