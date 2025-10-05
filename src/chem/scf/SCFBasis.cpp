#include "SCFBasis.hpp"
#include <iostream>
#include "../GaussianIntegral.hpp"
#include "../../lalib/LalibDeclarations.hpp"
#include "../../lalib/Lapack.hpp"
#include "SCFCommon.hpp"

#define MATRIX_BUFFER_SIZE 5

using namespace flo;

namespace scf {
	Molecule molecule;
	std::vector<AtomicBasis> atom_bases;

	void assignMolecule(const Molecule& _molecule, int charge, uint multiplicity) {
		molecule = _molecule;
		atom_bases.resize(molecule.size());
		uint total_electron_count = 0;
		nuclear_energy = 0.0;
		for (int i = 0; i < molecule.size(); ++i) {
			atom_bases[i].clear();
			total_electron_count += (uint)molecule[i].element;
			for (int j = i + 1; j < molecule.size(); ++j) {
				nuclear_energy += (int)molecule[i].element * (int)molecule[j].element / length(molecule[i].position - molecule[j].position);
			}
		}
		total_electron_count -= charge;
		if (total_electron_count % 2 == multiplicity % 2) {
			std::cerr << "ERROR: The specified multiplicity does not match the electron count. Molecules with odd numbers of electrons and have even multiplicity and vice versa!\n";
		}
		electron_count[0] = (total_electron_count + multiplicity) / 2;
		electron_count[1] = (total_electron_count + 1 - multiplicity) / 2;

		if (multiplicity == 1) spin_treatment = SpinTreatment::restricted;
		else spin_treatment = SpinTreatment::unrestricted;

		open_shell = multiplicity > 1;
	}

	void assignBasis(const BasisSet& basis_set, uint atom_index) {
		if (atom_index >= molecule.size()) return;
		atom_bases[atom_index] = basis_set.atomic_bases[(int)molecule[atom_index].element - 1];
	}

	void assignBasis(const BasisSet& basis_set) {
		for (int i = 0; i < molecule.size(); ++i) {
			assignBasis(basis_set, i);
		}
	}

	void assignOverlapMatrix() {
		overlap_matrix.resize(basis.size());
		for (int i = 0; i < basis.size(); ++i) {
			for (int j = 0; j <= i; ++j) {
				overlap_matrix.at(j, i) = overlapIntegral(basis[j], basis[i]);
			}
		}
	}

	void assignCoreHamiltonian() {
		core_hamiltonian.resize(basis.size());
		kinetic_energy_matrix.resize(basis.size());
		nuclear_attraction_matrix.resize(basis.size());
		for (int i = 0; i < basis.size(); ++i) {
			for (int j = 0; j <= i; ++j) {
				kinetic_energy_matrix.at(j, i) = kineticEnergyIntegral(basis[j], basis[i]);
			}
		}
		for (int i = 0; i < basis.size(); ++i) {
			for (int j = 0; j <= i; ++j) {
				nuclear_attraction_matrix.at(j, i) = 0.0;
				for (Atom& atom : molecule.atoms) {
					nuclear_attraction_matrix.at(j, i) -= (double)atom.charge * nuclearPotentialIntegral(basis[j], basis[i], atom.position);
				}
			}
		}
		core_hamiltonian = kinetic_energy_matrix + nuclear_attraction_matrix;
	}

	void assignOrthogonalization() {
		buffer[0] = overlap_matrix;
		VectorNd eigenvalues(basis.size());
		computeEigenvectors(buffer[0], asymmetric_buffer[0], eigenvalues);
		DiagonalMatrixNd eigenvalue_matrix(basis.size());

		for (int i = 0; i < basis.size(); ++i) {
			eigenvalue_matrix.at(i, i) = std::sqrt(eigenvalues[i]);
		}
		lowdin_matrix = (asymmetric_buffer[2] = asymmetric_buffer[0] * eigenvalue_matrix) * (asymmetric_buffer[1] = T(asymmetric_buffer[0]));

		for (int i = 0; i < basis.size(); ++i) {
			eigenvalue_matrix.at(i, i) = 1.0 / eigenvalue_matrix.at(i, i);
		}
		orthogonalization_matrix = (asymmetric_buffer[2] = asymmetric_buffer[0] * eigenvalue_matrix) * (asymmetric_buffer[1] = T(asymmetric_buffer[0]));
	}

	int p_permutations[3] = { 0, -1, 1 };
	int d_permutations[5] = { 0, -1, 1, -2, 2 };

	void constructBasis() {
		basis.clear();
		for (int i = 0; i < molecule.size(); ++i) {
			AtomicBasis& definition = atom_bases[i];
			if (!definition.size()) {
				std::cerr << "ERROR: No basis assigned to " << element_names[(int)molecule[i].element - 1] << " atom number " << (i + 1) << "!\n";
				return;
			}
		}

		atom_basis.reserve(molecule.size() + 1);
		atom_basis.push_back(0);

		for (int i = 0; i < molecule.size(); ++i) {
			vec3 center = molecule[i].position;
			AtomicBasis& definition = atom_bases[i];

			for (GTOShell& shell : definition) {
				uint l = shell.l;
				basis.reserve(basis.size() + 2 * l + 1);
				basis_atoms.reserve(basis.size() + 2 * l + 1);

				for (int m = -(int)l; m <= (int)l; ++m) {
					int ma = m;
					if (l == 1) ma = p_permutations[m + 1];
					if (l == 2) ma = d_permutations[m + 2];
					basis.push_back(ContractedGaussian(shell.primitives, l, ma, center));
					basis_atoms.push_back(i);
				}
			}

			atom_basis.push_back(basis.size());
		}

		matrix_size = basis.size();

		buffer.resize(MATRIX_BUFFER_SIZE);
		asymmetric_buffer.resize(MATRIX_BUFFER_SIZE);
		for (int i = 0; i < MATRIX_BUFFER_SIZE; ++i) {
			buffer[i].resize(basis.size());
			asymmetric_buffer[i].resize(basis.size(), basis.size());
		}

		assignOverlapMatrix();

		assignCoreHamiltonian();

		assignOrthogonalization();
	}
}