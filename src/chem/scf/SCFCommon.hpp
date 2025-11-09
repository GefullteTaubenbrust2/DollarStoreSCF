#pragma once
#include "../../lalib/Lalib.hpp"
#include "../GTO.hpp"
#include "SCFSolver.hpp"

namespace scf {
	extern bool open_shell;
	extern int net_charge;

	extern std::vector<flo::ContractedGaussian> basis;
	extern std::vector<uint> basis_atoms;
	extern std::vector<uint> atom_basis;

	extern flo::SymmetricMatrixNd kinetic_energy_matrix;
	extern flo::SymmetricMatrixNd nuclear_attraction_matrix;
	extern flo::SymmetricMatrixNd core_hamiltonian;
	extern flo::SymmetricMatrixNd two_electron_hamiltonian[2];
	extern flo::SymmetricMatrixNd fock_matrix[2];

	extern flo::SymmetricMatrixNd overlap_matrix;
	extern flo::SymmetricMatrixNd orthogonalization_matrix;
	extern flo::SymmetricMatrixNd lowdin_matrix;

	extern std::vector<flo::SymmetricMatrixNd> buffer;
	extern std::vector<flo::MatrixNd> asymmetric_buffer;
	extern std::vector<flo::VectorNd> vector_buffer;

	extern flo::VectorNd mo_levels[2];
	extern flo::MatrixNd coefficient_matrix[2];
	extern flo::SymmetricMatrixNd density_matrix[2];
	extern flo::SymmetricMatrixNd total_density_matrix;

	extern uint electron_count[2];

	extern uint matrix_size;

	extern SpinTreatment spin_treatment;
}