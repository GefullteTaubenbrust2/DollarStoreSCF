#include "SCFCommon.hpp"

namespace scf {
	bool open_shell = false;
	int net_charge = 0;

	std::vector<flo::ContractedGaussian> basis;
	std::vector<uint> basis_atoms;
	std::vector<uint> atom_basis;

	flo::SymmetricMatrixNd kinetic_energy_matrix;
	flo::SymmetricMatrixNd nuclear_attraction_matrix;
	flo::SymmetricMatrixNd core_hamiltonian;
	flo::SymmetricMatrixNd overlap_matrix;
	flo::SymmetricMatrixNd orthogonalization_matrix;
	flo::SymmetricMatrixNd lowdin_matrix;
	std::vector<flo::SymmetricMatrixNd> buffer;
	std::vector<flo::MatrixNd> asymmetric_buffer;

	flo::VectorNd mo_levels[2];
	flo::MatrixNd coefficient_matrix[2];
	flo::SymmetricMatrixNd fock_matrix[2];
	flo::SymmetricMatrixNd density_matrix[2];
	flo::SymmetricMatrixNd total_density_matrix;

	uint electron_count[2];

	uint matrix_size;

	double nuclear_energy;

	SpinTreatment spin_treatment;
}