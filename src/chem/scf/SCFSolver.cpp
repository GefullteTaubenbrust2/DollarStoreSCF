#include "SCFSolver.hpp"

#include "SCFCommon.hpp"
#include "ExactCoulomb.hpp"
#include "Damping.hpp"
#include "SCFDIIS.hpp"
#include "Energy.hpp"

#include "../GTO.hpp"

#include "../../util/FormattedStream.hpp"

#include "../../lalib/Lapack.hpp"
#include "../../lalib/PrintMatrix.h"

using namespace flo;

namespace scf {
	uint iteration = 0;
	
	double previous_electronic_energy = 0.0;
	double total_electronic_energy = 0.0;

	Array<double> diagonalization_buffer;

	double energy_threshold = 0.00000001;
	double rms = 0.0;

	void computeCoulombTerms(Spin spin) {
		for (uint mu = 0; mu < matrix_size; ++mu) {
			for (uint nu = mu; nu < matrix_size; ++nu) {
				two_electron_hamiltonian[(int)spin].at(mu, nu) = getExactCoulombMatrix(mu, nu);
			}
		}
	}

	void computeExchangeTerms(Spin spin) {
		for (uint mu = 0; mu < matrix_size; ++mu) {
			for (uint nu = mu; nu < matrix_size; ++nu) {
				two_electron_hamiltonian[(int)spin].at(mu, nu) -= getExactExchangeMatrix(mu, nu, spin);
			}
		}
	}

	void updateFockMatrix(Spin spin) {
		/*double damping = getZernerHehenbergerDamping();
		//double damping = 0.0;
		if (iteration == 1 || damping == 0.0) fock_matrix[(int)spin] = core_hamiltonian + two_electron_hamiltonian[(int)spin];
		else {
			fock_matrix[(int)spin] = (buffer[1] = (1.0 - damping) * (buffer[0] = core_hamiltonian + two_electron_hamiltonian[(int)spin])) + (buffer[2] = damping * fock_matrix[(int)spin]);
		}*/
		buffer[0] = core_hamiltonian + two_electron_hamiltonian[(int)spin];

		if (iteration == 1) fock_matrix[(int)spin] = buffer[0];
		else {
			diis::addPulayErrorVector(buffer[0], spin);
			diis::solve();
			fock_matrix[(int)spin] = diis::getResult(spin);
		}

		rms = diis::getRMS();
	}

	void solveRoothaanHall(Spin spin) {
		SymmetricMatrixNd& orthogonal_fock_matrix = buffer[0];
		orthogonal_fock_matrix = (asymmetric_buffer[0] = orthogonalization_matrix * fock_matrix[(int)spin]) * orthogonalization_matrix;

		diagonalization_buffer.resize(matrix_size * 3);

		buffer[1] = buffer[0];

		computeEigenvectors(orthogonal_fock_matrix, asymmetric_buffer[0], mo_levels[(int)spin], diagonalization_buffer.getPtr());

		coefficient_matrix[(int)spin] = orthogonalization_matrix * asymmetric_buffer[0];
	}

	void computeDensityMatrix(const Spin spin) {
		if (electron_count[(int)spin]) {
			for (uint mu = 0; mu < matrix_size; ++mu) {
				for (uint nu = mu; nu < matrix_size; ++nu) {
					density_matrix[(int)spin].at(mu, nu) = coefficient_matrix[(int)spin].at(mu, 0) * coefficient_matrix[(int)spin].at(nu, 0);
				}
			}
		}
		for (uint i = 1; i < electron_count[(int)spin]; ++i) {
			for (uint mu = 0; mu < matrix_size; ++mu) {
				for (uint nu = mu; nu < matrix_size; ++nu) {
					density_matrix[(int)spin].at(mu, nu) += coefficient_matrix[(int)spin].at(mu, i) * coefficient_matrix[(int)spin].at(nu, i);
				}
			}
		}
		if (spin_treatment == SpinTreatment::unrestricted) total_density_matrix = density_matrix[0] + density_matrix[1];
		else total_density_matrix = 2.0 * density_matrix[0];
	}

	void printSpinContamination() {
		double spin_overlap = trace(overlap_matrix * (asymmetric_buffer[0] = density_matrix[0] * (asymmetric_buffer[1] = overlap_matrix * (asymmetric_buffer[2] = density_matrix[1]))));
		uint spin = electron_count[0] - electron_count[1];
		double spin_squared = (double)(spin * spin + 2 * spin) * 0.25;
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 55);
		fout << '|' << '-' << '_' << "UHF/UKS spin contamination\n\nFor DFT calculations, this has less meaning than for Hartree-Fock." << '|';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 15);
		fout.addRow(NumberFormat(), TextAlignment::centered, 15);
		fout.addRow(NumberFormat(), TextAlignment::centered, 15);
		fout << '|' << '_' << "S(S+1)" << '|' << '_' << "<S^2>" << '|' << '_' << "Spin\ncontamination" << '|' << '\n';
		fout.resetRows();
		fout.addRow(NumberFormat::crudeFormatPositive(15, 5), TextAlignment::right, 15);
		fout.addRow(NumberFormat::crudeFormatPositive(15, 5), TextAlignment::right, 15);
		fout.addRow(NumberFormat::crudeFormatPositive(15, 5), TextAlignment::right, 15);
		fout << '|' << '_' << (double)spin_squared << '|' << '_' << (spin_squared + (double)electron_count[1] - spin_overlap) << '|' << '_' << ((double)electron_count[1] - spin_overlap) << '|' << '\n';

		fout.resetRows();
		fout << '\n';
	}

	double getTotalElectronicEnergy() {
		double energy = 0.0;
		if (spin_treatment == SpinTreatment::restricted) {
			for (int i = 0; i < electron_count[0]; ++i) {
				energy += 2.0 * mo_levels[0][i];
			}

			energy -= trace(density_matrix[0] * two_electron_hamiltonian[0]);

			return energy;
		}
		else {
			for (int i = 0; i < electron_count[0]; ++i) {
				energy += mo_levels[0][i];
			}
			for (int i = 0; i < electron_count[1]; ++i) {
				energy += mo_levels[1][i];
			}

			energy -= 0.5 * trace(density_matrix[0] * two_electron_hamiltonian[0]);

			energy -= 0.5 * trace(density_matrix[1] * two_electron_hamiltonian[1]);

			return energy;
		}
	}

	void runSCFIteration() {
		computeCoulombTerms(Spin::alpha);

		computeExchangeTerms(Spin::alpha);

		updateFockMatrix(Spin::alpha);

		if (spin_treatment == SpinTreatment::unrestricted) {
			computeCoulombTerms(Spin::beta);

			computeExchangeTerms(Spin::beta);

			updateFockMatrix(Spin::beta);
		}

		solveRoothaanHall(Spin::alpha);

		computeDensityMatrix(Spin::alpha);

		if (spin_treatment == SpinTreatment::unrestricted) {
			solveRoothaanHall(Spin::beta);

			computeDensityMatrix(Spin::beta);
		}

		previous_electronic_energy = total_electronic_energy;
		total_electronic_energy = getTotalElectronicEnergy();
	}

	bool checkConvergence() {
		if (!iteration) return false;
		return std::abs(total_electronic_energy - previous_electronic_energy) < energy_threshold;
	}

	void useCoreGuess() {
		for (uint mu = 0; mu < basis.size(); ++mu) {
			for (uint nu = mu; nu < basis.size(); ++nu) {
				density_matrix[0].at(mu, nu) = 0.0;
				total_density_matrix.at(mu, nu) = 0.0;
			}
		}

		if (spin_treatment == SpinTreatment::unrestricted) {
			for (uint mu = 0; mu < basis.size(); ++mu) {
				for (uint nu = mu; nu < basis.size(); ++nu) {
					density_matrix[1].at(mu, nu) = 0.0;
				}
			}
		}
	}

	void solveMOs() {
		iteration = 0;

		fout.resetRows();
		fout.addRow(NumberFormat(), flo::TextAlignment::left, 9);
		fout.addRow(NumberFormat(), flo::TextAlignment::centered, 16);
		fout.addRow(NumberFormat(), flo::TextAlignment::centered, 16);
		fout.addRow(NumberFormat(), flo::TextAlignment::centered, 16);
		fout.offsetRight(2);

		fout << '-' << '_' << '|' << "SCF\nIteration" << '|' << '-' << '_' << "Total elec-\ntronic energy" << '|' << '-' << '_' << "Delta E" << '|' << '-' << '_' << "[PF]^2" << '|';

		fout.resetRows();
		fout.addRow(NumberFormat::crudeFormatPositive(6, 6), flo::TextAlignment::left, 9);
		fout.addRow(NumberFormat::crudeFormat(16, 12), flo::TextAlignment::right, 16);
		fout.addRow(NumberFormat::crudeFormat(16, 12), flo::TextAlignment::right, 16);
		fout.addRow(NumberFormat::crudeFormat(16, 12), flo::TextAlignment::right, 16);

		diis::setIterationCount(8);

		total_electronic_energy = 0.0;
		previous_electronic_energy = 0.0;

		while (!checkConvergence()) {
			++iteration;
			runSCFIteration();
			fout << '|' << (i64)iteration << '|' << fout.formatRow(TextAlignment::right) << total_electronic_energy;
			fout << '|' << (total_electronic_energy - previous_electronic_energy) << '|' << rms << '|' << '\n';
		}
		fout << '-' << ',' << '-' << ',' << '-' << ',' << '-' << '\n';

		fout.resetRows();
		fout.addRow(NumberFormat(), flo::TextAlignment::centered, 72);
		fout << '|' << '_' << "Single point calculation has converged \\(^o^)/" << '|' << '\n';

		fout.resetRows();
		fout << '\n';

		if (spin_treatment == SpinTreatment::unrestricted) printSpinContamination();

		calculateEnergies();
	}
}