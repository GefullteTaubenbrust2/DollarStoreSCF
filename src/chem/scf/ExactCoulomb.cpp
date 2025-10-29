#include "ExactCoulomb.hpp"
#include "SCFCommon.hpp"

#include "../GaussianIntegral.hpp"

#include "../../util/FormattedStream.hpp"
#include "../../util/Time.hpp"

using namespace flo;

#define INTEGRAL_THRESHOLD 0.00000001

namespace scf {
	std::vector<double> repulsion_tensor;
	size_t diagonal_size = 0;
	size_t tensor_size = 0;

	void printExactRepulsionStart(const std::string& job_title) {
		fout.resetRows();
		fout.addRow(NumberFormat(), flo::TextAlignment::centered, 50);
		fout.centerTable(100);
		fout << '|' << '-' << " " << '|' << '\n';
		fout << '|' << job_title << '|' << '\n';
		fout << '|' << '_' << " " << '|' << '\n';
		fout.resetRows();

		global_stopclock.reset();
	}

	void printExactRepulsionEnd() {
		double runtime = global_stopclock.stop();

		fout.resetRows();
		fout.addRow(NumberFormat(), flo::TextAlignment::centered, 33);
		fout.centerTable(100);
		fout << '|' << '_' << '-' << "Done in" << runtime << " s!" << '|' << '\n';
		fout.resetRows();
		fout << '\n';
	}

	void assignExactRepulsionTensor() {
		printExactRepulsionStart("Computing exact four-center integrals... ");

		SymmetricMatrixNd diagonal_terms(basis.size());
		for (uint a = 0; a < basis.size(); ++a) {
			for (uint b = 0; b <= a; ++b) {
				diagonal_terms.at(b, a) = electronRepulsionIntegral(basis[a], basis[b], basis[a], basis[b]);
			}
		}
		diagonal_size = basis.size() * (basis.size() + 1) / 2;
		tensor_size = diagonal_size * (diagonal_size + 1) / 2;
		repulsion_tensor.resize(tensor_size);

		for (uint i = 0; i < basis.size(); ++i) {
			for (uint j = 0; j <= i; ++j) {
				for (uint k = 0; k < basis.size(); ++k) {
					for (uint l = 0; l <= k; ++l) {
						uint m = j + i * (i + 1) / 2;
						uint n = l + k * (k + 1) / 2;
						if (m > n) continue;
						size_t index = tensor_size - (diagonal_size - m) * (diagonal_size - m + 1) / 2 + n - m;
						
						if (diagonal_terms(i, j) * diagonal_terms(k, l) > INTEGRAL_THRESHOLD * INTEGRAL_THRESHOLD) {
							repulsion_tensor[index] = flo::electronRepulsionIntegral(basis[i], basis[k], basis[j], basis[l]);
						}
					}
				}
			}
		}

		printExactRepulsionEnd();
	}

	double sampleRepulsionTensor(uint i, uint k, uint j, uint l) {
		if (j > i) {
			uint x = i;
			i = j;
			j = x;
		}
		if (l > k) {
			uint x = k;
			k = l;
			l = x;
		}
		uint m = j + i * (i + 1) / 2;
		uint n = l + k * (k + 1) / 2;
		if (m > n) {
			uint x = i;
			i = k;
			k = x;

			x = j;
			j = l;
			l = x;

			m = j + i * (i + 1) / 2;
			n = l + k * (k + 1) / 2;
		}
		size_t index = tensor_size - (diagonal_size - m) * (diagonal_size - m + 1) / 2 + n - m;
		return repulsion_tensor[index];
	}

	double getExactCoulombMatrix(uint mu, uint nu) {
		double sum = 0.0;
		for (uint tau = 0; tau < basis.size(); ++tau) {
			for (uint sigma = 0; sigma <= tau; ++sigma) {
				if (tau == sigma) sum += total_density_matrix(sigma, tau) * sampleRepulsionTensor(mu, sigma, nu, tau);
				else sum += 2.0 * total_density_matrix(sigma, tau) * sampleRepulsionTensor(mu, sigma, nu, tau);
			}
		}
		return sum;
	}

	double getExactExchangeMatrix(uint mu, uint nu, Spin spin) {
		double sum = 0.0;
		SymmetricMatrixNd& density = density_matrix[(int)spin];
		for (uint sigma = 0; sigma < basis.size(); ++sigma) {
			for (uint tau = 0; tau < basis.size(); ++tau) {
				sum += density(sigma, tau) * sampleRepulsionTensor(mu, sigma, tau, nu);
			}
		}
		return sum;
	}

	void addExactTwoElectronGradient(VectorNd& result) {
		printExactRepulsionStart("Computing exact four-center gradients... ");

		for (uint i = 0; i < basis.size(); ++i) {
			uint atom_i = basis_atoms[i];
			for (uint j = 0; j <= i; ++j) {
				uint atom_j = basis_atoms[j];
				for (uint k = 0; k < basis.size(); ++k) {
					uint atom_k = basis_atoms[k];
					for (uint l = 0; l <= k; ++l) {
						uint atom_l = basis_atoms[l];

						if (atom_i == atom_j && atom_k == atom_l && atom_i == atom_k) continue;

						uint m = j + i * (i + 1) / 2;
						uint n = l + k * (k + 1) / 2;
						if (m > n) continue;

						auto integrals = electronRepulsionGradient(basis[i], basis[k], basis[j], basis[l]);

						double coulomb_density = total_density_matrix.at(j, i) * total_density_matrix.at(l, k);

						double exchange_density = density_matrix[0](j, k) * density_matrix[0](i, l);
						if (spin_treatment == SpinTreatment::unrestricted) exchange_density += density_matrix[1](j, k) * density_matrix[1](i, l);
						else exchange_density *= 2.0;

						if ((i != j) || (k != l)) {
							if (spin_treatment == SpinTreatment::unrestricted) exchange_density +=
								density_matrix[0](i, k) * density_matrix[0](j, l) +
								density_matrix[1](i, k) * density_matrix[1](j, l);
							else exchange_density += 2.0 * density_matrix[0](i, k) * density_matrix[0](j, l);
							exchange_density *= 0.5;
						}

						if (i != j) {
							coulomb_density *= 2.0;
							exchange_density *= 2.0;
						}
						if (k != l) {
							coulomb_density *= 2.0;
							exchange_density *= 2.0;
						}
						if (m != n) {
							coulomb_density *= 2.0;
							exchange_density *= 2.0;
						}

						double factor = 0.5 * (coulomb_density - exchange_density);

						result[atom_i * 3] += factor * integrals[0].x;
						result[atom_i * 3 + 1] += factor * integrals[0].y;
						result[atom_i * 3 + 2] += factor * integrals[0].z;

						result[atom_k * 3] += factor * integrals[1].x;
						result[atom_k * 3 + 1] += factor * integrals[1].y;
						result[atom_k * 3 + 2] += factor * integrals[1].z;

						result[atom_j * 3] += factor * integrals[2].x;
						result[atom_j * 3 + 1] += factor * integrals[2].y;
						result[atom_j * 3 + 2] += factor * integrals[2].z;

						result[atom_l * 3] += factor * integrals[3].x;
						result[atom_l * 3 + 1] += factor * integrals[3].y;
						result[atom_l * 3 + 2] += factor * integrals[3].z;
					}
				}
			}
		}

		printExactRepulsionEnd();
	}
}