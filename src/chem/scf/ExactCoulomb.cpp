#include "ExactCoulomb.hpp"
#include "SCFCommon.hpp"
#include "../GaussianIntegral.hpp"

using namespace flo;

#define INTEGRAL_THRESHOLD 0.00000001

namespace scf {
	std::vector<double> repulsion_tensor;
	size_t diagonal_size = 0;
	size_t tensor_size = 0;

	void assignExactRepulsionTensor() {
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
		if (spin_treatment == SpinTreatment::restricted) {
			for (uint tau = 0; tau < basis.size(); ++tau) {
				for (uint sigma = 0; sigma <= tau; ++sigma) {
					if (tau == sigma) sum += density_matrix[0](sigma, tau) * sampleRepulsionTensor(mu, sigma, nu, tau);
					else sum += 2.0 * density_matrix[0](sigma, tau) * sampleRepulsionTensor(mu, sigma, nu, tau);
				}
			}
			sum *= 2.0;
		}
		else {
			for (uint tau = 0; tau < basis.size(); ++tau) {
				for (uint sigma = 0; sigma <= tau; ++sigma) {
					if (tau == sigma) sum += (density_matrix[0](sigma, tau) + density_matrix[1](sigma, tau)) * sampleRepulsionTensor(mu, sigma, nu, tau);
					else sum += 2.0 * (density_matrix[0](sigma, tau) + density_matrix[1](sigma, tau)) * sampleRepulsionTensor(mu, sigma, nu, tau);
				}
			}
		}
		return sum;
	}

	double getExactExchangeMatrix(uint mu, uint nu, const Spin spin) {
		double sum = 0.0;
		SymmetricMatrixNd& density = density_matrix[(int)spin];
		for (uint sigma = 0; sigma < basis.size(); ++sigma) {
			for (uint tau = 0; tau < basis.size(); ++tau) {
				sum += density(sigma, tau) * sampleRepulsionTensor(mu, sigma, tau, nu);
			}
		}
		return sum;
	}
}