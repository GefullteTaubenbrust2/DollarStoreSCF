#include "DIIS.hpp"
#include "SCFCommon.hpp"
#include "../../lalib/Lapack.hpp"

using namespace flo;

namespace scf {
namespace diis {
	uint diis_iterations = 0;
	std::vector<MatrixNd> error_vectors[2];
	std::vector<SymmetricMatrixNd> fock_matrices[2];
	SymmetricMatrixNd b_matrix[2];
	VectorNd coefficients[2];

	bool has_rms = false;
	double rms = 0.0;

	void setIterationCount(uint iterations) {
		diis_iterations = iterations;

		error_vectors[0].reserve(iterations);
		error_vectors[0].resize(0);
		fock_matrices[0].reserve(iterations);
		fock_matrices[0].resize(0);
		b_matrix[0].resize(iterations + 1);
		b_matrix[0].resize(0);
		coefficients[0].resize(iterations + 1);
		coefficients[0].resize(0);

		if (spin_treatment == SpinTreatment::unrestricted) {
			error_vectors[1].reserve(iterations);
			error_vectors[1].resize(0);
			fock_matrices[1].reserve(iterations);
			fock_matrices[1].resize(0);
			b_matrix[1].resize(iterations + 1);
			b_matrix[1].resize(0);
			coefficients[1].resize(iterations + 1);
			coefficients[1].resize(0);
		}
	}

	void computeErrorVector(flo::MatrixNd& error_vector, const flo::SymmetricMatrixNd& fock_matrix, Spin spin) {
		MatrixNd& spf = asymmetric_buffer[1];

		spf = (asymmetric_buffer[2] = overlap_matrix * (asymmetric_buffer[3] = density_matrix[(int)spin] * fock_matrix));

		error_vector = spf - (asymmetric_buffer[2] = T(spf));
	}

	void addPulayErrorVector(const flo::SymmetricMatrixNd& fock_matrix, Spin spin) {
		MatrixNd& error_vector = asymmetric_buffer[0];
		
		computeErrorVector(error_vector, fock_matrix, spin);

		double new_rms = dot(error_vector, error_vector);

		rms = new_rms;
		has_rms = true;

		int index = 0;

		int vector_count = error_vectors[(int)spin].size();
		if (error_vectors[(int)spin].size() < diis_iterations) {
			error_vectors[(int)spin].resize(vector_count + 1);
			fock_matrices[(int)spin].resize(vector_count + 1);
			b_matrix[(int)spin].resize(vector_count + 2);
			coefficients[(int)spin].resize(vector_count + 2);

			index = vector_count;
			++vector_count;

			for (int i = 0; i < vector_count; ++i) {
				b_matrix[(int)spin].at(i, vector_count) = 1.0;
			}
			b_matrix[(int)spin].at(vector_count, vector_count) = 0.0;
		}
		else {
			double max_rms = 0.0;
			for (int i = 0; i < diis_iterations; ++i) {
				double error = b_matrix[(int)spin].at(i, i);
				if (error > max_rms) {
					max_rms = error;
					index = i;
				}
			}
		}

		error_vectors[(int)spin][index] = error_vector;
		fock_matrices[(int)spin][index] = fock_matrix;
		b_matrix[(int)spin].at(index, index) = new_rms;

		for (int i = 0; i < index; ++i) {
			b_matrix[(int)spin].at(i, index) = dot(error_vectors[(int)spin][i], error_vector);
		}
		for (int i = index + 1; i < vector_count; ++i) {
			b_matrix[(int)spin].at(index, i) = dot(error_vectors[(int)spin][i], error_vector);
		}
	}

	double getRMS() {
		if (has_rms) return rms;

		computeErrorVector(asymmetric_buffer[0], fock_matrix[0], Spin::alpha);

		rms = dot(asymmetric_buffer[0], asymmetric_buffer[0]);

		if (spin_treatment == SpinTreatment::unrestricted) {
			computeErrorVector(asymmetric_buffer[0], fock_matrix[1], Spin::beta);

			rms += dot(asymmetric_buffer[0], asymmetric_buffer[0]);
			rms * 0.5;
		}
		
		return rms;
	}

	void solve(Spin spin) {
		int vector_count = error_vectors[(int)spin].size();

		for (int i = 0; i < vector_count; ++i) {
			coefficients[(int)spin][i] = 0.0;
		}
		coefficients[(int)spin][vector_count] = 1.0;

		solveLinear(buffer[0] = b_matrix[(int)spin], coefficients[(int)spin], coefficients[(int)spin]);

		buffer[0].resize(matrix_size);
	}

	flo::SymmetricMatrixNd& getResult(Spin spin) {
		int vector_count = error_vectors[(int)spin].size();

		if (vector_count <= 1) return fock_matrices[(int)spin][0];

		SymmetricMatrixNd& total = buffer[0];
		total = 0.0;

		for (int i = 0; i < vector_count; ++i) {
			total = (buffer[2] = total) + (buffer[1] = coefficients[(int)spin][i] * fock_matrices[(int)spin][i]);
		}
		return total;
	}

	void solve() {
		solve(Spin::alpha);
		if (spin_treatment == SpinTreatment::unrestricted) solve(Spin::beta);
	}
}
}