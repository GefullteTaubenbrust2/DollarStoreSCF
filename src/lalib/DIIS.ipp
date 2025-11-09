#include "DIIS.hpp"
#include "Lapack.hpp"

namespace flo {
	template<typename T, class ErrorType, class ResultType>
	void DIISSolver<T, ErrorType, ResultType>::resize(uint iteration_count) {
		iterations = iteration_count;
		error_vectors.reserve(iterations);
		error_vectors.resize(0);
		result_vectors.reserve(iterations);
		result_vectors.resize(0);
		b_matrix.resize(iterations + 1);
		b_matrix.resize(0);
		coefficients.resize(iterations + 1);
		coefficients.resize(0);
	}

	template<typename T, class ErrorType, class ResultType>
	void DIISSolver<T, ErrorType, ResultType>::addErrorVector(const ResultType& result_vector, const ErrorType& error_vector) {
		double new_rms = dot(error_vector, error_vector);

		rms = new_rms;

		int index = 0;

		int vector_count = error_vectors.size();
		if (error_vectors.size() < iterations) {
			error_vectors.resize(vector_count + 1);
			result_vectors.resize(vector_count + 1);
			b_matrix.resize(vector_count + 2);
			coefficients.resize(vector_count + 2);

			index = vector_count;
			++vector_count;

			for (int i = 0; i < vector_count; ++i) {
				b_matrix.at(i, vector_count) = 1.0;
			}
			b_matrix.at(vector_count, vector_count) = 0.0;
		}
		else {
			double max_rms = 0.0;
			for (int i = 0; i < iterations; ++i) {
				double error = b_matrix.at(i, i);
				if (error > max_rms) {
					max_rms = error;
					index = i;
				}
			}
		}

		error_vectors[index] = error_vector;
		result_vectors[index] = result_vector;
		b_matrix.at(index, index) = new_rms;

		for (int i = 0; i < index; ++i) {
			b_matrix.at(i, index) = dot(error_vectors[i], error_vector);
		}
		for (int i = index + 1; i < vector_count; ++i) {
			b_matrix.at(index, i) = dot(error_vectors[i], error_vector);
		}
	}

	template<typename T, class ErrorType, class ResultType>
	void DIISSolver<T, ErrorType, ResultType>::solve(SymmetricMatrix<T>* buffer) {
		SymmetricMatrix<T> buf;
		if (!buffer) buffer = &buf;

		int vector_count = error_vectors.size();

		for (int i = 0; i < vector_count; ++i) {
			coefficients[i] = 0.0;
		}
		coefficients[vector_count] = 1.0;

		solveLinear(*buffer = b_matrix, coefficients, coefficients);
	}

	template<typename T, class ErrorType, class ResultType>
	T DIISSolver<T, ErrorType, ResultType>::getRMS() {
		return rms;
	}

	template<typename T, class ErrorType, class ResultType>
	void DIISSolver<T, ErrorType, ResultType>::calculateResult(ResultType& result) {
		int vector_count = error_vectors.size();

		if (!vector_count) return;

		if (vector_count <= 1) {
			result = result_vectors[0];
			return;
		}

		result = coefficients[0] * result_vectors[0];

		for (int i = 1; i < vector_count; ++i) {
			result = result + coefficients[i] * result_vectors[i];
		}
	}
}