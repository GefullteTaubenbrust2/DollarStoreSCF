#pragma once
#include "Lalib.hpp"

namespace flo {
	template<typename T, class ErrorType, class ResultType>
	struct DIISSolver {
	private:
		uint iterations = 0;
		std::vector<ErrorType> error_vectors;
		std::vector<ResultType> result_vectors;
		SymmetricMatrix<T> b_matrix;
		Vector<T> coefficients;

		double rms = 0.0;

	public:
		DIISSolver() = default;

		void resize(uint iterations);

		void addErrorVector(const ResultType& result_vector, const ErrorType& error_vector);

		void solve(SymmetricMatrix<T>* buffer = nullptr);

		T getRMS();

		void calculateResult(ResultType& result);
	};
}