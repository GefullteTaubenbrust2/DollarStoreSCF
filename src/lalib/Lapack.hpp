#pragma once
#include "../util/Types.hpp"
#include "LalibDeclarations.hpp"

namespace flo {
	void computeEigenvectors(SymmetricMatrixNd& A, MatrixNd& Q, VectorNd& R, double* work_buffer = nullptr);

	void computeEigenvalues(SymmetricMatrixNd& A, VectorNd& R, double* work_buffer = nullptr);

	void singularValueDecomposition(MatrixNd& A, MatrixNd& U, VectorNd& sigma, MatrixNd& VT, double* work_buffer = nullptr, int buffer_size = 0);

	void solveLinear(MatrixNd& A, VectorNd& B, VectorNd& x, int* ipiv = nullptr);

	void solveLinear(SymmetricMatrixNd& A, VectorNd& B, VectorNd& x, int* ipiv = nullptr);
}