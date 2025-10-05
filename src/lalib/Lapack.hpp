#pragma once
#include "../util/Types.hpp"
#include "LalibDeclarations.hpp"

namespace flo {
	void computeEigenvectors(SymmetricMatrixNd& A, MatrixNd& Q, VectorNd& R, double* work_buffer = nullptr);

	void solveLinear(MatrixNd& A, VectorNd& B, VectorNd& x, int* ipiv = nullptr);

	void solveLinear(SymmetricMatrixNd& A, VectorNd& B, VectorNd& x, int* ipiv = nullptr);
}