#pragma once
#include "Complex.hpp"

namespace flo {
	template<typename T> struct Vector;

	template<typename T> struct MatrixBase;
	template<typename T> struct Matrix;
	template<typename T> struct UpperTriMatrix;
	template<typename T> struct LowerTriMatrix;
	template<typename T> struct SymmetricMatrix;
	template<typename T> struct HermitianMatrix;
	template<typename T> struct DiagonalMatrix;

	enum MatrixSymmetryFlags;

	namespace intern {
		template<typename T> struct VectorAlgebra;
		template<typename T> struct VectorSum;
		template<typename T> struct VectorDifference;
		template<typename T> struct VectorScalarSum;
		template<typename T> struct VectorScalarDifference;
		template<typename T> struct VectorScalarProduct;

		template<typename T> struct MatrixAlgebra;
		template<typename T> struct MatrixProduct;
		template<typename T> struct MatrixVectorProduct;
		template<typename T> struct MatrixSum;
		template<typename T> struct MatrixDifference;
		template<typename T> struct MatrixScalarProduct;
		template<typename T> struct MatrixTranspose;
		template<typename T> struct MatrixAdjoint;
	}

	typedef Vector<double> VectorNd;
	typedef Vector<c64> VectorNc;

	typedef Matrix<double> MatrixNd;
	typedef Matrix<c64> MatrixNc;
	typedef UpperTriMatrix<double> UpperTriMatrixNd;
	typedef UpperTriMatrix<c64> UpperTriMatrixNc;
	typedef LowerTriMatrix<double> LowerTriMatrixNd;
	typedef LowerTriMatrix<c64> LowerTriMatrixNc;
	typedef SymmetricMatrix<double> SymmetricMatrixNd;
	typedef SymmetricMatrix<c64> SymmetricMatrixNc;
	typedef HermitianMatrix<double> HermitianMatrixNd;
	typedef HermitianMatrix<c64> HermitianMatrixNc;
	typedef DiagonalMatrix<double> DiagonalMatrixNd;
	typedef DiagonalMatrix<c64> DiagonalMatrixNc;
}