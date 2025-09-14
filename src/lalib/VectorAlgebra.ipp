#pragma once
#include "VectorAlgebra.hpp"

namespace flo {
namespace intern {
	template<typename T>
	VectorSum<T>::VectorSum(const Vector<T>& left, const Vector<T>& right) : left(left), right(right) {}

	template<typename T>
	void VectorSum<T>::evaluate(Vector<T>& target) const {
		if (left.size() != right.size()) return;
		target = left;
		target += right;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	VectorDifference<T>::VectorDifference(const Vector<T>& left, const Vector<T>& right) : left(left), right(right) {}

	template<typename T>
	void VectorDifference<T>::evaluate(Vector<T>& target) const {
		if (left.size() != right.size()) return;
		target = left;
		target -= right;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	VectorScalarSum<T>::VectorScalarSum(const Vector<T>& left, const T& right) : left(left), right(right) {}

	template<typename T>
	void VectorScalarSum<T>::evaluate(Vector<T>& target) const {
		target = left;
		target += right;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	VectorScalarDifference<T>::VectorScalarDifference(const Vector<T>& left, const T& right) : left(left), right(right) {}

	template<typename T>
	void VectorScalarDifference<T>::evaluate(Vector<T>& target) const {
		target = left;
		target -= right;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	VectorScalarProduct<T>::VectorScalarProduct(const Vector<T>& left, const T& right) : left(left), right(right) {}

	template<typename T>
	void VectorScalarProduct<T>::evaluate(Vector<T>& target) const {
		target = left;
		target *= right;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	MatrixVectorProduct<T>::MatrixVectorProduct(const MatrixBase<T>& matrix, const Vector<T>& vector, bool reverse_order) :
		matrix(matrix), vector(vector), reverse_order(reverse_order) {}

	template<typename T>
	void MatrixVectorProduct<T>::evaluate(Vector<T>& target) const {
		target.resize(vector.size());
		if (reverse_order) {
			if (vector.size() != matrix.getHeight()) return;
			for (uint x = 0; x < vector.size(); ++x) {
				T sum = (T)0.0;
				switch (matrix.getUpperSymmetry()) {
				case MatrixSymmetryFlags::asymmetric:
					for (uint i = 0; i < x; ++i) {
						sum += vector[i] * matrix.at(i, x);
					}
					break;
				case MatrixSymmetryFlags::transpose:
					for (uint i = 0; i < x; ++i) {
						sum += vector[i] * matrix.at(x, i);
					}
					break;
				case MatrixSymmetryFlags::adjoint:
					for (uint i = 0; i < x; ++i) {
						sum += vector[i] * cc(matrix.at(x, i));
					}
					break;
				}
				switch (matrix.getLowerSymmetry()) {
				case MatrixSymmetryFlags::asymmetric:
					for (uint i = x + 1; i < matrix.getHeight(); ++i) {
						sum += vector[i] * matrix.at(i, x);
					}
					break;
				case MatrixSymmetryFlags::transpose:
					for (uint i = x + 1; i < matrix.getHeight(); ++i) {
						sum += vector[i] * matrix.at(x, i);
					}
					break;
				case MatrixSymmetryFlags::adjoint:
					for (uint i = x + 1; i < matrix.getHeight(); ++i) {
						sum += vector[i] * cc(matrix.at(x, i));
					}
					break;
				}
				sum += vector[x] * matrix.at(x, x);
				target[x] = sum;
			}
		}
		else {
			if (vector.size() != matrix.getWidth()) return;
			for (uint x = 0; x < vector.size(); ++x) {
				T sum = (T)0.0;
				switch (matrix.getLowerSymmetry()) {
				case MatrixSymmetryFlags::asymmetric:
					for (uint i = 0; i < x; ++i) {
						sum += vector[i] * matrix.at(x, i);
					}
					break;
				case MatrixSymmetryFlags::transpose:
					for (uint i = 0; i < x; ++i) {
						sum += vector[i] * matrix.at(i, x);
					}
					break;
				case MatrixSymmetryFlags::adjoint:
					for (uint i = 0; i < x; ++i) {
						sum += vector[i] * cc(matrix.at(i, x));
					}
					break;
				}
				switch (matrix.getUpperSymmetry()) {
				case MatrixSymmetryFlags::asymmetric:
					for (uint i = x + 1; i < matrix.getWidth(); ++i) {
						sum += vector[i] * matrix.at(x, i);
					}
					break;
				case MatrixSymmetryFlags::transpose:
					for (uint i = x + 1; i < matrix.getWidth(); ++i) {
						sum += vector[i] * matrix.at(i, x);
					}
					break;
				case MatrixSymmetryFlags::adjoint:
					for (uint i = x + 1; i < matrix.getWidth(); ++i) {
						sum += vector[i] * cc(matrix.at(i, x));
					}
					break;
				}
				sum += vector[x] * matrix.at(x, x);
				target[x] = sum;
			}
		}
	}
}
}