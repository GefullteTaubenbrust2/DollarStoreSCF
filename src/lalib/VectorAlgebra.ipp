#pragma once
#include "VectorAlgebra.hpp"

namespace flo {
namespace intern {
	template<typename T>
	VectorSum<T>::VectorSum(const Vector<T>& left, const Vector<T>& right) : left(left), right(right) {}

	template<typename T>
	size_t VectorSum<T>::getSize() const {
		return min(left.size(), right.size());
	}

	template<typename T>
	void VectorSum<T>::evaluate(Vector<T>& target) const {
		if (left.size() != right.size()) return;
		target = left;
		target += right;
	}

	template<typename T>
	double VectorSum<T>::evaluate(size_t index) const {
		return left[index] + right[index];
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	VectorDifference<T>::VectorDifference(const Vector<T>& left, const Vector<T>& right) : left(left), right(right) {}

	template<typename T>
	size_t VectorDifference<T>::getSize() const {
		return min(left.size(), right.size());
	}

	template<typename T>
	void VectorDifference<T>::evaluate(Vector<T>& target) const {
		if (left.size() != right.size()) return;
		target = left;
		target -= right;
	}

	template<typename T>
	double VectorDifference<T>::evaluate(size_t index) const {
		return left[index] - right[index];
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	VectorScalarSum<T>::VectorScalarSum(const Vector<T>& left, const T& right) : left(left), right(right) {}

	template<typename T>
	size_t VectorScalarSum<T>::getSize() const {
		return left.size();
	}

	template<typename T>
	void VectorScalarSum<T>::evaluate(Vector<T>& target) const {
		target = left;
		target += right;
	}

	template<typename T>
	double VectorScalarSum<T>::evaluate(size_t index) const {
		return left[index] + right;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	VectorScalarDifference<T>::VectorScalarDifference(const Vector<T>& left, const T& right) : left(left), right(right) {}

	template<typename T>
	size_t VectorScalarDifference<T>::getSize() const {
		return left.size();
	}

	template<typename T>
	void VectorScalarDifference<T>::evaluate(Vector<T>& target) const {
		target = left;
		target -= right;
	}

	template<typename T>
	double VectorScalarDifference<T>::evaluate(size_t index) const {
		return left[index] - right;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	VectorScalarProduct<T>::VectorScalarProduct(const Vector<T>& left, const T& right) : left(left), right(right) {}

	template<typename T>
	size_t VectorScalarProduct<T>::getSize() const {
		return left.size();
	}

	template<typename T>
	void VectorScalarProduct<T>::evaluate(Vector<T>& target) const {
		target = left;
		target *= right;
	}

	template<typename T>
	double VectorScalarProduct<T>::evaluate(size_t index) const {
		return left[index] * right;
	}

	template<typename T>
	WeightedVectorSum<T> VectorScalarProduct<T>::operator+(const Vector<T>& other) const {
		return WeightedVectorSum<T>(left, right, other, 1.0);
	}

	template<typename T>
	WeightedVectorSum<T> VectorScalarProduct<T>::operator-(const Vector<T>& other) const {
		return WeightedVectorSum<T>(left, right, other, -1.0);
	}

	template<typename T>
	WeightedVectorSum<T> VectorScalarProduct<T>::operator+(const VectorScalarProduct<T>& other) const {
		return WeightedVectorSum<T>(left, right, other.left, other.right);
	}

	template<typename T>
	WeightedVectorSum<T> VectorScalarProduct<T>::operator-(const VectorScalarProduct<T>& other) const {
		return WeightedVectorSum<T>(left, right, other.left, -other.right);
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	MatrixVectorProduct<T>::MatrixVectorProduct(const MatrixBase<T>& matrix, const Vector<T>& vector, bool reverse_order) :
		matrix(matrix), vector(vector), reverse_order(reverse_order) {}

	template<typename T>
	size_t MatrixVectorProduct<T>::getSize() const {
		return reverse_order ? matrix.getWidth() : matrix.getHeight();
	}

	template<typename T>
	void MatrixVectorProduct<T>::evaluate(Vector<T>& target) const {
		if (reverse_order) {
			if (vector.size() != matrix.getHeight()) return;
			for (uint x = 0; x < matrix.getWidth(); ++x) {
				T sum = (T)0.0;
				switch (matrix.getUpperSymmetry()) {
				case MatrixSymmetryFlags::asymmetric:
					for (uint i = 0; i < x && i < matrix.getHeight(); ++i) {
						sum += vector[i] * matrix.at(i, x);
					}
					break;
				case MatrixSymmetryFlags::transpose:
					for (uint i = 0; i < x && i < matrix.getHeight(); ++i) {
						sum += vector[i] * matrix.at(x, i);
					}
					break;
				case MatrixSymmetryFlags::adjoint:
					for (uint i = 0; i < x && i < matrix.getHeight(); ++i) {
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
				if (x < matrix.getWidth() && x < matrix.getHeight()) sum += vector[x] * matrix.at(x, x);
				target[x] = sum;
			}
		}
		else {
			if (vector.size() != matrix.getWidth()) return;
			for (uint x = 0; x < matrix.getHeight(); ++x) {
				T sum = (T)0.0;
				switch (matrix.getLowerSymmetry()) {
				case MatrixSymmetryFlags::asymmetric:
					for (uint i = 0; i < x && i < matrix.getWidth(); ++i) {
						sum += vector[i] * matrix.at(x, i);
					}
					break;
				case MatrixSymmetryFlags::transpose:
					for (uint i = 0; i < x && i < matrix.getWidth(); ++i) {
						sum += vector[i] * matrix.at(i, x);
					}
					break;
				case MatrixSymmetryFlags::adjoint:
					for (uint i = 0; i < x && i < matrix.getWidth(); ++i) {
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
				if (x < matrix.getWidth() && x < matrix.getHeight()) sum += vector[x] * matrix.at(x, x);
				target[x] = sum;
			}
		}
	}

	template<typename T>
	double MatrixVectorProduct<T>::evaluate(size_t index) const {
		double sum = 0.0;
		if (reverse_order) {
			for (int i = 0; i < matrix.getHeight(); ++i) {
				sum += matrix(i, index) * vector[index];
			}
		}
		else {
			for (int i = 0; i < matrix.getWidth(); ++i) {
				sum += matrix(index, i) * vector[index];
			}
		}
		return sum;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	WeightedVectorSum<T>::WeightedVectorSum(const Vector<T>& left, T left_weight, const Vector<T>& right, T right_weight) : 
	left(left), right(right), left_weight(left_weight), right_weight(right_weight) {}

	template<typename T>
	size_t WeightedVectorSum<T>::getSize() const {
		return min(left.size(), right.size());
	}

	template<typename T>
	void WeightedVectorSum<T>::evaluate(Vector<T>& target) const {
		if (left.size() != right.size()) return;
		for (int i = 0; i < left.size(); ++i) {
			target[i] = left[i] * left_weight + right[i] * right_weight;
		}
	}

	template<typename T>
	double WeightedVectorSum<T>::evaluate(size_t index) const {
		return left[index] * left_weight + right[index] * right_weight;
	}
}
}