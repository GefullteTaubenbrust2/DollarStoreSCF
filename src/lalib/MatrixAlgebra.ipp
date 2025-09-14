#pragma once
#include "MatrixAlgebra.hpp"

namespace flo {
namespace intern {
	template<typename T>
	T evaluateDiagonal(const MatrixBase<T>& m1, const MatrixBase<T>& m2, MultiplicationFlags flags, uint y, uint x) {
		bool y_within_left = y < m1.getWidth();
		bool x_within_right = x < m2.getHeight();
		if (y == x && y_within_left && x_within_right) return m1.at(y, y) * m2.at(y, y);
		else {
			T sum = (T)0.0;
			if (x_within_right) {
				switch (flags & all_left) {
				case left_ordinary:
					sum += m1.at(y, x) * m2.at(x, x);
					break;
				case left_transpose:
					sum += m1.at(x, y) * m2.at(x, x);
					break;
				case left_adjoint:
					sum += cc(m1.at(x, y)) * m2.at(x, x);
					break;
				}
			}
			if (y_within_left) {
				switch (flags & all_right) {
				case right_ordinary:
					sum += m1.at(y, y) * m2.at(y, x);
					break;
				case right_transpose:
					sum += m1.at(y, y) * m2.at(x, y);
					break;
				case right_adjoint:
					sum += m1.at(y, y) * cc(m2.at(x, y));
					break;
				}
			}
			return sum;
		}
	}

	template<typename T>
	T evaluateSum(const MatrixBase<T>& m1, const MatrixBase<T>& m2, MultiplicationFlags flags, uint y, uint x, uint min, uint max) {
		T sum = (T)0.0;
		switch (flags) {
		case loro:
			for (uint i = min; i < max; ++i) {
				sum += m1.at(y, i) * m2.at(i, x);
			}
			break;
		case ltro:
			for (uint i = min; i < max; ++i) {
				sum += m1.at(i, y) * m2.at(i, x);
			}
			break;
		case laro:
			for (uint i = min; i < max; ++i) {
				sum += cc(m1.at(i, y)) * m2.at(i, x);
			}
			break;
		case lort:
			for (uint i = min; i < max; ++i) {
				sum += m1.at(y, i) * m2.at(x, i);
			}
			break;
		case ltrt:
			for (uint i = min; i < max; ++i) {
				sum += m1.at(i, y) * m2.at(x, i);
			}
			break;
		case lart:
			for (uint i = min; i < max; ++i) {
				sum += cc(m1.at(i, y)) * m2.at(x, i);
			}
			break;
		case lora:
			for (uint i = min; i < max; ++i) {
				sum += m1.at(y, i) * m2.at(x, i);
			}
			break;
		case ltra:
			for (uint i = min; i < max; ++i) {
				sum += m1.at(i, y) * m2.at(x, i);
			}
			break;
		case lara:
			for (uint i = min; i < max; ++i) {
				sum += cc(m1.at(i, y) * m2.at(x, i));
			}
			break;
		default:
			break;
		}
		return sum;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	MatrixProduct<T>::MatrixProduct(const MatrixBase<T>& m1, const MatrixBase<T>& m2) :
		left(m1), right(m2) {
		flags_uu = (MultiplicationFlags)(m1.getUpperSymmetry() | (m2.getUpperSymmetry() << 3));
		flags_ul = (MultiplicationFlags)(m1.getUpperSymmetry() | (m2.getLowerSymmetry() << 3));
		flags_lu = (MultiplicationFlags)(m1.getLowerSymmetry() | (m2.getUpperSymmetry() << 3));
		flags_ll = (MultiplicationFlags)(m1.getLowerSymmetry() | (m2.getLowerSymmetry() << 3));
	}

	template<typename T>
	T MatrixProduct<T>::evaluate(uint y, uint x) const {
		uint max_i = left.getWidth();
		if (right.getHeight() < max_i) max_i = right.getHeight();
		if (y > x) {
			return	evaluateSum<T>(left, right, flags_lu, y, x, 0, x) +
				evaluateSum<T>(left, right, flags_ll, y, x, x + 1, y) +
				evaluateSum<T>(left, right, flags_ul, y, x, y + 1, max_i) +
				evaluateDiagonal<T>(left, right, flags_ll, y, x);
		}
		else {
			return	evaluateSum<T>(left, right, flags_lu, y, x, 0, y) +
				evaluateSum<T>(left, right, flags_uu, y, x, y + 1, x) +
				evaluateSum<T>(left, right, flags_ul, y, x, x + 1, max_i) +
				evaluateDiagonal<T>(left, right, flags_uu, y, x);
		}
	}

	template<typename T>
	size_t MatrixProduct<T>::getWidth() const {
		return right.getWidth();
	}

	template<typename T>
	size_t MatrixProduct<T>::getHeight() const {
		return left.getHeight();
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	void MatrixAlgebra<T>::evaluate(MatrixBase<T>& target) {
		for (uint x = 0; x < target.getWidth(); ++x) {
			uint max_y = min(x, (uint)target.getHeight());
			if (target.getUpperSymmetry() == MatrixSymmetryFlags::asymmetric) {
				evaluateColumn(target, 0, max_y, x);
			}
			evaluateColumn(target, x, x + 1, x);
			if (target.getLowerSymmetry() == MatrixSymmetryFlags::asymmetric) {
				evaluateColumn(target, max_y + 1, target.getHeight(), x);
			}
		}
	}
	
	// ------------------------------------------------------------------------------------------------------------------------------------------- //
	
	template<typename T>
	MatrixSum<T>::MatrixSum(const MatrixBase<T>& left, const MatrixBase<T>& right) :
		left(left), right(right) {}

	template<typename T>
	void MatrixSum<T>::evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const {
		uint mid = min(min(x, ymax), (uint)MatrixSum<T>::left.getHeight());
		uint end = min(ymax, (uint)MatrixSum<T>::left.getHeight());

		if (x < MatrixSum<T>::left.getWidth()) {
			switch (MatrixSum<T>::left.getUpperSymmetry()) {
			case MatrixSymmetryFlags::asymmetric:
				for (uint y = ymin; y < mid; ++y) {
					target.at(y, x) = MatrixSum<T>::left.at(y, x);
				}
				break;
			case MatrixSymmetryFlags::transpose:
				for (uint y = ymin; y < mid; ++y) {
					target.at(y, x) = MatrixSum<T>::left.at(x, y);
				}
				break;
			case MatrixSymmetryFlags::adjoint:
				for (uint y = ymin; y < mid; ++y) {
					target.at(y, x) = cc(MatrixSum<T>::left.at(x, y));
				}
				break;
			default:
				for (uint y = ymin; y < mid; ++y) {
					target.at(y, x) = (T)0.0;
				}
				break;
			}
			switch (MatrixSum<T>::left.getLowerSymmetry()) {
			case MatrixSymmetryFlags::asymmetric:
				for (uint y = mid + 1; y < end; ++y) {
					target.at(y, x) = MatrixSum<T>::left.at(y, x);
				}
				break;
			case MatrixSymmetryFlags::transpose:
				for (uint y = mid + 1; y < end; ++y) {
					target.at(y, x) = MatrixSum<T>::left.at(y, x);
				}
				break;
			case MatrixSymmetryFlags::adjoint:
				for (uint y = mid + 1; y < end; ++y) {
					target.at(y, x) = MatrixSum<T>::left.at(y, x);
				}
				break;
			default:
				for (uint y = mid + 1; y < end; ++y) {
					target.at(y, x) = (T)0.0;
				}
				break;
			}
			if (x >= ymin && x < ymax && x < MatrixSum<T>::left.getHeight()) {
				target.at(x, x) = MatrixSum<T>::left.at(x, x);
			}
		}
		else {
			for (uint y = ymin; y < end; ++y) {
				target.at(y, x) = (T)0.0;
			}
		}

		mid = min(min(x, ymax), (uint)MatrixSum<T>::right.getHeight());
		end = min(ymax, (uint)MatrixSum<T>::right.getHeight());

		switch (MatrixSum<T>::right.getUpperSymmetry()) {
		case MatrixSymmetryFlags::asymmetric:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) += MatrixSum<T>::right.at(y, x);
			}
			break;
		case MatrixSymmetryFlags::transpose:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) += MatrixSum<T>::right.at(x, y);
			}
			break;
		case MatrixSymmetryFlags::adjoint:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) += cc(MatrixSum<T>::right.at(x, y));
			}
			break;
		}
		switch (MatrixSum<T>::right.getLowerSymmetry()) {
		case MatrixSymmetryFlags::asymmetric:
			for (uint y = mid + 1; y < end; ++y) {
				target.at(y, x) += MatrixSum<T>::right.at(y, x);
			}
			break;
		case MatrixSymmetryFlags::transpose:
			for (uint y = mid + 1; y < end; ++y) {
				target.at(y, x) += MatrixSum<T>::right.at(y, x);
			}
			break;
		case MatrixSymmetryFlags::adjoint:
			for (uint y = mid + 1; y < end; ++y) {
				target.at(y, x) += MatrixSum<T>::right.at(y, x);
			}
			break;
		}
		if (x >= ymin && x < ymax && x < MatrixSum<T>::right.getWidth() && x < MatrixSum<T>::right.getHeight()) {
			target.at(x, x) += MatrixSum<T>::right.at(x, x);
		}
	}

	template<typename T>
	size_t MatrixSum<T>::getWidth() const {
		return max(left.getWidth(), right.getWidth());
	}

	template<typename T>
	size_t MatrixSum<T>::getHeight() const {
		return max(left.getHeight(), right.getHeight());
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	MatrixDifference<T>::MatrixDifference(const MatrixBase<T>& left, const MatrixBase<T>& right) :
		left(left), right(right) {}

	template<typename T>
	void MatrixDifference<T>::evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const {
		uint mid = min(min(x, ymax), (uint)MatrixDifference<T>::left.getHeight());
		uint end = min(ymax, (uint)MatrixDifference<T>::left.getHeight());

		if (x < MatrixDifference<T>::left.getWidth()) {
			switch (MatrixDifference<T>::left.getUpperSymmetry()) {
			case MatrixSymmetryFlags::asymmetric:
				for (uint y = ymin; y < mid; ++y) {
					target.at(y, x) = MatrixDifference<T>::left.at(y, x);
				}
				break;
			case MatrixSymmetryFlags::transpose:
				for (uint y = ymin; y < mid; ++y) {
					target.at(y, x) = MatrixDifference<T>::left.at(x, y);
				}
				break;
			case MatrixSymmetryFlags::adjoint:
				for (uint y = ymin; y < mid; ++y) {
					target.at(y, x) = cc(MatrixDifference<T>::left.at(x, y));
				}
				break;
			default:
				for (uint y = ymin; y < mid; ++y) {
					target.at(y, x) = (T)0.0;
				}
				break;
			}
			switch (MatrixDifference<T>::left.getLowerSymmetry()) {
			case MatrixSymmetryFlags::asymmetric:
				for (uint y = mid + 1; y < end; ++y) {
					target.at(y, x) = MatrixDifference<T>::left.at(y, x);
				}
				break;
			case MatrixSymmetryFlags::transpose:
				for (uint y = mid + 1; y < end; ++y) {
					target.at(y, x) = MatrixDifference<T>::left.at(y, x);
				}
				break;
			case MatrixSymmetryFlags::adjoint:
				for (uint y = mid + 1; y < end; ++y) {
					target.at(y, x) = MatrixDifference<T>::left.at(y, x);
				}
				break;
			default:
				for (uint y = mid + 1; y < end; ++y) {
					target.at(y, x) = (T)0.0;
				}
				break;
			}
			if (x >= ymin && x < ymax && x < MatrixDifference<T>::left.getHeight()) {
				target.at(x, x) = MatrixDifference<T>::left.at(x, x);
			}
		}
		else {
			for (uint y = ymin; y < end; ++y) {
				target.at(y, x) = (T)0.0;
			}
		}

		mid = min(min(x, ymax), (uint)MatrixDifference<T>::right.getHeight());
		end = min(ymax, (uint)MatrixDifference<T>::right.getHeight());

		switch (MatrixDifference<T>::right.getUpperSymmetry()) {
		case MatrixSymmetryFlags::asymmetric:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) -= MatrixDifference<T>::right.at(y, x);
			}
			break;
		case MatrixSymmetryFlags::transpose:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) -= MatrixDifference<T>::right.at(x, y);
			}
			break;
		case MatrixSymmetryFlags::adjoint:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) -= cc(MatrixDifference<T>::right.at(x, y));
			}
			break;
		}
		switch (MatrixDifference<T>::right.getLowerSymmetry()) {
		case MatrixSymmetryFlags::asymmetric:
			for (uint y = mid + 1; y < end; ++y) {
				target.at(y, x) -= MatrixDifference<T>::right.at(y, x);
			}
			break;
		case MatrixSymmetryFlags::transpose:
			for (uint y = mid + 1; y < end; ++y) {
				target.at(y, x) -= MatrixDifference<T>::right.at(y, x);
			}
			break;
		case MatrixSymmetryFlags::adjoint:
			for (uint y = mid + 1; y < end; ++y) {
				target.at(y, x) -= MatrixDifference<T>::right.at(y, x);
			}
			break;
		}
		if (x >= ymin && x < ymax && x < MatrixDifference<T>::right.getWidth() && x < MatrixDifference<T>::right.getHeight()) {
			target.at(x, x) -= MatrixDifference<T>::right.at(x, x);
		}
	}

	template<typename T>
	size_t MatrixDifference<T>::getWidth() const {
		return max(left.getWidth(), right.getWidth());
	}

	template<typename T>
	size_t MatrixDifference<T>::getHeight() const {
		return max(left.getHeight(), right.getHeight());
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	MatrixScalarProduct<T>::MatrixScalarProduct(const MatrixBase<T>& matrix, const T& scalar) :
		matrix(matrix), scalar(scalar) {}

	template<typename T>
	void MatrixScalarProduct<T>::evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const {
		for (uint y = ymin; y < ymax; ++y) {
			target.at(y, x) = scalar * matrix.at(y, x);
		}
	}

	template<typename T>
	size_t MatrixScalarProduct<T>::getWidth() const {
		return matrix.getWidth();
	}

	template<typename T>
	size_t MatrixScalarProduct<T>::getHeight() const {
		return matrix.getHeight();
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	MatrixTranspose<T>::MatrixTranspose(const MatrixBase<T>& matrix) : matrix(matrix) {}

	template<typename T>
	void MatrixTranspose<T>::evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const {
		uint mid = clamp(x, ymin, ymax - 1);
		switch (matrix.getLowerSymmetry()) {
		case MatrixSymmetryFlags::asymmetric:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) = matrix.at(x, y);
			}
			break;
		case MatrixSymmetryFlags::transpose:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) = matrix.at(y, x);
			}
			break;
		case MatrixSymmetryFlags::adjoint:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) = cc(matrix.at(y, x));
			}
			break;
		}
		if (x < ymax && x >= ymin) target.at(x, x) = matrix.at(x, x);
		switch (matrix.getLowerSymmetry()) {
		case MatrixSymmetryFlags::asymmetric:
			for (uint y = mid + 1; y < ymax; ++y) {
				target.at(y, x) = matrix.at(x, y);
			}
			break;
		case MatrixSymmetryFlags::transpose:
			for (uint y = mid + 1; y < ymax; ++y) {
				target.at(y, x) = matrix.at(y, x);
			}
			break;
		case MatrixSymmetryFlags::adjoint:
			for (uint y = mid + 1; y < ymax; ++y) {
				target.at(y, x) = cc(matrix.at(y, x));
			}
			break;
		}
	}

	template<typename T>
	size_t MatrixTranspose<T>::getWidth() const {
		return matrix.getHeight();
	}

	template<typename T>
	size_t MatrixTranspose<T>::getHeight() const {
		return matrix.getWidth();
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	MatrixAdjoint<T>::MatrixAdjoint(const MatrixBase<T>& matrix) : matrix(matrix) {}

	template<typename T>
	void MatrixAdjoint<T>::evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const {
		uint mid = clamp(x, ymin, ymax - 1);
		switch (matrix.getLowerSymmetry()) {
		case MatrixSymmetryFlags::asymmetric:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) = cc(matrix.at(x, y));
			}
			break;
		case MatrixSymmetryFlags::transpose:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) = cc(matrix.at(y, x));
			}
			break;
		case MatrixSymmetryFlags::adjoint:
			for (uint y = ymin; y < mid; ++y) {
				target.at(y, x) = matrix.at(y, x);
			}
			break;
		}
		if (x < ymax && x >= ymin) target.at(x, x) = matrix.at(x, x);
		switch (matrix.getLowerSymmetry()) {
		case MatrixSymmetryFlags::asymmetric:
			for (uint y = mid + 1; y < ymax; ++y) {
				target.at(y, x) = cc(matrix.at(x, y));
			}
			break;
		case MatrixSymmetryFlags::transpose:
			for (uint y = mid + 1; y < ymax; ++y) {
				target.at(y, x) = cc(matrix.at(y, x));
			}
			break;
		case MatrixSymmetryFlags::adjoint:
			for (uint y = mid + 1; y < ymax; ++y) {
				target.at(y, x) = matrix.at(y, x);
			}
			break;
		}
	}

	template<typename T>
	size_t MatrixAdjoint<T>::getWidth() const {
		return matrix.getHeight();
	}

	template<typename T>
	size_t MatrixAdjoint<T>::getHeight() const {
		return matrix.getWidth();
	}
}
}