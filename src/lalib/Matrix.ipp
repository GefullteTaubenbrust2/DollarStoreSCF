#pragma once
#include "Matrix.hpp"

namespace flo {
	template<typename T>
	flo::MatrixBase<T>::MatrixBase(size_t data_size) : data(data_size) {}

	template<typename T>
	MatrixSymmetryFlags MatrixBase<T>::getLowerSymmetry() const { return MatrixSymmetryFlags::asymmetric; }

	template<typename T>
	MatrixSymmetryFlags MatrixBase<T>::getUpperSymmetry() const { return MatrixSymmetryFlags::asymmetric; }

	template<typename T>
	intern::MatrixProduct<T> MatrixBase<T>::operator*(const MatrixBase<T>& other) const {
		return intern::MatrixProduct<T>(*this, other);
	}

	template<typename T>
	intern::MatrixSum<T> MatrixBase<T>::operator+(const MatrixBase<T>& other) const {
		return intern::MatrixSum<T>(*this, other);
	}

	template<typename T>
	intern::MatrixDifference<T> MatrixBase<T>::operator-(const MatrixBase<T>& other) const {
		return intern::MatrixDifference<T>(*this, other);
	}

	template<typename T>
	intern::MatrixScalarProduct<T> MatrixBase<T>::operator*(const T& other) const {
		return intern::MatrixScalarProduct<T>(*this, other);
	}

	template<typename T>
	intern::MatrixScalarProduct<T> operator*(const T& scalar, const MatrixBase<T>& matrix) {
		return intern::MatrixScalarProduct<T>(matrix, scalar);
	}

	template<typename T>
	intern::MatrixVectorProduct<T> MatrixBase<T>::operator*(const Vector<T>& other) const {
		return intern::MatrixVectorProduct<T>(*this, other, false);
	}

	template<typename T>
	T* MatrixBase<T>::getRawData() const {
		return data.getPtr();
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	Matrix<T>::Matrix(size_t width, size_t height) :
		width(width), height(height), MatrixBase<T>(width* height) {
	}

	template<typename T>
	Matrix<T>::Matrix(size_t width, size_t height, const std::initializer_list<T>& init) :
		width(width), height(height), MatrixBase<T>(width* height) {
		uint i = 0;
		for (T it : init) {
			Matrix<T>::data[i] = it;
			++i;
			if (i >= Matrix<T>::data.size()) break;
		}
	}

	template<typename T>
	Matrix<T>::Matrix(const Matrix& copy) {
		*this = copy;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator=(const Matrix& other) {
		if (&other == this) return *this;

		width = other.width;
		height = other.height;
		Matrix<T>::data = other.data;
		return *this;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator=(const UpperTriMatrix<T>& other) {
		resize(other.size, other.size);
		for (uint x = 0; x < width; ++x) {
			for (uint y = 0; y <= x; ++y) {
				at(y, x) = other.at(y, x);
			}
			for (uint y = x + 1; y < height; ++y) {
				at(y, x) = (T)0.0;
			}
		}
		return *this;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator=(const LowerTriMatrix<T>& other) {
		resize(other.getWidth(), other.getWidth());
		for (uint x = 0; x < width; ++x) {
			for (uint y = x; y <= height; ++y) {
				at(y, x) = other.at(y, x);
			}
			for (uint y = 0; y < x; ++y) {
				at(y, x) = (T)0.0;
			}
		}
		return *this;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator=(const DiagonalMatrix<T>& other) {
		resize(other.getWidth(), other.getWidth());
		for (uint x = 0; x < width; ++x) {
			for (uint y = 0; y <= height; ++y) {
				at(y, x) = (T)0.0;
			}
			at(x, x) = other.at(x, x);
		}
		return *this;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator=(const intern::MatrixProduct<T>& expression) {
		resize(expression.getWidth(), expression.getHeight());
		for (uint x = 0; x < width; ++x) {
			for (uint y = 0; y < height; ++y) {
				at(y, x) = expression.evaluate(y, x);
			}
		}
		return *this;
	}

	template<typename T>
	Matrix<T>& Matrix<T>::operator=(const intern::MatrixAlgebra<T>& expression) {
		resize(expression.getWidth(), expression.getHeight());
		for (uint x = 0; x < width; ++x) {
			expression.evaluateColumn(*this, 0, height, x);
		}
		return *this;
	}

	template<typename T>
	size_t Matrix<T>::getWidth() const {
		return width;
	}

	template<typename T>
	size_t Matrix<T>::getHeight() const {
		return height;
	}

	template<typename T>
	T Matrix<T>::operator()(uint y, uint x) const {
		if (x >= width || y >= height) return (T)0.0;
		return Matrix<T>::data[x * height + y];
	}

	template<typename T>
	T& Matrix<T>::at(uint y, uint x) {
		return Matrix<T>::data[x * height + y];
	}

	template<typename T>
	const T& Matrix<T>::at(uint y, uint x) const {
		return Matrix<T>::data[x * height + y];
	}

	template<typename T>
	void Matrix<T>::resize(size_t _width, size_t _height) {
		width = _width;
		height = _height;
		Matrix<T>::data.resize(width * height);
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	UpperTriMatrix<T>::UpperTriMatrix(size_t size) {
		resize(size);
	}

	template<typename T>
	UpperTriMatrix<T>::UpperTriMatrix(size_t size, const std::initializer_list<T>& init) {
		resize(size);
		uint i = 0;
		for (T it : init) {
			UpperTriMatrix<T>::data[i] = it;
			++i;
			if (i >= UpperTriMatrix<T>::data.size()) break;
		}
	}

	template<typename T>
	UpperTriMatrix<T>::UpperTriMatrix(const UpperTriMatrix<T>& other) {
		*this = other;
	}

	template<typename T>
	UpperTriMatrix<T>& UpperTriMatrix<T>::operator=(const UpperTriMatrix<T>& other) {
		if (&other == this) return *this;

		resize(other.size);
		UpperTriMatrix<T>::data = other.data;
		return *this;
	}

	template<typename T>
	UpperTriMatrix<T>& UpperTriMatrix<T>::operator=(const Matrix<T>& other) {
		resize(min(other.getWidth(), other.getHeight()));
		for (uint x = 0; x < size; ++x) {
			for (uint y = 0; y <= x; ++y) {
				at(y, x) = other.at(y, x);
			}
		}
		return *this;
	}

	template<typename T>
	UpperTriMatrix<T>& UpperTriMatrix<T>::operator=(const LowerTriMatrix<T>& other) {
		resize(other.getWidth());
		for (uint x = 0; x < size; ++x) {
			for (uint y = 0; y <= x; ++y) {
				at(y, x) = (T)0.0;
			}
			at(x, x) = other.at(x, x);
		}
		return *this;
	}

	template<typename T>
	UpperTriMatrix<T>& UpperTriMatrix<T>::operator=(const DiagonalMatrix<T>& other) {
		resize(other.getWidth());
		for (uint x = 0; x < size; ++x) {
			for (uint y = 0; y <= x; ++y) {
				at(y, x) = (T)0.0;
			}
			at(x, x) = other.at(x, x);
		}
		return *this;
	}

	template<typename T>
	UpperTriMatrix<T>& UpperTriMatrix<T>::operator=(const intern::MatrixProduct<T>& expression) {
		resize(min(expression.getWidth(), expression.getHeight()));
		for (uint x = 0; x < size; ++x) {
			for (uint y = 0; y <= x; ++y) {
				at(y, x) = expression.evaluate(y, x);
			}
		}
		return *this;
	}

	template<typename T>
	UpperTriMatrix<T>& UpperTriMatrix<T>::operator=(const intern::MatrixAlgebra<T>& expression) {
		resize(min(expression.getWidth(), expression.getHeight()));
		for (uint x = 0; x < size; ++x) {
			expression.evaluateColumn(*this, 0, x + 1, x);
		}
		return *this;
	}

	template<typename T>
	void UpperTriMatrix<T>::resize(size_t _size) {
		size = _size;
		UpperTriMatrix<T>::data.resize((size * size + size) / 2);
		columns.resize(size);
		for (uint i = 0; i < size; ++i) {
			columns[i] = &UpperTriMatrix<T>::data[(i * i + i) / 2];
		}
	}

	template<typename T>
	size_t UpperTriMatrix<T>::getWidth() const { return size; }

	template<typename T>
	size_t UpperTriMatrix<T>::getHeight() const { return size; }

	template<typename T>
	T& UpperTriMatrix<T>::at(uint y, uint x) {
		return columns[x][y];
	}

	template<typename T>
	const T& UpperTriMatrix<T>::at(uint y, uint x) const {
		return columns[x][y];
	}

	template<typename T>
	T UpperTriMatrix<T>::operator()(uint y, uint x) const {
		if (x >= size || y >= size || x < y) return (T)0.0;
		return columns[x][y];
	}

	template<typename T>
	MatrixSymmetryFlags UpperTriMatrix<T>::getLowerSymmetry() const { return MatrixSymmetryFlags::vanishes; }

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	LowerTriMatrix<T>::LowerTriMatrix(size_t size) {
		resize(size);
	}

	template<typename T>
	LowerTriMatrix<T>::LowerTriMatrix(const LowerTriMatrix<T>& other) {
		*this = other;
	}

	template<typename T>
	LowerTriMatrix<T>::LowerTriMatrix(size_t size, const std::initializer_list<T>& init) {
		resize(size);
		uint i = 0;
		for (T it : init) {
			UpperTriMatrix<T>::data[i] = it;
			++i;
			if (i >= UpperTriMatrix<T>::data.size()) break;
		}
	}

	template<typename T>
	LowerTriMatrix<T>& LowerTriMatrix<T>::operator=(const LowerTriMatrix<T>& other) {
		resize(other.size);
		UpperTriMatrix<T>::data = other.data;
		return *this;
	}

	template<typename T>
	LowerTriMatrix<T>& LowerTriMatrix<T>::operator=(const Matrix<T>& other) {
		resize(min(other.getWidth(), other.getHeight()));
		size_t& size = LowerTriMatrix<T>::size;
		for (uint x = 0; x < size; ++x) {
			for (uint y = x; y <= size; ++y) {
				at(y, x) = other.at(y, x);
			}
		}
		return *this;
	}

	template<typename T>
	LowerTriMatrix<T>& LowerTriMatrix<T>::operator=(const UpperTriMatrix<T>& other) {
		resize(other.getWidth());
		for (uint x = 0; x < LowerTriMatrix<T>::size; ++x) {
			for (uint y = x; y <= LowerTriMatrix<T>::size; ++y) {
				at(y, x) = (T)0.0;
			}
			at(x, x) = other.at(x, x);
		}
		return *this;
	}

	template<typename T>
	LowerTriMatrix<T>& LowerTriMatrix<T>::operator=(const DiagonalMatrix<T>& other) {
		resize(other.getWidth());
		for (uint x = 0; x < LowerTriMatrix<T>::size; ++x) {
			for (uint y = x; y <= LowerTriMatrix<T>::size; ++y) {
				at(y, x) = (T)0.0;
			}
			at(x, x) = other.at(x, x);
		}
		return *this;
	}

	template<typename T>
	LowerTriMatrix<T>& LowerTriMatrix<T>::operator=(const intern::MatrixProduct<T>& expression) {
		resize(min(expression.getWidth(), expression.getHeight()));
		for (uint x = 0; x < LowerTriMatrix<T>::size; ++x) {
			for (uint y = x; y < LowerTriMatrix<T>::size; ++y) {
				at(y, x) = expression.evaluate(y, x);
			}
		}
		return *this;
	}

	template<typename T>
	LowerTriMatrix<T>& LowerTriMatrix<T>::operator=(const intern::MatrixAlgebra<T>& expression) {
		resize(min(expression.getWidth(), expression.getHeight()));
		for (uint x = 0; x < LowerTriMatrix<T>::size; ++x) {
			expression.evaluateColumn(*this, x, LowerTriMatrix<T>::size, x);
		}
		return *this;
	}

	template<typename T>
	void LowerTriMatrix<T>::resize(size_t _size) {
		auto& m_size = LowerTriMatrix<T>::size;
		auto& m_data = LowerTriMatrix<T>::data;
		auto& m_columns = LowerTriMatrix<T>::columns;

		m_size = _size;
		m_data.resize((m_size * m_size + m_size) / 2);
		m_columns.resize(m_size);
		if (m_size) m_columns[0] = &m_data[0];
		for (uint i = 1; i < m_size; ++i) {
			m_columns[i] = m_columns[i - 1] + (m_size - i + 1);
		}
	}

	template<typename T>
	T& LowerTriMatrix<T>::at(uint y, uint x) {
		return LowerTriMatrix<T>::columns[x][y - x];
	}

	template<typename T>
	const T& LowerTriMatrix<T>::at(uint y, uint x) const {
		return LowerTriMatrix<T>::columns[x][y - x];
	}

	template<typename T>
	T LowerTriMatrix<T>::operator()(uint y, uint x) const {
		if (x >= LowerTriMatrix<T>::size || y >= LowerTriMatrix<T>::size || x > y) return (T)0.0;
		return LowerTriMatrix<T>::columns[x][y - x];
	}

	template<typename T>
	MatrixSymmetryFlags LowerTriMatrix<T>::getLowerSymmetry() const { return MatrixSymmetryFlags::asymmetric; }

	template<typename T>
	MatrixSymmetryFlags LowerTriMatrix<T>::getUpperSymmetry() const { return MatrixSymmetryFlags::vanishes; }

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	SymmetricMatrix<T>::SymmetricMatrix(size_t size) : UpperTriMatrix<T>(size) {}

	template<typename T>
	SymmetricMatrix<T>::SymmetricMatrix(size_t size, const std::initializer_list<T>& init) : UpperTriMatrix<T>(size, init) {}

	template<typename T>
	SymmetricMatrix<T>::SymmetricMatrix(const SymmetricMatrix<T>& copy) {
		*this = copy;
	}

	template<typename T>
	SymmetricMatrix<T>& SymmetricMatrix<T>::operator=(const Matrix<T>& other) { UpperTriMatrix<T>::operator=(other); return *this; }

	template<typename T>
	SymmetricMatrix<T>& SymmetricMatrix<T>::operator=(const UpperTriMatrix<T>& other) { UpperTriMatrix<T>::operator=(other); return *this; }

	template<typename T>
	SymmetricMatrix<T>& SymmetricMatrix<T>::operator=(const LowerTriMatrix<T>& other) { UpperTriMatrix<T>::operator=(other); return *this; }

	template<typename T>
	SymmetricMatrix<T>& SymmetricMatrix<T>::operator=(const intern::MatrixProduct<T>& expression) { UpperTriMatrix<T>::operator=(expression); return *this; }

	template<typename T>
	SymmetricMatrix<T>& SymmetricMatrix<T>::operator=(const intern::MatrixAlgebra<T>& expression) { UpperTriMatrix<T>::operator=(expression); return *this; }

	template<typename T>
	T SymmetricMatrix<T>::operator()(uint y, uint x) const {
		if (x >= SymmetricMatrix<T>::size || y >= SymmetricMatrix<T>::size) return (T)0.0;
		if (x < y) return SymmetricMatrix<T>::columns[y][x];
		return SymmetricMatrix<T>::columns[x][y];
	}

	template<typename T>
	MatrixSymmetryFlags SymmetricMatrix<T>::getLowerSymmetry() const { return MatrixSymmetryFlags::transpose; }

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	HermitianMatrix<T>::HermitianMatrix(size_t size) : UpperTriMatrix<T>(size) {}

	template<typename T>
	HermitianMatrix<T>::HermitianMatrix(size_t size, const std::initializer_list<T>& init) : UpperTriMatrix<T>(size, init) {}

	template<typename T>
	HermitianMatrix<T>::HermitianMatrix(const HermitianMatrix<T>& copy) {
		*this = copy;
	}

	template<typename T>
	HermitianMatrix<T>& HermitianMatrix<T>::operator=(const Matrix<T>& other) { UpperTriMatrix<T>::operator=(other); return *this; }

	template<typename T>
	HermitianMatrix<T>& HermitianMatrix<T>::operator=(const UpperTriMatrix<T>& other) { UpperTriMatrix<T>::operator=(other); return *this; }

	template<typename T>
	HermitianMatrix<T>& HermitianMatrix<T>::operator=(const LowerTriMatrix<T>& other) { UpperTriMatrix<T>::operator=(other); return *this; }

	template<typename T>
	HermitianMatrix<T>& HermitianMatrix<T>::operator=(const intern::MatrixProduct<T>& expression) { UpperTriMatrix<T>::operator=(expression); return *this; }

	template<typename T>
	HermitianMatrix<T>& HermitianMatrix<T>::operator=(const intern::MatrixAlgebra<T>& expression) { UpperTriMatrix<T>::operator=(expression); return *this; }

	template<typename T>
	T HermitianMatrix<T>::operator()(uint y, uint x) const {
		if (x >= HermitianMatrix<T>::size || y >= HermitianMatrix<T>::size) return (T)0.0;
		if (x < y) return cc(HermitianMatrix<T>::columns[y][x]);
		return HermitianMatrix<T>::columns[x][y];
	}

	template<typename T>
	MatrixSymmetryFlags HermitianMatrix<T>::getLowerSymmetry() const { return MatrixSymmetryFlags::adjoint; }

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	DiagonalMatrix<T>::DiagonalMatrix(size_t size) {
		resize(size);
	}

	template<typename T>
	DiagonalMatrix<T>::DiagonalMatrix(size_t size, const std::initializer_list<T>& init) {
		resize(size);
		uint i = 0;
		for (T it : init) {
			DiagonalMatrix<T>::data[i] = it;
			++i;
			if (i >= DiagonalMatrix<T>::data.size()) break;
		}
	}

	template<typename T>
	DiagonalMatrix<T>::DiagonalMatrix(const DiagonalMatrix<T>& other) {
		*this = other;
	}

	template<typename T>
	DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(const Matrix<T>& other) {
		resize(min(other.getWidth(), other.getHeight()));
		for (uint i = 0; i < size; ++i) {
			DiagonalMatrix<T>::data[i] = other.at(i, i);
		}
		return *this;
	}

	template<typename T>
	DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(const UpperTriMatrix<T>& other) {
		resize(other.getWidth());
		for (uint i = 0; i < size; ++i) {
			DiagonalMatrix<T>::data[i] = other.at(i, i);
		}
		return *this;
	}

	template<typename T>
	DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(const LowerTriMatrix<T>& other) {
		resize(other.getWidth());
		for (uint i = 0; i < size; ++i) {
			DiagonalMatrix<T>::data[i] = other.at(i, i);
		}
		return *this;
	}

	template<typename T>
	DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(const DiagonalMatrix<T>& other) {
		resize(other.size);
		for (uint i = 0; i < size; ++i) {
			DiagonalMatrix<T>::data[i] = other.data[i];
		}
		return *this;
	}

	template<typename T>
	DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(const intern::MatrixProduct<T>& expression) {
		resize(min(expression.getWidth(), expression.getHeight()));
		for (uint i = 0; i < size; ++i) {
			DiagonalMatrix<T>::data[i] = expression.evaluate(i, i);
		}
		return *this;
	}

	template<typename T>
	DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(const intern::MatrixAlgebra<T>& expression) {
		resize(min(expression.getWidth(), expression.getHeight()));
		for (uint i = 0; i < size; ++i) {
			expression.evaluateColumn(*this, i, i + 1, i);
		}
		return *this;
	}

	template<typename T>
	size_t DiagonalMatrix<T>::getWidth() const {
		return size;
	}

	template<typename T>
	size_t DiagonalMatrix<T>::getHeight() const {
		return size;
	}

	template<typename T>
	T DiagonalMatrix<T>::operator()(uint y, uint x) const {
		if (x >= size || y >= size || x != y) return (T)0.0;
		return DiagonalMatrix<T>::data[x];
	}

	template<typename T>
	T& DiagonalMatrix<T>::at(uint y, uint x) {
		return DiagonalMatrix<T>::data[x];
	}

	template<typename T>
	const T& DiagonalMatrix<T>::at(uint y, uint x) const {
		return DiagonalMatrix<T>::data[x];
	}

	template<typename T>
	void DiagonalMatrix<T>::resize(size_t _size) {
		size = _size;
		DiagonalMatrix<T>::data.resize(_size);
	}

	template<typename T>
	MatrixSymmetryFlags DiagonalMatrix<T>::getLowerSymmetry() const { return MatrixSymmetryFlags::vanishes; }

	template<typename T>
	MatrixSymmetryFlags DiagonalMatrix<T>::getUpperSymmetry() const { return MatrixSymmetryFlags::vanishes; }

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	template<typename T>
	T trace(const MatrixBase<T>& matrix) {
		uint min_size = min(matrix.getWidth(), matrix.getHeight);
		T result = (T)0.0;
		for (uint i = 0; i < min_size; ++i) {
			result += matrix(i, i);
		}
		return result;
	}

	template<typename T>
	T trace(const intern::MatrixProduct<T>& matrix) {
		uint min_size = min(matrix.getWidth(), matrix.getHeight());
		T result = (T)0.0;
		for (uint i = 0; i < min_size; ++i) {
			result += matrix.evaluate(i, i);
		}
		return result;
	}

	template<typename t>
	intern::MatrixTranspose<t> T(const MatrixBase<t>& matrix) {
		return intern::MatrixTranspose<t>(matrix);
	}

	template<typename T>
	intern::MatrixAdjoint<T> dagger(const MatrixBase<T>& matrix) {
		return intern::MatrixAdjoint<T>(matrix);
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& stream, MatrixBase<T>& matrix) {
		stream << '[';
		for (uint y = 0; y < matrix.getHeight(); ++y) {
			stream << '[';
			for (uint x = 0; x < matrix.getWidth() - 1; ++x) {
				stream << matrix(y, x) << ", ";
			}
			if (matrix.getWidth() > 0) stream << matrix(y, matrix.getWidth() - 1);
			stream << ']';
			if (y < matrix.getHeight() - 1) stream << ", ";
		}
		stream << ']';
		return stream;
	}
}