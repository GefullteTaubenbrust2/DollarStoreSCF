#pragma once
#include "Complex.hpp"
#include "LalibDeclarations.hpp"
#include "../util/Pointer.hpp"
#include "MathUtil.hpp"

namespace flo {
	enum MatrixSymmetryFlags {
		vanishes = 0,
		asymmetric = 1,
		transpose = 2,
		adjoint = 4,
	};

	template<typename T>
	intern::MatrixScalarProduct<T> operator*(const T& scalar, const MatrixBase<T>& matrix);

	template<typename T>
	struct MatrixBase {
	protected:
		Array<T> data;

	public:
		MatrixBase<T>() = default;

		MatrixBase<T>(size_t data_size);

		virtual size_t getWidth() const = 0;

		virtual size_t getHeight() const = 0;

		virtual T operator()(uint y, uint x) const = 0;

		virtual T& at(uint y, uint x) = 0;

		virtual const T& at(uint y, uint x) const = 0;

		virtual MatrixSymmetryFlags getLowerSymmetry() const;

		virtual MatrixSymmetryFlags getUpperSymmetry() const;

		virtual intern::MatrixProduct<T> operator*(const MatrixBase<T>& other) const;

		virtual intern::MatrixSum<T> operator+(const MatrixBase<T>& other) const;

		virtual intern::MatrixDifference<T> operator-(const MatrixBase<T>& other) const;

		virtual intern::MatrixScalarProduct<T> operator*(const T& other) const;

		virtual intern::MatrixVectorProduct<T> operator*(const Vector<T>& other) const;

		virtual intern::WeightedMatrixSum<T> operator+(const intern::MatrixScalarProduct<T>& other) const;

		virtual intern::WeightedMatrixSum<T> operator-(const intern::MatrixScalarProduct<T>& other) const;

		T* getRawData() const;
	};

	template<typename T>
	T trace(const MatrixBase<T>& matrix);

	template<typename T>
	T trace(const intern::MatrixProduct<T>& matrix);

	template<typename t>
	intern::MatrixTranspose<t> T(const MatrixBase<t>& matrix);

	template<typename T>
	intern::MatrixAdjoint<T> dagger(const MatrixBase<T>& matrix);

	template<typename T>
	struct Matrix : public MatrixBase<T> {
	protected:
		size_t width = 0, height = 0;

	public:
		Matrix<T>() = default;

		Matrix<T>(size_t width, size_t height);

		Matrix<T>(size_t width, size_t height, const std::initializer_list<T>& init);

		Matrix<T>(const Matrix& copy);

		Matrix<T>& operator=(const Matrix& other);

		Matrix<T>& operator=(const UpperTriMatrix<T>& other);

		Matrix<T>& operator=(const LowerTriMatrix<T>& other);

		Matrix<T>& operator=(const SymmetricMatrix<T>& other);

		Matrix<T>& operator=(const HermitianMatrix<T>& other);

		Matrix<T>& operator=(const DiagonalMatrix<T>& other);

		Matrix<T>& operator=(const intern::MatrixProduct<T>& expression);

		Matrix<T>& operator=(const intern::MatrixAlgebra<T>& expression);

		Matrix<T>& operator=(T value);

		virtual size_t getWidth() const override;

		virtual size_t getHeight() const override;

		virtual T operator()(uint y, uint x) const override;

		virtual T& at(uint y, uint x) override;

		virtual const T& at(uint y, uint x) const override;

		void resize(size_t _width, size_t _height);
	};

	template<typename T>
	struct UpperTriMatrix : public MatrixBase<T> {
	protected:
		Array<T*> columns;
		size_t size = 0;

	public:
		UpperTriMatrix<T>() = default;

		UpperTriMatrix<T>(size_t size);

		UpperTriMatrix<T>(size_t size, const std::initializer_list<T>& init);

		UpperTriMatrix<T>(const UpperTriMatrix<T>& other);

		UpperTriMatrix<T>& operator=(const UpperTriMatrix<T>& other);

		UpperTriMatrix<T>& operator=(const Matrix<T>& other);

		UpperTriMatrix<T>& operator=(const LowerTriMatrix<T>& other);

		UpperTriMatrix<T>& operator=(const DiagonalMatrix<T>& other);

		UpperTriMatrix<T>& operator=(const intern::MatrixProduct<T>& expression);

		UpperTriMatrix<T>& operator=(const intern::MatrixAlgebra<T>& expression);

		UpperTriMatrix<T>& operator=(T value);

		void resize(size_t _size);

		virtual size_t getWidth() const override;

		virtual size_t getHeight() const override;

		virtual T& at(uint y, uint x) override;

		virtual const T& at(uint y, uint x) const override;

		virtual T operator()(uint y, uint x) const override;

		virtual MatrixSymmetryFlags getLowerSymmetry() const override;
	};

	template<typename T>
	struct LowerTriMatrix : public UpperTriMatrix<T> {
		LowerTriMatrix<T>() = default;

		LowerTriMatrix<T>(size_t size);

		LowerTriMatrix<T>(const LowerTriMatrix<T>& other);

		LowerTriMatrix<T>(size_t size, const std::initializer_list<T>& init);

		LowerTriMatrix<T>& operator=(const LowerTriMatrix<T>& other);

		LowerTriMatrix<T>& operator=(const Matrix<T>& other);

		LowerTriMatrix<T>& operator=(const UpperTriMatrix<T>& other);

		LowerTriMatrix<T>& operator=(const DiagonalMatrix<T>& other);

		LowerTriMatrix<T>& operator=(const intern::MatrixProduct<T>& expression);

		LowerTriMatrix<T>& operator=(const intern::MatrixAlgebra<T>& expression);

		LowerTriMatrix<T>& operator=(T value);

		void resize(size_t _size);

		virtual T& at(uint y, uint x) override;

		virtual const T& at(uint y, uint x) const override;

		virtual T operator()(uint y, uint x) const override;

		virtual MatrixSymmetryFlags getLowerSymmetry() const override;

		virtual MatrixSymmetryFlags getUpperSymmetry() const override;
	};

	template<typename T>
	struct SymmetricMatrix : public UpperTriMatrix<T> {
		SymmetricMatrix<T>() = default;

		SymmetricMatrix<T>(size_t size);

		SymmetricMatrix<T>(size_t size, const std::initializer_list<T>& init);

		SymmetricMatrix<T>(const SymmetricMatrix<T>& copy);

		SymmetricMatrix<T>& operator=(const Matrix<T>& other);

		SymmetricMatrix<T>& operator=(const UpperTriMatrix<T>& other);

		SymmetricMatrix<T>& operator=(const LowerTriMatrix<T>& other);

		SymmetricMatrix<T>& operator=(const intern::MatrixProduct<T>& expression);

		SymmetricMatrix<T>& operator=(const intern::MatrixAlgebra<T>& expression);

		SymmetricMatrix<T>& operator=(T value);

		virtual T operator()(uint y, uint x) const override;

		virtual MatrixSymmetryFlags getLowerSymmetry() const override;
	};

	template<typename T>
	struct HermitianMatrix : public UpperTriMatrix<T> {
		HermitianMatrix<T>() = default;

		HermitianMatrix<T>(size_t size);

		HermitianMatrix<T>(size_t size, const std::initializer_list<T>& init);

		HermitianMatrix<T>(const HermitianMatrix<T>& copy);

		HermitianMatrix<T>& operator=(const Matrix<T>& other);

		HermitianMatrix<T>& operator=(const UpperTriMatrix<T>& other);

		HermitianMatrix<T>& operator=(const LowerTriMatrix<T>& other);

		HermitianMatrix<T>& operator=(const intern::MatrixProduct<T>& expression);

		HermitianMatrix<T>& operator=(const intern::MatrixAlgebra<T>& expression);

		HermitianMatrix<T>& operator=(T value);

		virtual T operator()(uint y, uint x) const override;

		virtual MatrixSymmetryFlags getLowerSymmetry() const override;
	};

	template<typename T>
	struct DiagonalMatrix : public MatrixBase<T> {
	protected:
		size_t size = 0;

	public:
		DiagonalMatrix<T>() = default;

		DiagonalMatrix<T>(size_t size);

		DiagonalMatrix<T>(size_t size, const std::initializer_list<T>& init);

		DiagonalMatrix<T>(const DiagonalMatrix<T>& other);

		DiagonalMatrix<T>& operator=(const Matrix<T>& other);

		DiagonalMatrix<T>& operator=(const UpperTriMatrix<T>& other);

		DiagonalMatrix<T>& operator=(const LowerTriMatrix<T>& other);

		DiagonalMatrix<T>& operator=(const DiagonalMatrix<T>& other);

		DiagonalMatrix<T>& operator=(const intern::MatrixProduct<T>& expression);

		DiagonalMatrix<T>& operator=(const intern::MatrixAlgebra<T>& expression);

		DiagonalMatrix<T>& operator=(T value);

		virtual size_t getWidth() const override;

		virtual size_t getHeight() const override;

		virtual T operator()(uint y, uint x) const override;

		virtual T& at(uint y, uint x) override;

		virtual const T& at(uint y, uint x) const override;

		void resize(size_t _size);

		virtual MatrixSymmetryFlags getLowerSymmetry() const override;

		virtual MatrixSymmetryFlags getUpperSymmetry() const override;
	};

	template<typename T>
	std::ostream& operator<<(std::ostream& stream, const MatrixBase<T>& matrix);
}