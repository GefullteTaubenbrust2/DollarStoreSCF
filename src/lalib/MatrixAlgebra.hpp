#pragma once
#include "LalibDeclarations.hpp"
#include "MathUtil.hpp"
#include <functional>

namespace flo {
namespace intern {
	enum MultiplicationFlags {
		left_ordinary = 1,
		left_transpose = 2,
		left_adjoint = 4,
		right_ordinary = 8,
		right_transpose = 16,
		right_adjoint = 32,

		all_left = 7,
		all_right = 56,

		loro = 9,
		lort = 17,
		lora = 33,
		ltro = 10,
		ltrt = 18,
		ltra = 34,
		laro = 12,
		lart = 20,
		lara = 36,
	};

	template<typename T>
	T evaluateDiagonal(const MatrixBase<T>& m1, const MatrixBase<T>& m2, MultiplicationFlags flags, uint y, uint x);

	template<typename T>
	T evaluateSum(const MatrixBase<T>& m1, const MatrixBase<T>& m2, MultiplicationFlags flags, uint y, uint x, uint min, uint max);

	template<typename T>
	struct MatrixProduct {
		const MatrixBase<T>& left;
		const MatrixBase<T>& right;
		// Definition of regimes: 
		// uu: y < i | x < i
		// ul: y < i | x > i
		// lu: y > i | x < i
		// ll: y > i | x > i
		MultiplicationFlags flags_uu;
		MultiplicationFlags flags_ul;
		MultiplicationFlags flags_lu;
		MultiplicationFlags flags_ll;

		MatrixProduct(const MatrixBase<T>& m1, const MatrixBase<T>& m2);

		T evaluate(uint y, uint x) const;

		size_t getWidth() const;

		size_t getHeight() const;
	};

	template<typename T>
	struct MatrixAlgebra {
		virtual void evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const = 0;

		void evaluate(MatrixBase<T>& target);

		virtual size_t getWidth() const = 0;

		virtual size_t getHeight() const = 0;
	};

	template<typename T>
	struct MatrixSum : public MatrixAlgebra<T> {
		const MatrixBase<T>& left;
		const MatrixBase<T>& right;

		MatrixSum<T>(const MatrixBase<T>& left, const MatrixBase<T>& right);

		virtual void evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const override;

		virtual size_t getWidth() const override;

		virtual size_t getHeight() const override;
	};

	template<typename T>
	struct MatrixDifference : public MatrixAlgebra<T> {
		const MatrixBase<T>& left;
		const MatrixBase<T>& right;

		MatrixDifference<T>(const MatrixBase<T>& left, const MatrixBase<T>& right);

		virtual void evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const override;

		virtual size_t getWidth() const override;

		virtual size_t getHeight() const override;
	};

	template<typename T>
	struct MatrixScalarProduct : public MatrixAlgebra<T> {
		const MatrixBase<T>& matrix;
		const T& scalar;

		MatrixScalarProduct(const MatrixBase<T>& matrix, const T& scalar);

		virtual void evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const override;

		virtual size_t getWidth() const override;

		virtual size_t getHeight() const override;
	};

	template<typename T>
	struct MatrixTranspose : public MatrixAlgebra<T> {
		const MatrixBase<T>& matrix;

		MatrixTranspose<T>(const MatrixBase<T>& matrix);

		virtual void evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const override;

		virtual size_t getWidth() const override;

		virtual size_t getHeight() const override;
	};

	template<typename T>
	struct MatrixAdjoint : public MatrixAlgebra<T> {
		const MatrixBase<T>& matrix;

		MatrixAdjoint<T>(const MatrixBase<T>& matrix);

		virtual void evaluateColumn(MatrixBase<T>& target, uint ymin, uint ymax, uint x) const override;

		virtual size_t getWidth() const override;

		virtual size_t getHeight() const override;
	};
}
}