#pragma once
#include "LalibDeclarations.hpp"
#include "MathUtil.hpp"

namespace flo {
namespace intern {
	template<typename T>
	struct VectorAlgebra {
		virtual size_t getSize() const = 0;

		virtual void evaluate(Vector<T>& target) const = 0;

		virtual double evaluate(size_t index) const = 0;
	};

	template<typename T>
	struct VectorSum : public VectorAlgebra<T> {
		const Vector<T>& left;
		const Vector<T>& right;

		VectorSum<T>(const Vector<T>& left, const Vector<T>& right);

		virtual size_t getSize() const override;

		virtual void evaluate(Vector<T>& target) const override;

		virtual double evaluate(size_t index) const override;
	};

	template<typename T>
	struct VectorDifference : public VectorAlgebra<T> {
		const Vector<T>& left;
		const Vector<T>& right;

		VectorDifference<T>(const Vector<T>& left, const Vector<T>& right);

		virtual size_t getSize() const override;

		virtual void evaluate(Vector<T>& target) const override;

		virtual double evaluate(size_t index) const override;
	};

	template<typename T>
	struct VectorScalarSum : public VectorAlgebra<T> {
		const Vector<T>& left;
		const T& right;

		VectorScalarSum<T>(const Vector<T>& left, const T& right);

		virtual size_t getSize() const override;

		virtual void evaluate(Vector<T>& target) const override;

		virtual double evaluate(size_t index) const override;
	};

	template<typename T>
	struct VectorScalarDifference : public VectorAlgebra<T> {
		const Vector<T>& left;
		const T& right;

		VectorScalarDifference<T>(const Vector<T>& left, const T& right);

		virtual size_t getSize() const override;

		virtual void evaluate(Vector<T>& target) const override;

		virtual double evaluate(size_t index) const override;
	};

	template<typename T>
	struct VectorScalarProduct : public VectorAlgebra<T> {
		const Vector<T>& left;
		const T& right;

		VectorScalarProduct<T>(const Vector<T>& left, const T& right);

		virtual size_t getSize() const override;

		virtual void evaluate(Vector<T>& target) const override;

		virtual double evaluate(size_t index) const override;

		WeightedVectorSum<T> operator+(const Vector<T>& other) const;

		WeightedVectorSum<T> operator-(const Vector<T>& other) const;

		WeightedVectorSum<T> operator+(const VectorScalarProduct<T>& other) const;

		WeightedVectorSum<T> operator-(const VectorScalarProduct<T>& other) const;
	};

	template<typename T>
	struct MatrixVectorProduct : public VectorAlgebra<T> {
		const MatrixBase<T>& matrix;
		const Vector<T>& vector;
		bool reverse_order = false;

		MatrixVectorProduct<T>(const MatrixBase<T>& matrix, const Vector<T>& vector, bool reverse_order);

		virtual size_t getSize() const override;

		virtual void evaluate(Vector<T>& target) const override;

		virtual double evaluate(size_t index) const override;
	};

	template<typename T>
	struct WeightedVectorSum : public VectorAlgebra<T> {
		const Vector<T>& left;
		const Vector<T>& right;
		T left_weight = 1.0;
		T right_weight = 1.0;

		WeightedVectorSum<T>(const Vector<T>& left, T left_weight, const Vector<T>& right, T right_weight);

		virtual size_t getSize() const override;

		virtual void evaluate(Vector<T>& target) const override;

		virtual double evaluate(size_t index) const override;
	};
}
}