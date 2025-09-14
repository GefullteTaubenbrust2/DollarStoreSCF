#pragma once
#include "LalibDeclarations.hpp"
#include "MathUtil.hpp"

namespace flo {
namespace intern {
	template<typename T>
	struct VectorAlgebra {
		virtual void evaluate(Vector<T>& target) const = 0;
	};

	template<typename T>
	struct VectorSum : public VectorAlgebra<T> {
		const Vector<T>& left;
		const Vector<T>& right;

		VectorSum<T>(const Vector<T>& left, const Vector<T>& right);

		virtual void evaluate(Vector<T>& target) const override;
	};

	template<typename T>
	struct VectorDifference : public VectorAlgebra<T> {
		const Vector<T>& left;
		const Vector<T>& right;

		VectorDifference<T>(const Vector<T>& left, const Vector<T>& right);

		virtual void evaluate(Vector<T>& target) const override;
	};

	template<typename T>
	struct VectorScalarSum : public VectorAlgebra<T> {
		const Vector<T>& left;
		const T& right;

		VectorScalarSum<T>(const Vector<T>& left, const T& right);

		virtual void evaluate(Vector<T>& target) const override;
	};

	template<typename T>
	struct VectorScalarDifference : public VectorAlgebra<T> {
		const Vector<T>& left;
		const T& right;

		VectorScalarDifference<T>(const Vector<T>& left, const T& right);

		virtual void evaluate(Vector<T>& target) const override;
	};

	template<typename T>
	struct VectorScalarProduct : public VectorAlgebra<T> {
		const Vector<T>& left;
		const T& right;

		VectorScalarProduct<T>(const Vector<T>& left, const T& right);

		virtual void evaluate(Vector<T>& target) const override;
	};

	template<typename T>
	struct MatrixVectorProduct : public VectorAlgebra<T> {
		const MatrixBase<T>& matrix;
		const Vector<T>& vector;
		bool reverse_order = false;

		MatrixVectorProduct<T>(const MatrixBase<T>& matrix, const Vector<T>& vector, bool reverse_order);

		virtual void evaluate(Vector<T>& target) const override;
	};
}
}