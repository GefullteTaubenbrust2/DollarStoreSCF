#pragma once
#include <vector>
#include <initializer_list>
#include <ostream>
#include "LalibDeclarations.hpp"
#include "../util/Types.hpp"
#include "../util/Pointer.hpp"

namespace flo {
	template<typename T>
	struct Vector {
	private:
		Array<T> entries;

	public:
		Vector<T>() = default;

		Vector<T>(size_t size);

		Vector<T>(const std::initializer_list<T>& init);

		Vector<T>(const Vector<T>& copy);

		Vector<T>(Vector<T>& reference, size_t index, size_t size);

		void reference(Vector<T>& reference, size_t index, size_t size);

		Vector<T>& operator=(const Vector<T>& other);

		Vector<T>& operator=(const T& value);

		Vector<T>& operator=(const intern::VectorAlgebra<T>& expression);

		size_t size() const;

		T& operator[](uint index);

		const T& operator[](uint index) const;

		void resize(size_t _size);

		T operator*(const Vector<T>& b) const;

		Vector<T>& operator+=(const Vector<T>& b);

		Vector<T>& operator-=(const Vector<T>& b);

		Vector<T>& operator+=(const intern::VectorAlgebra<T>& expression);

		Vector<T>& operator-=(const intern::VectorAlgebra<T>& expression);

		Vector<T>& operator+=(const T& b);

		Vector<T>& operator-=(const T& b);

		Vector<T>& operator*=(const T& b);

		Vector<T>& operator/=(const T& b);

		intern::VectorSum<T> operator+(const Vector<T>& other) const;

		intern::VectorDifference<T> operator-(const Vector<T>& other) const;

		intern::VectorScalarSum<T> operator+(const T& other) const;

		intern::VectorScalarDifference<T> operator-(const T& other) const;

		intern::VectorScalarProduct<T> operator*(const T& other) const;

		intern::VectorScalarProduct<T> operator/(const T& other) const;

		intern::MatrixVectorProduct<T> operator*(const MatrixBase<T>& other) const;

		intern::WeightedVectorSum<T> operator+(const intern::VectorScalarProduct<T>& other) const;

		intern::WeightedVectorSum<T> operator-(const intern::VectorScalarProduct<T>& other) const;

		T* getRawData() const;
	};

	template<typename T>
	T length2(const Vector<T>& x);

	template<typename T>
	T length(const Vector<T>& x);

	template<typename T>
	T dot(const Vector<T>& x, const Vector<T>& y);

	template<typename T>
	Vector<T>& cc(const Vector<T>& x, Vector<T>& target);

	template<typename T>
	std::ostream& operator<<(std::ostream& outs, const Vector<T>& v);
}