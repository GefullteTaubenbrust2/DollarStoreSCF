#pragma once
#include "Vector.hpp"

namespace flo {
	template<typename T>
	Vector<T>::Vector(size_t size) :
		entries(size) {
	}

	template<typename T>
	Vector<T>::Vector(const std::initializer_list<T>& init) : entries(init) {
	}

	template<typename T>
	Vector<T>::Vector(const Vector<T>& copy) {
		*this = copy;
	}

	template<typename T>
	Vector<T>::Vector<T>(Vector<T>& reference, size_t index, size_t size) {
		if (index < reference.size()) entries.reference(reference.getRawData() + index, min(size, reference.entries.capacity() - index), reference.entries.capacity() - index);
	}

	template<typename T>
	void Vector<T>::reference(Vector<T>& reference, size_t index, size_t size) {
		if (index < reference.size()) entries.reference(reference.getRawData() + index, min(size, reference.entries.capacity() - index), reference.entries.capacity() - index);
	}

	template<typename T>
	Vector<T>& Vector<T>::operator=(const Vector<T>& other) {
		if (&other == this) return *this;
		resize(other.size());
		entries = other.entries;
		return *this;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator=(const T& value) {
		for (int i = 0; i < entries.size(); ++i) {
			entries[i] = value;
		}
		return *this;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator=(const intern::VectorAlgebra<T>& expression) {
		resize(expression.getSize());
		expression.evaluate(*this);
		return *this;
	}

	template<typename T>
	size_t Vector<T>::size() const {
		return entries.size();
	}

	template<typename T>
	T& Vector<T>::operator[](uint index) {
		return entries[index];
	}

	template<typename T>
	const T& Vector<T>::operator[](uint index) const {
		return entries[index];
	}

	template<typename T>
	void Vector<T>::resize(size_t _size) {
		entries.resize(_size);
	}

	template<typename T>
	T Vector<T>::operator*(const Vector<T>& b) const {
		if (size() != b.size()) return T(0.0);
		T result = T(0.0);
		for (uint i = 0; i < size(); ++i) {
			result += cc(entries[i]) * b[i];
		}
		return result;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator+=(const Vector<T>& b) {
		if (size() != b.size()) return *this;
		for (uint i = 0; i < size(); ++i) {
			entries[i] += b[i];
		}
		return *this;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator-=(const Vector<T>& b) {
		if (size() != b.size()) return *this;
		for (uint i = 0; i < size(); ++i) {
			entries[i] -= b[i];
		}
		return *this;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator+=(const intern::VectorAlgebra<T>& expression) {
		if (size() != expression.getSize()) return *this;
		for (uint i = 0; i < size(); ++i) {
			entries[i] += expression.evaluate(i);
		}
		return *this;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator-=(const intern::VectorAlgebra<T>& expression) {
		if (size() != expression.getSize()) return *this;
		for (uint i = 0; i < size(); ++i) {
			entries[i] -= expression.evaluate(i);
		}
		return *this;

	}

	template<typename T>
	Vector<T>& Vector<T>::operator+=(const T& b) {
		for (uint i = 0; i < size(); ++i) {
			entries[i] += b;
		}
		return *this;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator-=(const T& b) {
		for (uint i = 0; i < size(); ++i) {
			entries[i] -= b;
		}
		return *this;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator*=(const T& b) {
		for (uint i = 0; i < size(); ++i) {
			entries[i] *= b;
		}
		return *this;
	}

	template<typename T>
	Vector<T>& Vector<T>::operator/=(const T& b) {
		for (uint i = 0; i < size(); ++i) {
			entries[i] /= b;
		}
		return *this;
	}

	template<typename T>
	intern::VectorSum<T> Vector<T>::operator+(const Vector<T>& other) const {
		return intern::VectorSum<T>(*this, other);
	}
	template<typename T>

	intern::VectorDifference<T> Vector<T>::operator-(const Vector<T>& other) const {
		return intern::VectorDifference<T>(*this, other);
	}

	template<typename T>
	intern::VectorScalarSum<T> Vector<T>::operator+(const T& other) const {
		return intern::VectorScalarSum<T>(*this, other);
	}

	template<typename T>
	intern::VectorScalarDifference<T> Vector<T>::operator-(const T& other) const {
		return intern::VectorScalarDifference<T>(*this, other);
	}

	template<typename T>
	intern::VectorScalarProduct<T> Vector<T>::operator*(const T& other) const {
		return intern::VectorScalarProduct<T>(*this, other);
	}

	template<typename T>
	intern::VectorScalarProduct<T> Vector<T>::operator/(const T& other) const {
		return intern::VectorScalarProduct<T>(*this, 1.0 / other);
	}

	template<typename T>
	intern::VectorScalarSum<T> operator+(const T& scalar, const Vector<T>& vector) {
		return intern::VectorScalarSum<T>(vector, scalar);
	}

	template<typename T>
	intern::VectorScalarDifference<T> operator-(const T& scalar, const Vector<T>& vector) {
		return intern::VectorScalarDifference<T>(vector, scalar);
	}

	template<typename T>
	intern::VectorScalarProduct<T> operator*(const T& scalar, const Vector<T>& vector) {
		return intern::VectorScalarProduct<T>(vector, scalar);
	}

	template<typename T>
	intern::VectorScalarProduct<T> operator/(const T& scalar, const Vector<T>& vector) {
		return intern::VectorScalarProduct<T>(vector, 1.0 / scalar);
	}

	template<typename T>
	intern::MatrixVectorProduct<T> Vector<T>::operator*(const MatrixBase<T>& other) const {
		return intern::MatrixVectorProduct<T>(other, *this, true);
	}

	template<typename T>
	intern::WeightedVectorSum<T> Vector<T>::operator+(const intern::VectorScalarProduct<T>& other) const {
		return intern::WeightedVectorSum<T>(*this, 1.0, other.left, other.right);
	}

	template<typename T>
	intern::WeightedVectorSum<T> Vector<T>::operator-(const intern::VectorScalarProduct<T>& other) const {
		return intern::WeightedVectorSum<T>(*this, 1.0, other.left, -other.right);
	}

	template<typename T>
	T* Vector<T>::getRawData() const {
		return entries.getPtr();
	}

	template<typename T>
	T length2(const Vector<T>& x) {
		T sum = (T)0.0;
		for (int i = 0; i < x.size(); ++i) {
			sum += cc(x[i]) * x[i];
		}
		return sum;
	}

	template<typename T>
	T length(const Vector<T>& x) {
		return std::sqrt(length2(x));
	}

	template<typename T>
	T dot(const Vector<T>& x, const Vector<T>& y) {
		T sum = (T)0.0;
		for (int i = 0; i < x.size(); ++i) {
			sum += x[i] * y[i];
		}
		return sum;
	}

	template<typename T>
	Vector<T>& cc(const Vector<T>& x, Vector<T>& target) {
		target.resize(x.size());
		for (uint i = 0; i < x.size(); ++i) {
			target[i] = cc(x[i]);
		}
		return target;
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& outs, const Vector<T>& v) {
		outs << '[';
		for (int i = 0; i < (int)v.size() - 1; ++i) {
			outs << v[i] << ", ";
		}
		if (v.size() > 0) outs << v[v.size() - 1];
		return outs << ']';
	}
}