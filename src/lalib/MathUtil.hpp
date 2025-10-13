#pragma once
#include <cmath>
#include "../util/Types.hpp"

namespace flo {
	#define PI 3.14159265358979

	template<typename T>
	constexpr T max(const T& a, const T& b) {
		return a > b ? a : b;
	}

	template<typename T>
	constexpr T min(const T& a, const T& b) {
		return a < b ? a : b;
	}

	template<typename T>
	constexpr T clamp(const T& x, const T& min, const T& max) {
		return x > max ? max : (x < min ? min : x);
	}

	template<typename T>
	constexpr T factorial(T x) {
		T result = (T)1.0;
		for (T i = x; i > 0; --i) {
			result *= i;
		}
		return result;
	}

	template<typename T>
	constexpr T doubleFactorial(T x) {
		T result = (T)1.0;
		for (T i = x; i > 0; i -= 2) {
			result *= i;
		}
		return result;
	}

	template<typename T>
	constexpr T binomial(T n, T k) {
		T result = (T)1.0;
		for (T i = k; i > 0; --i) {
			result /= i;
		}
		for (T i = n; i > n - k; --i) {
			result *= i;
		}
		return result;
	}

	template<typename T>
	constexpr T pow(T x, uint p) {
		if (p == 1) return x;
		else if (!p) return (T)1.0;
		T y = x;
		uint l = 2;
		for (; l < p; l *= 2) {
			y = y * y;
		}
		return y * pow(x, p - l / 2);
	}

	template<typename T>
	struct Vector2 {
		T x = (T)0.0, y = (T)0.0;

		Vector2<T>() = default;

		Vector2<T>(const T& x) : x(x), y(x) {}

		Vector2<T>(const T& x, const T& y) : x(x), y(y) {}

		Vector2<T> operator+(const Vector2<T>& other) const {
			return Vector2<T>(x + other.x, y + other.y);
		}

		Vector2<T> operator+(const T& other) const {
			return Vector2<T>(x + other, y + other);
		}

		friend Vector2<T> operator+(const T& a, const Vector2<T>& b) {
			return Vector2<T>(a + b.x, a + b.y);
		}

		Vector2<T> operator-(const Vector2<T>& other) const {
			return Vector2<T>(x - other.x, y - other.y);
		}

		Vector2<T> operator-(const T& other) const {
			return Vector2<T>(x - other, y - other);
		}

		friend Vector2<T> operator-(const T& a, const Vector2<T>& b) {
			return Vector2<T>(a - b.x, a - b.y);
		}

		Vector2<T> operator*(const Vector2<T>& other) const {
			return Vector2<T>(x * other.x, y * other.y);
		}

		Vector2<T> operator*(const T& other) const {
			return Vector2<T>(x * other, y * other);
		}

		friend Vector2<T> operator*(const T& a, const Vector2<T>& b) {
			return Vector2<T>(a * b.x, a * b.y);
		}

		Vector2<T> operator/(const Vector2<T>& other) const {
			return Vector2<T>(x / other.x, y / other.y);
		}

		Vector2<T> operator/(const T& other) const {
			return Vector2<T>(x / other, y / other);
		}

		friend Vector2<T> operator/(const T& a, const Vector2<T>& b) {
			return Vector2<T>(a / b.x, a / b.y);
		}

		T& operator[](int i) {
			T* arr[2] = { &x, &y };
			return *arr[i];
		}

		const T& operator[](int i) const {
			const T* arr[2] = { &x, &y };
			return *arr[i];
		}
	};

	template<typename T>
	T dot(const Vector2<T>& a, const Vector2<T>& b) {
		return a.x * b.x + a.y * b.y;
	}

	template<typename T>
	T cross(const Vector2<T>& a, const Vector2<T>& b) {
		return a.x * b.y - a.y * b.x;
	}

	template<typename T>
	T length2(const Vector2<T>& x) {
		return dot(x, x);
	}

	template<typename T>
	T length(const Vector2<T>& x) {
		return std::sqrt(dot(x, x));
	}

	template<typename T>
	struct Vector3 {
		T x = (T)0.0, y = (T)0.0, z = (T)0.0;

		Vector3<T>() = default;

		Vector3<T>(const T& x) : x(x), y(x), z(x) {}

		Vector3<T>(const T& x, const T& y, const T& z) : x(x), y(y), z(z) {}

		Vector3<T> operator+(const Vector3<T>& other) const {
			return Vector3<T>(x + other.x, y + other.y, z + other.z);
		}

		Vector3<T> operator+(const T& other) const {
			return Vector3<T>(x + other, y + other, z + other);
		}

		friend Vector3<T> operator+(const T& a, const Vector3<T>& b) {
			return Vector3<T>(a + b.x, a + b.y, a + b.z);
		}

		Vector3<T> operator-(const Vector3<T>& other) const {
			return Vector3<T>(x - other.x, y - other.y, z - other.z);
		}

		Vector3<T> operator-(const T& other) const {
			return Vector3<T>(x - other, y - other, z - other);
		}

		friend Vector3<T> operator-(const T& a, const Vector3<T>& b) {
			return Vector3<T>(a - b.x, a - b.y, a - b.z);
		}

		Vector3<T> operator*(const Vector3<T>& other) const {
			return Vector3<T>(x * other.x, y * other.y, z * other.z);
		}

		Vector3<T> operator*(const T& other) const {
			return Vector3<T>(x * other, y * other, z * other);
		}

		friend Vector3<T> operator*(const T& a, const Vector3<T>& b) {
			return Vector3<T>(a * b.x, a * b.y, a * b.z);
		}

		Vector3<T> operator/(const Vector3<T>& other) const {
			return Vector3<T>(x / other.x, y / other.y, z / other.z);
		}

		Vector3<T> operator/(const T& other) const {
			return Vector3<T>(x / other, y / other, z / other);
		}

		friend Vector3<T> operator/(const T& a, const Vector3<T>& b) {
			return Vector3<T>(a / b.x, a / b.y, a / b.z);
		}

		T& operator[](int i) {
			T* arr[3] = { &x, &y, &z };
			return *arr[i];
		}

		const T& operator[](int i) const {
			const T* arr[3] = { &x, &y, &z };
			return *arr[i];
		}
	};

	template<typename T>
	T dot(const Vector3<T>& a, const Vector3<T>& b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	template<typename T>
	Vector3<T> cross(const Vector3<T>& a, const Vector3<T>& b) {
		return Vector3<T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
	}

	template<typename T>
	T length2(const Vector3<T>& x) {
		return dot(x, x);
	}

	template<typename T>
	T length(const Vector3<T>& x) {
		return std::sqrt(dot(x, x));
	}

	typedef Vector2<double> vec2;
	typedef Vector3<double> vec3;

	typedef Vector2<int> vec2i;
	typedef Vector3<int> vec3i;
}