#pragma once
#include <ostream>
#include <cmath>
#include "../util/Types.hpp"

namespace flo {
	template<typename T>
	struct Complex {
		T real = (T)0.0, imag = (T)0.0;

		Complex<T>() = default;

		explicit Complex<T>(T x) : real(x), imag(0) {}
		explicit Complex<T>(double x) : real(T(x)), imag(0) {}

		Complex<T>(const T& real, const T& imag) :
			real(real), imag(imag)
		{}

		Complex<T> operator+(const Complex<T>& b) const {
			return Complex<T>(real + b.real, imag + b.imag);
		}

		Complex<T> operator-(const Complex<T>& b) const {
			return Complex<T>(real - b.real, imag - b.imag);
		}

		Complex<T> operator*(const Complex<T>& b) const {
			return Complex<T>(real * b.real - imag * b.imag, real * b.imag + imag * b.real);
		}

		Complex<T> operator/(const Complex<T>& b) const {
			T denom = b.real * b.real + b.imag * b.imag;
			return Complex<T>((real * b.real + imag * b.imag) / denom, (imag * b.real - real * b.imag) / denom);
		}

		Complex<T> operator+(const T& b) const {
			return Complex<T>(real + b, imag);
		}

		Complex<T> operator-(const T& b) const {
			return Complex<T>(real - b, imag);
		}

		Complex<T> operator*(const T& b) const {
			return Complex<T>(real * b, imag * b);
		}

		Complex<T> operator/(const T& b) const {
			return Complex<T>(real / b, imag / b);
		}

		void operator+=(const Complex<T>& b) {
			real += b.real; 
			imag += b.imag;
		}

		void operator-=(const Complex<T>& b) {
			real -= b.real;
			imag -= b.imag;
		}

		void operator*=(const Complex<T>& b) {
			*this = *this * b;
		}

		void operator/=(const Complex<T>& b) {
			*this = *this / b;
		}

		void operator+=(const T& b) {
			real += b;
		}

		void operator-=(const T& b) {
			real -= b;
		}

		void operator*=(const T& b) {
			real *= b;
			imag *= b;
		}

		void operator/=(const T& b) {
			real /= b;
			imag /= b;
		}
	};
	
	typedef Complex<f32> c32;
	typedef Complex<f64> c64;

	template<typename T>
	Complex<T> operator+(const T& a, const Complex<T>& b) {
		return Complex<T>(a + b.real, b.imag);
	}

	template<typename T>
	Complex<T> operator-(const T& a, const Complex<T>& b) {
		return Complex<T>(a - b.real, b.imag);
	}

	template<typename T>
	Complex<T> operator*(const T& a, const Complex<T>& b) {
		return Complex<T>(a * b.real, a * b.imag);
	}

	template<typename T>
	Complex<T> operator/(const T& a, const Complex<T>& b) {
		T denom = b.real * b.real + b.imag * b.imag;
		return Complex<T>((a * b.real) / denom, (-a * b.imag) / denom);
	}

	template<typename T>
	T abs2(const Complex<T>& a) {
		return a.real * a.real + a.imag * a.imag;
	}

	template<typename T>
	T abs(const Complex<T>& a) {
		return std::sqrt(a.real * a.real + a.imag * a.imag);
	}

	template<typename T>
	T phase(const Complex<T>& a) {
		return std::atan2(a.imag, a.real);
	}

	template<typename T>
	Complex<T> ln(const Complex<T>& a) {
		return Complex<T>(std::log(abs(a)), phase(a));
	}

	template<typename T>
	Complex<T> exp(const Complex<T>& a) {
		T factor = std::exp(a.real);
		return Complex<T>(factor * std::cos(a.imag), factor * std::sin(a.imag));
	}

	template<typename T>
	Complex<T> log(const Complex<T>& a, const Complex<T>& b) {
		return log(b) / log(a);
	}

	template<typename T>
	Complex<T> pow(const Complex<T>& a, const Complex<T>& b) {
		return exp(b * ln(a));
	}

	template<typename T>
	Complex<T> cc(const Complex<T>& x) {
		return Complex<T>(x.real, -x.imag);
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& outs, const Complex<T>& c) {
		return outs << c.real << '+' << c.imag << 'i';
	}

	float cc(float x);
	double cc(double x);
}