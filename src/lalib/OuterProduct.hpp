#pragma once
#include "Matrix.h"

namespace flo {
	template<typename T>
	Matrix<T>& outerProduct(const Vector<T>& a, const Vector<T>& b, Matrix<T>& target) {
		target.resize(b.size(), a.size());
		for (uint x = 0; x < b.size(); ++x) {
			for (uint y = 0; y < a.size(); ++y) {
				target.at(y, x) = a[y] * b[x];
			}
		}
		return target;
	}

	template<typename T>
	HermitianMatrix<T>& outerProduct(const Vector<T>& a, HermitianMatrix<T>& target) {
		target.resize(a.size(), a.size());
		for (uint x = 0; x < b.size(); ++x) {
			for (uint y = 0; y <= x; ++y) {
				target.at(y, x) = a[y] * a[x];
			}
		}
		return target;
	}

	template<typename T>
	SymmetricMatrix<T>& outerProduct(const Vector<T>& a, SymmetricMatrix<T>& target) {
		target.resize(a.size(), a.size());
		for (uint x = 0; x < b.size(); ++x) {
			for (uint y = 0; y <= x; ++y) {
				target.at(y, x) = a[y] * a[x];
			}
		}
		return target;
	}
}