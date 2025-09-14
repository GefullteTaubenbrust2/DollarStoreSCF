#pragma once
#include <string>
#include <vector>
#include "../lalib/MathUtil.hpp"
#include "../util/Types.hpp"

namespace flo {
	struct CartesianPolynomialTerm {
		double weight = 0.0;
		uint x = 0, y = 0, z = 0;

		CartesianPolynomialTerm() = default;

		CartesianPolynomialTerm(double weight, uint x, uint y, uint z);

		bool operator&&(const CartesianPolynomialTerm& other) const;
	};

	struct GTOPrimitive {
		double zeta = 0.0;
		double weight = 0.0;

		GTOPrimitive() = default;

		GTOPrimitive(double zeta, double weight);
	};

	struct ContractedGaussian {
		vec3 center;
		std::vector<GTOPrimitive> primitives;
		std::vector<CartesianPolynomialTerm> spherical_harmonic;
		uint l;

		ContractedGaussian() = default;

		ContractedGaussian(const std::vector<GTOPrimitive>& primitives, uint l, int m);

		ContractedGaussian(const std::vector<GTOPrimitive>& primitives, uint l, int m, const vec3& center);

		double operator()(const vec3& position) const;
	};

	struct GTOShell {
		std::vector<GTOPrimitive> primitives;
		uint l = 0;

		GTOShell() = default;

		GTOShell(uint l, const std::vector<GTOPrimitive>& primitives);
	};

	typedef std::vector<GTOShell> AtomicBasis;

	struct BasisSet {
		AtomicBasis atomic_bases[118];
	};

	std::string getOrbitalName(uint l, int m);
}