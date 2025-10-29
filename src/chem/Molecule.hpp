#pragma once
#include <vector>
#include "../util/Types.hpp"
#include "Element.hpp"
#include "../lalib/MathUtil.hpp"
#include "../lalib/Lalib.hpp"

namespace flo {
	struct Atom {
		vec3 position;
		uint charge = 1;
		Element element = Element::hydrogen;
		double mass = 1.0;

		Atom() = default;

		Atom(const vec3& position, Element element, double mass);

		Atom(const vec3& position, Element element);

		vec3 getRepulsionGradient(const vec3& position);
	};

	struct Molecule {
		std::vector<Atom> atoms;
		MatrixNd displacements;

		size_t size();

		Atom& operator[](size_t index);

		const Atom& operator[](size_t index) const;
	};
}