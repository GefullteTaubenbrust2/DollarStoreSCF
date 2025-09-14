#include "Molecule.hpp"
#include "Constants.hpp"

namespace flo {
	Atom::Atom(const vec3& position, Element element, double mass) :
	position(position), element(element), charge((uint)element), mass(mass) {
	}

	Atom::Atom(const vec3& position, Element element) :
	position(position), element(element), charge((uint)element), mass(common_atomic_masses[(int)element] * Da_to_me) {
	}

	size_t Molecule::size() {
		return atoms.size();
	}
	Atom& Molecule::operator[](size_t index) {
		return atoms[index];
	}

	const Atom& Molecule::operator[](size_t index) const {
		return atoms[index];
	}
}