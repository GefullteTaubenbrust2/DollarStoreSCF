#pragma once
#include "../GTO.hpp"
#include "../Molecule.hpp"

namespace scf {
	void assignMolecule(const flo::Molecule& molecule, int charge, uint multiplicity);

	void assignBasis(const flo::BasisSet& basis_set, uint atom_index);

	void assignBasis(const flo::BasisSet& basis_set, flo::Element element);

	void assignBasis(const flo::BasisSet& basis_set);

	void constructBasis();
}