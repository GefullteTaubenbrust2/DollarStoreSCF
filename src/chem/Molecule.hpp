#pragma once
#include <vector>
#include "Element.hpp"
#include "../util/Types.hpp"
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
	private:
		std::vector<uint> internal_coordinate_indices;
		std::vector<int> affected_indices;
		MatrixNd displacement_matrix, wilson_matrix;
		vec3 base_position;
		vec3 base_direction;
		vec3 base_normal;

	public:
		std::vector<Atom> atoms;

		Molecule() = default;

		size_t size() const;

		Atom& operator[](size_t index);

		const Atom& operator[](size_t index) const;

		std::string getAtomName(size_t index) const;

		void generateInternalCoordinates(std::vector<uint>* bond_map = nullptr);

		size_t getInternalCoordinateCount() const;

		void calculateInternalCoordinates(VectorNd& internal_coordinates) const;

		void assignInternalCoordinates(const VectorNd& internal_coordinates, double* step_size = nullptr);

		// The 'displacement matrix' is equivalent to the pseudoinverse of Wilson's B matrix.
		const MatrixNd& calculateDisplacementMatrix();

		const MatrixNd& calculateWilsonMatrix();

		// The second derivatives of the internal coordinates as described in
		// [1] B. Vebjorn, H. Trygve, The efficient optimization of molecular geometries using redunant 
		//     internal coordinates, Chem. Phys. 117, 2002, 10.1063/1.1515483.
		void calculateCoordinateHessian(MatrixNd& matrix, size_t coordinate_index);

		void printXYZData(const std::string& table_title) const;

		void printInternalCoordinates(const std::string& table_title) const;
	};
}