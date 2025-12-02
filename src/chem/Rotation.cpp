#include "Rotation.hpp"
#include "../lalib/Lapack.hpp"

namespace flo {
	vec3 calculateMomentOfInertia(const Molecule& mol, MatrixNd& vectors) {
		SymmetricMatrixNd inertia_matrix(3);
		inertia_matrix = 0.0;

		vec3 center_of_mass(0.0);
		double total_mass = 0.0;

		for (int i = 0; i < mol.size(); ++i) {
			center_of_mass += mol[i].position * mol[i].mass;
			total_mass += mol[i].mass;
		}

		center_of_mass /= total_mass;

		for (int i = 0; i < mol.size(); ++i) {
			vec3 p = mol[i].position - center_of_mass;
			double mass = mol[i].mass;

			inertia_matrix.at(0, 0) += mass * (p.y * p.y + p.z * p.z);
			inertia_matrix.at(1, 1) += mass * (p.x * p.x + p.z * p.z);
			inertia_matrix.at(2, 2) += mass * (p.x * p.x + p.y * p.y);

			inertia_matrix.at(0, 1) += -mass * p.x * p.y;
			inertia_matrix.at(0, 2) += -mass * p.x * p.z;
			inertia_matrix.at(1, 2) += -mass * p.y * p.z;
		}

		VectorNd eigenvalues(3);

		vectors.resize(3, 3);

		computeEigenvectors(inertia_matrix, vectors, eigenvalues);

		return vec3(eigenvalues[0], eigenvalues[1], eigenvalues[2]);
	}
}