#include "Optimization.h"

#include "../../lalib/Lalib.hpp"

#include "../Molecule.hpp"

#include "SCFSolver.hpp"
#include "SCFBasis.hpp"
#include "Gradient.hpp"
#include "Output.hpp"

using namespace flo;

namespace scf {
	extern Molecule molecule;

	VectorNd gradient;
	VectorNd internal_gradient;
	VectorNd previous_gradient;
	VectorNd descent_step;
	VectorNd internal_coordinates;

	int optimization_iteration = 0;

	double step_size = 0.1;

	bool optimizationConverged() {
		if (optimization_iteration <= 1) return false;
		return step_size < 0.0001;
	}

	void doGradientDescentStep() {
		if (optimization_iteration > 1) {
			if (previous_gradient * internal_gradient < 0.0) step_size *= 0.5;
		}

		descent_step = (-step_size / length(internal_gradient)) * internal_gradient;

		internal_coordinates += descent_step;

		previous_gradient = internal_gradient;
	}

	void runStructureOptimization() {
		gradient.resize(molecule.size() * 3);
		internal_gradient.resize(molecule.size() * 3);
		internal_coordinates.resize(molecule.size() * 3);

		for (int i = 0; i < molecule.size(); ++i) {
			internal_coordinates[3 * i]     = molecule[i].position.x;
			internal_coordinates[3 * i + 1] = molecule[i].position.y;
			internal_coordinates[3 * i + 2] = molecule[i].position.z;
		}

		optimization_iteration = 1;

		while (!optimizationConverged()) {

			constructBasis();

			solveMOs();

			printEnergyContributions();

			getEnergyGradient(gradient);
			
			std::cout << gradient << '\n';

			internal_gradient = gradient;

			doGradientDescentStep();

			for (int i = 0; i < molecule.size(); ++i) {
				molecule[i].position.x = internal_coordinates[3 * i];
				molecule[i].position.y = internal_coordinates[3 * i + 1];
				molecule[i].position.z = internal_coordinates[3 * i + 2];
			}

			++optimization_iteration;
		}
	}
}