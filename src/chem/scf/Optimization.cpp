#include "Optimization.h"

#include "../../lalib/Lalib.hpp"
#include "../../lalib/LBFGS.hpp"
#include "../../lalib/DIIS.ipp"

#include "../../util/FormattedStream.hpp"
#include "../../util/Time.hpp"

#include "../Molecule.hpp"

#include "SCFCommon.hpp"
#include "SCFSolver.hpp"
#include "SCFBasis.hpp"
#include "Energy.hpp"
#include "Gradient.hpp"
#include "Output.hpp"

using namespace flo;

namespace scf {
	extern Molecule molecule;

namespace optimization {
	VectorNd gradient;
	VectorNd internal_gradient;
	VectorNd previous_gradient;
	VectorNd descent_step;
	VectorNd internal_coordinates;
	VectorNd previous_coordinates;
	VectorNd error_vector;

	int iteration = 0;

	double step_size = 0.1;
	double delta_E = 0.0;

	uint diis_iterations = 8;
	uint lbfgs_iterations = 8;

	double previous_energy = 0.0;
	double step_threshold = 0.001;
	double energy_threshold = 0.000001;

	DIISSolver<double, VectorNd, VectorNd> optimizer_diis;

	bool use_internal_coordinates = true;

	struct ConvergenceInfo {
		double energy = 0.0;
		double delta_E = 0.0;
		double delta_x = 0.0;
		double gradient = 0.0;

		ConvergenceInfo() = default;

		ConvergenceInfo(double energy, double delta_E, double delta_x, double gradient) :
		energy(energy), delta_E(delta_E), delta_x(delta_x), gradient(gradient) {}
	};

	std::vector<ConvergenceInfo> convergence_info;

	bool optimizationConverged() {
		if (iteration <= 2) return false;

		return (step_size < step_threshold && std::abs(delta_E) < energy_threshold);
	}

	void doGradientDescentStep() {
		if (iteration > 1) {
			if (previous_gradient * internal_gradient < 0.0) step_size *= 0.5;
		}

		descent_step = (-step_size / length(internal_gradient)) * internal_gradient;

		internal_coordinates += descent_step;

		previous_gradient = internal_gradient;
	}

	void calculateGradient() {
		getEnergyGradient(gradient);

		if (use_internal_coordinates) {
			const MatrixNd& displacement_matrix = molecule.calculateDisplacementMatrix();
			internal_gradient = gradient * displacement_matrix;
		}
		else {
			internal_gradient = gradient;
		}
	}

	void applyCoordinates() {
		if (use_internal_coordinates) {
			molecule.assignInternalCoordinates(internal_coordinates, &step_size);
		}
		else {
			step_size = length(descent_step);
			for (int i = 0; i < molecule.size(); ++i) {
				molecule[i].position.x = internal_coordinates[3 * i];
				molecule[i].position.y = internal_coordinates[3 * i + 1];
				molecule[i].position.z = internal_coordinates[3 * i + 2];
			}
		}
	}

	void printConvergence() {
		fout.resetRows();
		fout.offsetRight(2);
		fout.addRow(NumberFormat(), TextAlignment::centered, 93);
		fout << '|' << '-' << '_' << "Structure optimization convergence" << '|' << '\n';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 9);
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout << '|' << '_' << "Iteration" << '|' << '_' << "Energy" << '|' << '_' << "Delta E" << '|' << '_' << "Step size" << '|' << '_' << "Gradient" << '|' << '\n';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::left, 9);
		fout.addRow(NumberFormat::crudeFormat(16, 12), TextAlignment::right, 16);
		fout.addRow(NumberFormat::crudeFormat(16, 12), TextAlignment::right, 16);
		fout.addRow(NumberFormat::crudeFormat(16, 12), TextAlignment::right, 16);
		fout.addRow(NumberFormat::crudeFormat(16, 12), TextAlignment::right, 16);

		for (int i = 0; i < convergence_info.size(); ++i) {
			const ConvergenceInfo& info = convergence_info[i];
			fout << '|' << (i64)(i + 1) << '|' << info.energy << '|' << info.delta_E << '|' << info.delta_x << '|' << info.gradient << '|' << '\n';
		}
		fout << '-' << ',' << '-' << ',' << '-' << ',' << '-' << ',' << '-' << '\n';
		fout.resetRows();
		fout << '\n';
	}
}

	using namespace optimization;

	void runStructureOptimization() {
		if (use_internal_coordinates) {
			molecule.generateInternalCoordinates();
			internal_coordinates.resize(molecule.getInternalCoordinateCount());
			internal_gradient.resize(molecule.getInternalCoordinateCount());
		}
		else {
			internal_coordinates.resize(molecule.size() * 3);
			internal_gradient.resize(molecule.size() * 3);
		}
		gradient.resize(molecule.size() * 3);
		gradient = 0.0;

		if (use_internal_coordinates) {
			molecule.calculateInternalCoordinates(internal_coordinates);
		}
		else {
			for (int i = 0; i < molecule.size(); ++i) {
				internal_coordinates[3 * i] = molecule[i].position.x;
				internal_coordinates[3 * i + 1] = molecule[i].position.y;
				internal_coordinates[3 * i + 2] = molecule[i].position.z;
			}
		}

		iteration = 1;

		useCoreGuess();

		Stopclock clock;

		LBFGSSolver solver(lbfgs_iterations, internal_coordinates.size());

		optimizer_diis.resize(diis_iterations);

		convergence_info.clear();

		previous_energy = 0.0;

		while (true) {

			if (iteration > 1) molecule.printXYZData("New cartesian coordinates");

			constructBasis();

			solveMOs();

			printEnergyContributions();

			delta_E = getTotalEnergy() - previous_energy;
			previous_energy = getTotalEnergy();

			convergence_info.push_back(ConvergenceInfo(previous_energy, delta_E, step_size, length(gradient)));

			if (optimizationConverged()) break;

			calculateGradient();

			error_vector = gradient;

			optimizer_diis.addErrorVector(internal_coordinates, error_vector);

			if (iteration % 2) {
				solver.doStep(internal_gradient, internal_coordinates, descent_step, getTotalEnergy());
			}
			else {
				optimizer_diis.solve(&buffer[0]);

				previous_coordinates = internal_coordinates;

				optimizer_diis.calculateResult(internal_coordinates);

				descent_step = internal_coordinates - previous_coordinates;
			}

			applyCoordinates();

			++iteration;
		}

		fout.resetRows();
		fout.addRow(NumberFormat::crudeFormat(10, 8), TextAlignment::centered, 50);
		fout.centerTable(100);
		fout << '|' << '_' << '-' << "\nStructure optimization has converged \\(^o^)/\nTime taken: " << clock.stop().asSeconds() << " s\n" << '|' << '\n';
		fout.resetRows();
		fout << '\n';

		printConvergence();

		if (use_internal_coordinates) molecule.printInternalCoordinates("Final internal coordinates");
		molecule.printXYZData("Final cartesian coordinates");
	}
}