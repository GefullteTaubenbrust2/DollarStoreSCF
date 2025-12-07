#include "Frequency.hpp"

#include "../../lalib/Lapack.hpp"

#include "../../util/FormattedStream.hpp"

#include "../Molecule.hpp"
#include "../Constants.hpp"

#include "SCFBasis.hpp"
#include "SCFCommon.hpp"
#include "SCFSolver.hpp"
#include "Gradient.hpp"
#include "Output.hpp"
#include "Energy.hpp"

using namespace flo;

namespace scf {
	extern Molecule molecule;
	VectorNd frequencies;

namespace freq {
	SymmetricMatrixNd nuclear_hessian, k_matrix;
	VectorNd base_gradient;
	VectorNd gradient;
	VectorNd base_coordinates;
	VectorNd coordinates;
	VectorNd base_internal_gradient;
	const MatrixNd* displacement_matrix = nullptr;
	const MatrixNd* wilson_matrix = nullptr;

	double finite_difference = 0.001;
	double step_size = 0.0;

	double equilibrium_total_energy = 0.0;

	bool use_internal_coordinates = true;

	void applyOffset(int coordinate_index) {
		coordinates = base_coordinates;

		if (use_internal_coordinates) {
			double cartesian_offset = 0.0;

			for (int i = 0; i < molecule.size() * 3; ++i) {
				cartesian_offset += displacement_matrix->at(i, coordinate_index) * displacement_matrix->at(i, coordinate_index);
			}

			step_size = finite_difference / std::sqrt(cartesian_offset);

			coordinates[coordinate_index] += step_size;

			molecule.assignInternalCoordinates(coordinates);
		}
		else {
			coordinates[coordinate_index] += finite_difference;

			for (int j = 0; j < molecule.size(); ++j) {
				molecule[j].position.x = coordinates[j * 3];
				molecule[j].position.y = coordinates[j * 3 + 1];
				molecule[j].position.z = coordinates[j * 3 + 2];
			}
		}
	}

	void updateHessian(int coordinate_index) {
		if (use_internal_coordinates) {
			for (int i = 0; i < molecule.size() * 3; ++i) {
				for (int j = 0; j < molecule.size() * 3; ++j) {
					if (i >= j) nuclear_hessian.at(j, i) += 0.5 * gradient[i] * wilson_matrix->at(coordinate_index, j) / step_size;
					if (i <= j) nuclear_hessian.at(i, j) += 0.5 * gradient[i] * wilson_matrix->at(coordinate_index, j) / step_size;
				}
			}
		}
		else {
			for (int j = 0; j <= coordinate_index; ++j) {
				nuclear_hessian.at(j, coordinate_index) += 0.5 * gradient[j] / finite_difference;
			}
			for (int j = coordinate_index; j < gradient.size(); ++j) {
				nuclear_hessian.at(coordinate_index, j) += 0.5 * gradient[j] / finite_difference;
			}
		}
	}
}

using namespace freq;

	// Analytical frequency calculations are faster and more reliable, but also extremely painful to implement.
	// Just think about this - CPHF equations aside, we need the second and mixed derivates of:
	// Overlap, kinetic energy, nuclear potential, 2-center, 3-center and 4-center repulsion integrals.
	// That's 10 new integrals just for this one feature. And these would be more complicated than the gradients,
	// of which the 4-center repulsion gradient is already an utter mess.
	void runNumericalFrequencyCalculation() {
		if (use_internal_coordinates) {
			base_coordinates.resize(molecule.getInternalCoordinateCount());
			molecule.calculateInternalCoordinates(base_coordinates);
		}
		else {
			base_coordinates.resize(molecule.size() * 3);
			for (int i = 0; i < molecule.size(); ++i) {
				base_coordinates[i * 3] = molecule[i].position.x;
				base_coordinates[i * 3 + 1] = molecule[i].position.y;
				base_coordinates[i * 3 + 2] = molecule[i].position.z;
			}
		}
		nuclear_hessian.resize(molecule.size() * 3);
		nuclear_hessian = 0.0;

		displacement_matrix = &molecule.calculateDisplacementMatrix();
		wilson_matrix = &molecule.calculateWilsonMatrix();

		getEnergyGradient(base_gradient);

		base_internal_gradient = base_gradient * *displacement_matrix;

		equilibrium_total_energy = getTotalEnergy();

		uint iteration_count = use_internal_coordinates ? molecule.getInternalCoordinateCount() : molecule.size() * 3;

		for (int i = 0; i < iteration_count; ++i) {
			applyOffset(i);

			constructBasis();

			solveMOs();

			getEnergyGradient(gradient);

			gradient -= base_gradient;

			updateHessian(i);
		}

		MatrixNd& projector = asymmetric_buffer[0];
		projector = *displacement_matrix * *wilson_matrix;
		MatrixNd& projector_T = asymmetric_buffer[1];
		projector_T = T(projector);

		nuclear_hessian = (asymmetric_buffer[2] = (asymmetric_buffer[3] = projector_T * nuclear_hessian) * projector);

		for (int i = 0; i < 4; ++i) asymmetric_buffer[i].resize(matrix_size, matrix_size);

		buffer[0].resize(nuclear_hessian.getWidth());
		asymmetric_buffer[0].resize(nuclear_hessian.getWidth(), nuclear_hessian.getWidth());

		for (int i = 0; i < nuclear_hessian.getWidth(); ++i) {
			for (int j = i; j < nuclear_hessian.getWidth(); ++j) {
				buffer[0].at(i, j) = nuclear_hessian.at(i, j) / std::sqrt(molecule[i / 3].mass * molecule[j / 3].mass);
			}
		}

		frequencies.resize(nuclear_hessian.getWidth());
		
		computeEigenvectors(buffer[0], asymmetric_buffer[0], frequencies);

		buffer[0].resize(matrix_size);
		asymmetric_buffer[0].resize(matrix_size, matrix_size);

		for (int i = 0; i < frequencies.size(); ++i) {
			if (frequencies[i] > 0.0) frequencies[i] = std::sqrt(frequencies[i]);
			else frequencies[i] = -std::sqrt(-frequencies[i]);
		}

		fout.resetRows();
		fout.offsetRight(2);
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout << '|' << '-' << '_' << "Vibrational\nFrequencies" << '|' << '\n';
		if (frequencies[6] < 0.0) {
			fout << '|' << '-' << '_' << "Imaginary Modes" << '|' << '\n';
			fout.resetRows();
			fout.addRow(NumberFormat::crudeFormatPositive(11, -2), TextAlignment::right, 16);
			for (int i = 6; i < frequencies.size(); ++i) {
				if (frequencies[i] > 0.0) break;
				fout << '|' << (-frequencies[i] * au_to_per_cm) << "i 1/cm" << '|' << '\n';
			}
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout << '|' << '-' << '_' << "Real Modes" << '|' << '\n';
		}
		fout.resetRows();
		fout.addRow(NumberFormat::crudeFormatPositive(11, -2), TextAlignment::right, 16);
		for (int i = 6; i < frequencies.size(); ++i) {
			if (frequencies[i] < 0.0) continue;
			fout << '|' << (frequencies[i] * au_to_per_cm) << " 1/cm" << '|' << '\n';
		}
		fout << '-' << '\n';
		fout.resetRows();
		fout << '\n';
	}
}