#include "Population.hpp"

#include "SCFCommon.hpp"
#include "../Molecule.hpp"
#include "../../util/FormattedStream.hpp"

#include <sstream>

using namespace flo;

namespace scf {
	extern Molecule molecule;

	DiagonalMatrixNd mulliken_populations[2];
	DiagonalMatrixNd lowdin_populations[2];

	void computeOrbitalPopulations(std::vector<double>& populations, std::vector<std::string>& labels, DiagonalMatrixNd& population_matrix, bool do_labels) {
		uint atom_index = 9999999999;
		std::vector<double> buffer;
		uint l = 999999999;
		int m = 0;
		for (int i = 0; i <= matrix_size; ++i) {
			bool update = i == matrix_size;
			if (!update) {
				update |= basis_atoms[i] != atom_index;
			}
			if (update) {
				populations.insert(populations.end(), buffer.begin(), buffer.end());

				if (do_labels) {
					labels.reserve(labels.size() + buffer.size());
					for (uint li = 0;; ++li) {
						bool skip = false;
						for (int mi = -(int)li; mi <= (int)li; ++mi) {
							uint index = li * li + (mi + li);
							if (index >= buffer.size()) {
								skip = true;
								break;
							}
							std::ostringstream oss;
							oss << element_symbols[(int)molecule[atom_index].element - 1] << (atom_index + 1) << ' ' << getOrbitalName(li, mi);
							labels.push_back(oss.str());
						}
						if (skip) break;
					}
				}

				if (i == matrix_size) break;

				atom_index = basis_atoms[i];

				buffer.clear();
			}
			if (basis[i].l != l) {
				l = basis[i].l;
				m = -(int)l;
			}
			if (m > l) m = -(int)l;

			uint orbital_index = l * l + (m + l);
			if (buffer.size() <= orbital_index) buffer.resize((l + 1) * (l + 1));

			buffer[orbital_index] += population_matrix.at(i, i);

			++m;
		}
	}

	void computeAtomicPopulations(std::vector<double>& populations, DiagonalMatrixNd& population_matrix) {
		populations.resize(molecule.size());
		for (int i = 0; i < populations.size(); ++i) populations[i] = 0.0;
		for (int i = 0; i < matrix_size; ++i) {
			populations[basis_atoms[i]] += population_matrix.at(i, i);
		}
	}

	void printPopulations(DiagonalMatrixNd& alpha_populations, DiagonalMatrixNd& beta_populations) {
		std::vector<double> alpha_values, beta_values;
		std::vector<std::string> labels;

		if (spin_treatment == SpinTreatment::unrestricted) {
			fout << '\n' << '|' << '_' << '-' << "By atomic orbitals" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout << '|' << '_' << "Orbital" << '|' << '_' << "Population" << '|' << '_' << "Spin population" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::left, 16);
			fout.addRow(NumberFormat::scientificFormat(16, 8), TextAlignment::right, 16);
			fout.addRow(NumberFormat::scientificFormat(16, 8), TextAlignment::right, 16);

			computeOrbitalPopulations(alpha_values, labels, alpha_populations, true);
			computeOrbitalPopulations(beta_values, labels, beta_populations, false);
			for (int i = 0; i < labels.size(); ++i) {
				if (i) fout << '\n';
				fout << '|' << labels[i] << '|' << (alpha_values[i] + beta_values[i]) << '|' << (alpha_values[i] - beta_values[i]) << '|';
			}

			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 58);
			fout << '|' << '-' << '_' << "By atoms" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout << '|' << '_' << "Orbital" << '|' << '_' << "Population" << '|' << '_' << "Spin population" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::left, 16);
			fout.addRow(NumberFormat::scientificFormat(16, 8), TextAlignment::right, 16);
			fout.addRow(NumberFormat::scientificFormat(16, 8), TextAlignment::right, 16);

			computeAtomicPopulations(alpha_values, alpha_populations);
			computeAtomicPopulations(beta_values, beta_populations);
			for (int i = 0; i < alpha_values.size(); ++i) {
				std::ostringstream oss;
				oss << element_symbols[(int)molecule[i].element - 1] << (i + 1);
				fout << '|' << oss.str() << '|' << (alpha_values[i] + beta_values[i]) << '|' << (alpha_values[i] - beta_values[i]) << '|' << '\n';
			}
			fout << '-' << ',' << '-' << ',' << '-' << '\n';
			fout.resetRows();
			fout << '\n';
		}
		else {
			fout << '\n' << '|' << '_' << '-' << "By atomic orbitals" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout << '|' << '_' << "Orbital" << '|' << '_' << "Population" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::left, 16);
			fout.addRow(NumberFormat::scientificFormat(16, 8), TextAlignment::right, 16);

			computeOrbitalPopulations(alpha_values, labels, alpha_populations, true);
			for (int i = 0; i < labels.size(); ++i) {
				if (i) fout << '\n';
				fout << '|' << labels[i] << '|' << alpha_values[i] << '|';
			}

			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 37);
			fout << '|' << '-' << '_' << "By atoms" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout.addRow(NumberFormat(), TextAlignment::centered, 16);
			fout << '|' << '_' << "Orbital" << '|' << '_' << "Population" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::left, 16);
			fout.addRow(NumberFormat::scientificFormat(16, 8), TextAlignment::right, 16);

			computeAtomicPopulations(alpha_values, alpha_populations);
			for (int i = 0; i < alpha_values.size(); ++i) {
				std::ostringstream oss;
				oss << element_symbols[(int)molecule[i].element - 1] << (i + 1);
				fout << '|' << oss.str() << '|' << alpha_values[i] << '|' << '\n';
			}
			fout << '-' << ',' << '-' << '\n';
			fout.resetRows();
			fout << '\n';
		}
	}

	void calculateMullikenPopulations() {
		mulliken_populations[0].resize(matrix_size);
		if (spin_treatment == SpinTreatment::unrestricted) mulliken_populations[1].resize(matrix_size);

		mulliken_populations[0] = density_matrix[0] * overlap_matrix;
		if (spin_treatment == SpinTreatment::unrestricted) mulliken_populations[1] = density_matrix[1] * overlap_matrix;
		else mulliken_populations[0] = 2.0 * mulliken_populations[0];
	}

	void printMullikenPopulations() {
		fout.offsetRight(2);
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, spin_treatment == SpinTreatment::unrestricted ? 58 : 37);
		fout << '|' << '_' << '-' << "Mulliken Population Analysis" << '|';

		calculateMullikenPopulations();

		printPopulations(mulliken_populations[0], mulliken_populations[1]);
	}
	void calculateLowdinPopulations() {
		lowdin_populations[0].resize(matrix_size);
		if (spin_treatment == SpinTreatment::unrestricted) lowdin_populations[1].resize(matrix_size);

		lowdin_populations[0] = (asymmetric_buffer[0] = lowdin_matrix * density_matrix[0]) * lowdin_matrix;
		if (spin_treatment == SpinTreatment::unrestricted) lowdin_populations[1] = (asymmetric_buffer[0] = lowdin_matrix * density_matrix[1]) * lowdin_matrix;
		else lowdin_populations[0] = 2.0 * lowdin_populations[0];
	}


	void printLowdinPopulations() {
		fout.offsetRight(2);
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, spin_treatment == SpinTreatment::unrestricted ? 58 : 37);
		fout << '|' << '_' << '-' << "Lowdin Population Analysis" << '|';

		calculateLowdinPopulations();

		printPopulations(lowdin_populations[0], lowdin_populations[1]);
	}
}