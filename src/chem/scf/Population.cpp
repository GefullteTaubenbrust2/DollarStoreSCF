#include "Population.hpp"

#include "SCFCommon.hpp"
#include "../Molecule.hpp"
#include "../../util/FormattedStream.hpp"

#include <functional>

using namespace flo;

namespace scf {
	extern Molecule molecule;

	DiagonalMatrixNd mulliken_populations[2];
	DiagonalMatrixNd lowdin_populations[2];

	void computeOrbitalPopulations(std::vector<double>& populations, std::vector<std::string>& labels, DiagonalMatrixNd& population_matrix, bool do_labels) {
		uint atom_index = 0;
		std::vector<double> buffer;
		uint l = 99999999;
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
							labels.push_back(molecule.getAtomName(atom_index) + getOrbitalName(li, mi));
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
				fout << '|' << molecule.getAtomName(i) << '|' << (alpha_values[i] + beta_values[i]) << '|' << (alpha_values[i] - beta_values[i]) << '|' << '\n';
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
			fout << '|' << '_' << "Atom" << '|' << '_' << "Population" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::left, 16);
			fout.addRow(NumberFormat::scientificFormat(16, 8), TextAlignment::right, 16);

			computeAtomicPopulations(alpha_values, alpha_populations);
			for (int i = 0; i < alpha_values.size(); ++i) {
				fout << '|' << molecule.getAtomName(i) << '|' << alpha_values[i] << '|' << '\n';
			}
			fout << '-' << ',' << '-' << '\n';
			fout.resetRows();
			fout << '\n';
		}
	}

	void printBondOrders(std::function<double(uint, uint)> calculateBondOrder, const std::string& title) {
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 37);
		fout << '|' << '-' << '_' << title << '|';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout << '|' << '_' << "Bond" << '|' << '_' << "Order" << '|';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::left, 16);
		fout.addRow(NumberFormat::scientificFormat(16, 8), TextAlignment::right, 16);

		for (int i = 0; i < molecule.atoms.size(); ++i) {
			for (int j = 0; j < i; ++j) {
				double distance = length(molecule[i].position - molecule[j].position);
				if (distance < 10.0) {
					fout << '|' << (molecule.getAtomName(j) + '-' + molecule.getAtomName(i)) << '|' << calculateBondOrder(j, i) << '|' << '\n';
				}
			}
		}

		fout << '-' << ',' << '-' << '\n';
		fout.resetRows();
		fout << '\n';
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

	double calculateMullikenBondOrders(uint atom_a, uint atom_b) {
		double bond_order = 0.0;

		for (int mu = atom_basis[atom_a]; mu < atom_basis[atom_a + 1]; ++mu) {
			for (int nu = atom_basis[atom_b]; nu < atom_basis[atom_b + 1]; ++nu) {
				if (spin_treatment == SpinTreatment::unrestricted) bond_order += (density_matrix[0].at(mu, nu) + density_matrix[1].at(mu, nu)) * overlap_matrix.at(mu, nu);
				else bond_order += 2.0 * density_matrix[0].at(mu, nu) * overlap_matrix.at(mu, nu);
			}
		}

		return 2.0 * bond_order;
	}

	void printMullikenBondOrders() {
		printBondOrders(calculateMullikenBondOrders, "Mulliken bond orders");
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

	double calculateWibergBondOrders(uint atom_a, uint atom_b) {
		double bond_order = 0.0;

		SymmetricMatrixNd& orthogonal_density_matrix = buffer[0];

		for (int mu = atom_basis[atom_a]; mu < atom_basis[atom_a + 1]; ++mu) {
			for (int nu = atom_basis[atom_b]; nu < atom_basis[atom_b + 1]; ++nu) {
				bond_order += orthogonal_density_matrix.at(mu, nu) * orthogonal_density_matrix.at(mu, nu);
			}
		}

		return bond_order;
	}

	void printWibergBondOrders() {
		if (spin_treatment == SpinTreatment::unrestricted) buffer[0] = (asymmetric_buffer[0] = lowdin_matrix * (asymmetric_buffer[1] = density_matrix[0] + density_matrix[1])) * lowdin_matrix;
		else buffer[0] = (asymmetric_buffer[0] = lowdin_matrix * density_matrix[0]) * (buffer[1] = 2.0 * lowdin_matrix);
		printBondOrders(calculateWibergBondOrders, "Wiberg bond orders");
	}

	void printMayerValence() {
		calculateMullikenPopulations();

		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 71);
		fout << '|' << '-' << '_' << "Mayer valence indices" << '|';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 8);
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout.addRow(NumberFormat(), TextAlignment::centered, 16);
		fout << '|' << '_' << "Atom" << '|' << '_' << "Total valence" << '|' << '_' << "Bonded valence" << '|' << '_' << "Free valence" << '|';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::left, 8);
		fout.addRow(NumberFormat::crudeFormat(16, 8), TextAlignment::right, 16);
		fout.addRow(NumberFormat::crudeFormat(16, 8), TextAlignment::right, 16);
		fout.addRow(NumberFormat::crudeFormat(16, 8), TextAlignment::right, 16);

		if (spin_treatment == SpinTreatment::unrestricted) {
			asymmetric_buffer[0] = (buffer[0] = density_matrix[0] + density_matrix[1]) * overlap_matrix;
			asymmetric_buffer[1] = (buffer[0] = density_matrix[0] - density_matrix[1]) * overlap_matrix;
		}
		else {
			asymmetric_buffer[0] = (buffer[0] = 2.0 * density_matrix[0]) * overlap_matrix;
			asymmetric_buffer[1] = 0.0;
		}

		MatrixNd& lowdin_total_population = asymmetric_buffer[0];
		MatrixNd& lowdin_spin_population = asymmetric_buffer[1];

		SymmetricMatrixNd bond_order_matrix(molecule.size());

		for (int atom_a = 0; atom_a < molecule.size(); ++atom_a) {
			for (int atom_b = 0; atom_b < atom_a; ++atom_b) {
				bond_order_matrix = calculateWibergBondOrders(atom_b, atom_a);
			}
		}

		for (int atom_a = 0; atom_a < molecule.size(); ++atom_a) {
			double gross_population = 0.0;
			double self_bonding = 0.0;
			double total_bond_order = 0.0;

			for (int mu = atom_basis[atom_a]; mu < atom_basis[atom_a + 1]; ++mu) {
				gross_population += mulliken_populations[0].at(mu, mu);
				if (spin_treatment == SpinTreatment::unrestricted) gross_population += mulliken_populations[1].at(mu, mu);

				for (int nu = atom_basis[atom_a]; nu < atom_basis[atom_a + 1]; ++nu) {
					self_bonding += lowdin_total_population(mu, nu) * lowdin_total_population(nu, mu);
				}
			}

			for (int atom_b = 0; atom_b < molecule.size(); ++atom_b) {
				if (atom_a != atom_b) total_bond_order += bond_order_matrix(atom_a, atom_b);
			}

			double total_valence = 2.0 * gross_population - self_bonding;
			double bonded_valence = total_valence - total_bond_order;
			
			fout << '|' << molecule.getAtomName(atom_a) << '|' << total_valence << '|' << bonded_valence << '|' << total_bond_order << '|' << '\n';
		}

		fout << '-' << ',' << '-' << ',' << '-' << ',' << '-';
		fout.resetRows();
		fout << '\n';
	}

	double calculateMayerBondOrders(uint atom_a, uint atom_b) {
		double bond_order = 0.0;

		MatrixNd& lowdin_total_population = asymmetric_buffer[0];
		MatrixNd& lowdin_spin_population = asymmetric_buffer[1];

		for (int mu = atom_basis[atom_a]; mu < atom_basis[atom_a + 1]; ++mu) {
			for (int nu = atom_basis[atom_b]; nu < atom_basis[atom_b + 1]; ++nu) {
				bond_order += lowdin_total_population.at(mu, nu) * lowdin_total_population(nu, mu) + lowdin_spin_population.at(mu, nu) * lowdin_spin_population(nu, mu);
			}
		}

		return bond_order;
	}

	void printMayerBondOrders() {
		if (spin_treatment == SpinTreatment::unrestricted) {
			asymmetric_buffer[0] = (buffer[0] = density_matrix[0] + density_matrix[1]) * overlap_matrix;
			asymmetric_buffer[1] = (buffer[0] = density_matrix[0] - density_matrix[1]) * overlap_matrix;
		}
		else {
			asymmetric_buffer[0] = (buffer[0] = 2.0 * density_matrix[0]) * overlap_matrix;
			asymmetric_buffer[1] = 0.0;
		}

		printBondOrders(calculateMayerBondOrders, "Mayer bond orders");
	}
}