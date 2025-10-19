#include "WriteOrbitals.hpp"

#include <fstream>

#include "SCFCommon.hpp"

#include "../Molecule.hpp"
#include "../Constants.hpp"

#include "../../util/FormattedStream.hpp"

using namespace flo;

namespace scf {
	extern Molecule molecule;

	void writeOrbitalsMolden(const std::string& path, bool write_virtual) {
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 40);
		fout.centerTable(100);
		fout << '|' << '_' << '-' << "Writing orbitals to molden file" << '|' << '\n';
		fout.resetRows();
		fout << '\n';

		const uint molden_permutations[5][9] = { 
			{ 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			{ 0, 2, 1, 0, 0, 0, 0, 0, 0 },
			{ 2, 1, 3, 0, 4, 0, 0, 0, 0 },
			{ 3, 2, 4, 1, 5, 0, 6, 0, 0 },
			{ 4, 3, 5, 2, 6, 1, 7, 0, 8 }
		};
		const char shell_labels[5] = { 's', 'p', 'd', 'f', 'g' };

		std::ofstream file(path, std::ios::binary);

		file << "[Molden Format]\n";
		file << "[5D]\n";
		file << "[9G]\n";
		file << "[Atoms] AU\n";

		FormattedStream fstream;
		fstream.write_to = &file;

		fstream.addRow(NumberFormat(), TextAlignment::left, 2);
		fstream.addRow(NumberFormat::crudeFormatPositive(4, 4), TextAlignment::right, 4);
		fstream.addRow(NumberFormat::crudeFormatPositive(4, 4), TextAlignment::right, 4);
		fstream.addRow(NumberFormat::crudeFormat(8, 8), TextAlignment::right, 8);
		fstream.addRow(NumberFormat::crudeFormat(8, 8), TextAlignment::right, 8);
		fstream.addRow(NumberFormat::crudeFormat(8, 8), TextAlignment::right, 8);

		for (int i = 0; i < molecule.size(); ++i) {
			Atom& atom = molecule[i];
			fstream << element_symbols[(int)atom.element - 1] << ',' << (i64)(i + 1) << ',' << (i64)atom.element << ',' << atom.position.x << ',' << atom.position.y << ',' << atom.position.z << '\n';
		}

		fstream.resetRows();

		fstream.addRow(NumberFormat::scientificFormat(14, 8), TextAlignment::right, 14);
		fstream.addRow(NumberFormat::scientificFormat(14, 8), TextAlignment::right, 14);

		file << "[GTO]\n";

		for (int i = 0; i < molecule.size(); ++i) {
			file << (i + 1) << " 0\n";

			for (int j = atom_basis[i]; j < atom_basis[i + 1];) {

				ContractedGaussian& cgf = basis[j];

				if (cgf.l < 5) {
					file << shell_labels[cgf.l] << ' ' << cgf.primitives.size() << " 1.00\n";

					for (int k = 0; k < cgf.primitives.size(); ++k) {
						fstream << cgf.primitives[k].zeta << ',' << cgf.primitives[k].unnormalized_weight << '\n';
					}
				}

				j += 2 * cgf.l + 1;
			}

			file << '\n';
		}

		file << "[MO]\n";

		fstream.resetRows();
		fstream.addRow(NumberFormat::crudeFormatPositive(5, 5), TextAlignment::right, 5);
		fstream.addRow(NumberFormat::scientificFormat(14, 8), TextAlignment::right, 14);
		basis;
		for (int i = 0; i < (write_virtual ? matrix_size : electron_count[0]); ++i) {
			file << "Sym=   a" << i << '\n';
			file << "Ene=   " << mo_levels[0][i] << '\n';
			file << "Spin=  Alpha\n";
			file << "Occup= " << ((i < electron_count[0]) * (1 + (spin_treatment == SpinTreatment::restricted))) << '\n';

			int ao_index = 0;
			for (int j = 0; j < molecule.size(); ++j) {
				for (int k = atom_basis[j]; k < atom_basis[j + 1];) {
					ContractedGaussian& cgf = basis[k];
					if (cgf.l < 5) {
						for (int m = 0; m <= 2 * cgf.l; ++m) {
							++ao_index;

							fstream << (i64)ao_index << ',' << coefficient_matrix[0](k + molden_permutations[cgf.l][m], i) << '\n';
						}
					}

					k += 2 * cgf.l + 1;
				}
			}
		}

		if (spin_treatment == SpinTreatment::unrestricted) {
			for (int i = 0; i < (write_virtual ? matrix_size : electron_count[1]); ++i) {
				file << "Sym=   b" << i << '\n';
				file << "Ene=   " << mo_levels[1][i] << '\n';
				file << "Spin=  Beta\n";
				file << "Occup= " << (i < electron_count[1]) << '\n';

				int ao_index = 0;
				for (int j = 0; j < molecule.size(); ++j) {
					for (int k = atom_basis[j]; k < atom_basis[j + 1];) {
						ContractedGaussian& cgf = basis[k];
						if (cgf.l < 5) {
							for (int m = 0; m <= 2 * cgf.l; ++m) {
								++ao_index;

								fstream << (i64)ao_index << ',' << coefficient_matrix[1](k + molden_permutations[cgf.l][m], i) << '\n';
							}
						}

						k += 2 * cgf.l + 1;
					}
				}
			}
		}
	}

	struct WFXPrimitive {
		int kx = 0, ky = 0, kz = 0;
		double weight = 0, exponent = 0;
		int center = 0;
		int type = 1;

		WFXPrimitive(int kx, int ky, int kz, double exponent, int center, int type) : 
		kx(kx), ky(ky), kz(kz), weight(weight), exponent(exponent), center(center), type(type) {}

		WFXPrimitive(int l, int m, double exponent, int center) :
			weight(weight), exponent(exponent), center(center) {
			for (int j = 1; j <= l; ++j) {
				type += (j * j + j) / 2;
			}
			type += m;
			switch (l) {
			case 0:
				kx = 0;
				ky = 0;
				kz = 0;
				break;
			case 1:
				kx = m == 0;
				ky = m == 1;
				kz = m == 2;
				break;
			case 2: {
				const int x[6] = { 2, 0, 0, 1, 1, 0 };
				const int y[6] = { 0, 2, 0, 1, 0, 1 };
				const int z[6] = { 0, 0, 2, 0, 1, 1 };
				kx = x[m];
				ky = y[m];
				kz = z[m];
				break;
			}
			case 3: {
				const int x[10] = { 3, 0, 0, 2, 2, 0, 1, 1, 0, 1 };
				const int y[10] = { 0, 3, 0, 1, 0, 2, 2, 0, 1, 1 };
				const int z[10] = { 0, 0, 3, 0, 1, 1, 0, 2, 2, 1 };
				kx = x[m];
				ky = y[m];
				kz = z[m];
				break;
			}
			case 4: {
				const int x[15] = { 4, 0, 0, 3, 3, 1, 0, 1, 0, 2, 2, 0, 2, 1, 1 };
				const int y[15] = { 0, 4, 0, 1, 0, 3, 3, 0, 1, 2, 0, 2, 1, 2, 1 };
				const int z[15] = { 0, 0, 4, 0, 1, 0, 1, 3, 3, 0, 2, 2, 1, 1, 2 };
				kx = x[m];
				ky = y[m];
				kz = z[m];
				break;
			}
			}
		}
	};

	void writeOrbitalsWFX(const std::string& path, bool write_virtual) {
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 30);
		fout.centerTable(100);
		fout << '|' << '_' << '-' << "Writing orbitals to .wfx file" << '|' << '\n';
		fout.resetRows();
		fout << '\n';

		std::ofstream file(path, std::ios::binary);

		FormattedStream fstream;
		fstream.write_to = &file;

		file << "<Title>\n";
		file << "Dollar Store SCF output\n";
		file << "</Title>\n\n";

		file << "<Keywords>\n";
		file << "GTO\n";
		file << "</Keywords>\n\n";

		file << "<Number of Nuclei>\n";
		file << molecule.size() << '\n';
		file << "</Number of Nuclei>\n\n";

		file << "<Net Charge>\n";
		file << net_charge << '\n';
		file << "</Net Charge>\n\n";

		file << "<Number of Electrons>\n";
		file << (electron_count[0] + electron_count[1]) << '\n';
		file << "</Number of Electrons>\n\n";

		file << "<Number of Alpha Electrons>\n";
		file << electron_count[0] << '\n';
		file << "</Number of Alpha Electrons>\n\n";

		file << "<Number of Beta Electrons>\n";
		file << electron_count[1] << '\n';
		file << "</Number of Beta Electrons>\n\n";

		file << "<Electronic Spin Multiplicity>\n";
		file << (electron_count[0] - electron_count[1] + 1) << '\n';
		file << "</Electronic Spin Multiplicity>\n\n";

		file << "<Nuclear Names>\n";
		for (int i = 0; i < molecule.size(); ++i) {
			file << element_symbols[(int)molecule[i].element - 1] << (i + 1) << '\n';
		}
		file << "</Nuclear Names>\n\n";

		file << "<Atomic Numbers>\n";
		for (int i = 0; i < molecule.size(); ++i) {
			file << (int)molecule[i].element << '\n';
		}
		file << "</Atomic Numbers>\n\n";

		fstream.resetRows();
		fstream.addRow(NumberFormat::crudeFormat(10, 8), TextAlignment::right, 10);
		fstream.addRow(NumberFormat::crudeFormat(10, 8), TextAlignment::right, 10);
		fstream.addRow(NumberFormat::crudeFormat(10, 8), TextAlignment::right, 10);

		file << "<Nuclear Cartesian Coordinates>\n";
		for (int i = 0; i < molecule.size(); ++i) {
			fstream << molecule[i].position.x << ',' << molecule[i].position.y << ',' << molecule[i].position.z << '\n';
		}
		file << "</Nuclear Cartesian Coordinates>\n\n";

		std::vector<WFXPrimitive> primitives;

		for (int i = 0; i < molecule.size(); ++i) {
			for (int j = atom_basis[i]; j < atom_basis[i + 1]; j += 2 * basis[j].l + 1) {
				ContractedGaussian& cgf = basis[j];

				int l1 = cgf.l + 1;

				primitives.reserve(primitives.size() + (l1 * l1 + l1) / 2);

				if (cgf.l < 5) {
					for (int p = 0; p < cgf.primitives.size(); ++p) {
						for (int m = 0; m < (l1 * l1 + l1) / 2; ++m) {
							primitives.push_back(WFXPrimitive(cgf.l, m, cgf.primitives[p].zeta, i + 1));
						}
					}
				}
				else {
					int type = 1;
					for (int k = 1; k < l1; ++k) type += (k * k + k) / 2;

					for (int p = 0; p < cgf.primitives.size(); ++p) {
						for (int x = 0; x <= cgf.l; ++x) {
							for (int y = 0; y <= cgf.l - x; ++y) {
								primitives.push_back(WFXPrimitive(x, y, cgf.l - x - y, cgf.primitives[p].zeta, i + 1, type));
								++type;
							}
						}
					}
				}
			}
		}

		file << "<Number of Primitives>\n";
		file << primitives.size() << '\n';
		file << "</Number of Primitives>\n\n";

		file << "<Primitive Centers>\n";

		fstream.resetRows();
		for (int i = 0; i < 10; ++i) {
			fstream.addRow(NumberFormat(), TextAlignment::right, 4);
		}

		for (int i = 0; i < primitives.size(); ++i) {
			fstream << (i64)primitives[i].center << (i % 10 == 9 ? '\n' : ',');
		}
		if (primitives.size() % 10) fstream << '\n';
		file << "</Primitive Centers>\n\n";

		fstream.resetRows();
		for (int i = 0; i < 10; ++i) {
			fstream.addRow(NumberFormat(), TextAlignment::right, 4);
		}

		file << "<Primitive Types>\n";
		for (int i = 0; i < primitives.size(); ++i) {
			fstream << (i64)primitives[i].type << (i % 10 == 9 ? '\n' : ',');
		}
		if (primitives.size() % 10) fstream << '\n';
		file << "</Primitive Types>\n\n";

		fstream.resetRows();
		for (int i = 0; i < 5; ++i) {
			fstream.addRow(NumberFormat::scientificFormat(14, 8), TextAlignment::right, 14);
		}

		file << "<Primitive Exponents>\n";
		for (int i = 0; i < primitives.size(); ++i) {
			fstream << primitives[i].exponent << (i % 5 == 4 ? '\n' : ',');
		}
		if (primitives.size() % 5) fstream << '\n';
		file << "</Primitive Exponents>\n\n";

		fstream.resetRows();
		fstream.addRow(NumberFormat::scientificFormat(14, 8), TextAlignment::right, 14);

		int alpha_orbitals = (write_virtual ? matrix_size : electron_count[0]);
		int beta_orbitals = (write_virtual ? matrix_size : electron_count[1]) * (spin_treatment == SpinTreatment::unrestricted);

		file << "<Molecular Orbital Occupation Numbers>\n";
		for (int i = 0; i < alpha_orbitals; ++i) {
			fstream << (i64)((i < electron_count[0]) * (1 + (spin_treatment == SpinTreatment::restricted))) << '\n';
		}
		for (int i = 0; i < beta_orbitals; ++i) {
			fstream << (i64)(i < electron_count[0]) << '\n';
		}
		file << "</Molecular Orbital Occupation Numbers>\n\n";

		file << "<Molecular Orbital Energies>\n";
		for (int i = 0; i < alpha_orbitals; ++i) {
			fstream << mo_levels[0][i] << '\n';
		}
		for (int i = 0; i < beta_orbitals; ++i) {
			fstream << mo_levels[1][i] << '\n';
		}
		file << "</Molecular Orbital Energies>\n\n";

		file << "<Molecular Orbital Spin Types>\n";
		if (spin_treatment == SpinTreatment::unrestricted) {
			for (int i = 0; i < alpha_orbitals; ++i) {
				file << "Alpha\n";
			}
			for (int i = 0; i < beta_orbitals; ++i) {
				file << "Beta\n";
			}
		}
		else {
			for (int i = 0; i < alpha_orbitals; ++i) {
				file << "Alpha and Beta\n";
			}
		}
		file << "</Molecular Orbital Spin Types>\n\n";

		fstream.resetRows();
		for (int i = 0; i < 5; ++i) {
			fstream.addRow(NumberFormat::scientificFormat(14, 8), TextAlignment::right, 14);
		}

		file << "<Molecular Orbital Primitive Coefficients>\n";
		for (int i = 0; i < alpha_orbitals; ++i) {
			file << "<MO Number>\n";
			file << (i + 1) << '\n';
			file << "</MO Number>\n";

			for (int j = 0; j < primitives.size(); ++j) {
				primitives[j].weight = 0.0;
			}

			int primitive_index = 0;
			for (int j = 0; j < basis.size(); j += 2 * basis[j].l + 1) {
				int l1 = basis[j].l + 1;
				int p_count = (l1 * l1 + l1) / 2;

				for (int p = 0; p < basis[j].primitives.size(); ++p) {
					for (int k = 0; k < 2 * basis[j].l + 1; ++k) {
						for (int s = 0; s < basis[j + k].spherical_harmonic.size(); ++s) {
							CartesianPolynomialTerm& term = basis[j + k].spherical_harmonic[s];
							for (int t = primitive_index; t < primitive_index + p_count; ++t) {
								WFXPrimitive& prim = primitives[t];
								if (prim.kx == term.x && prim.ky == term.y && prim.kz == term.z) {
									prim.weight += coefficient_matrix[0](j + k, i) * term.weight * basis[j].primitives[p].weight;
									break;
								}
							}
						}
					}
					primitive_index += p_count;
				}
			}

			for (int i = 0; i < primitives.size(); ++i) {
				fstream << primitives[i].weight << (i % 5 == 4 ? '\n' : ',');
			}
			if (primitives.size() % 5) fstream << '\n';
		}

		for (int i = 0; i < beta_orbitals; ++i) {
			file << "<MO Number>\n";
			file << (i + 1 + alpha_orbitals) << '\n';
			file << "</MO Number>\n";

			for (int j = 0; j < primitives.size(); ++j) {
				primitives[j].weight = 0.0;
			}

			int primitive_index = 0;
			for (int j = 0; j < basis.size(); j += 2 * basis[j].l + 1) {
				int l1 = basis[j].l + 1;
				int p_count = (l1 * l1 + l1) / 2;

				for (int p = 0; p < basis[j].primitives.size(); ++p) {
					for (int k = 0; k < 2 * basis[j].l + 1; ++k) {
						for (int s = 0; s < basis[j + k].spherical_harmonic.size(); ++s) {
							CartesianPolynomialTerm& term = basis[j + k].spherical_harmonic[s];
							for (int t = primitive_index; t < primitive_index + p_count; ++t) {
								WFXPrimitive& prim = primitives[t];
								if (prim.kx == term.x && prim.ky == term.y && prim.kz == term.z) {
									prim.weight += coefficient_matrix[1](j + k, i) * term.weight * basis[j].primitives[p].weight;
									break;
								}
							}
						}
					}
					primitive_index += p_count;
				}
			}

			for (int i = 0; i < primitives.size(); ++i) {
				fstream << primitives[i].weight << (i % 5 == 4 ? '\n' : ',');
			}
			if (primitives.size() % 5) fstream << '\n';
		}
		file << "</Molecular Orbital Primitive Coefficients>\n\n";
	}
}