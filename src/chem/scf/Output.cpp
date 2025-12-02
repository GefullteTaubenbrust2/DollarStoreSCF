#include "Output.hpp"
#include <iostream>

#include "SCFCommon.hpp"
#include "Energy.hpp"

#include "../../lalib/Lalib.hpp"
#include "../../lalib/PrintMatrix.h"

#include "../../util/FormattedStream.hpp"

using namespace flo;

namespace scf {
	void printTitle() {
        std::cout << "\
___________________________________________________________________________________________________\n\
    _____   _______ __      __      _______ _____           _______ _______ _______ _____   _______\n\
   / ___ \\ / ___  // /     / /     / ___  // __  \\         / _____//_   __// ___  // __  \\ / _____/\n\
  / /  / // /  / // /     / /     / /__/ // /_/ _/        /_/____   /  /  / /  / // /_/ _// /____\n\
 / /__/ // /__/ // /____ / /____ / ___  // ___ \\         _____/ /  /  /  / /__/ // ___ \\ / /____\n\
/______//______//______//______//_/  /_//_/  /_/        /______/  /__/  /______//_/  /_//______/\n\
___________________________________________________________________________________________________\n\
                            ____________    ____________   ______________\n\
                           /    ____    \\  /     __     \\ |              |\n\
                          /    /    \\____\\|     /  \\_____||     _________|\n\
                          \\    \\________  |    |          |    |_______\n\
                           \\________    \\ |    |     _____|     _______|\n\
                          _____     \\    \\|    |    |    ||    |\n\
                          \\    \\____/    /|     \\__/     ||    |\n\
                           \\____________/  \\____________/ |____|\n\
___________________________________________________________________________________________________\n\
\n\
                                       Gefuellte Taubenbrust\n\
___________________________________________________________________________________________________\n\
";
	}

	void printMOLevels() {
		if (spin_treatment == SpinTreatment::unrestricted) {
			fout.resetRows();
			fout.offsetRight(2);

			fout.addRow(NumberFormat(), TextAlignment::centered, 20 + 65);
			fout << '|' << '-' << '_' << "MO levels" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 20 + 25);
			fout.addRow(NumberFormat(), TextAlignment::centered, 20 + 15);
			fout << '|' << '-' << '_' << "Alpha" << '|' << '-' << '_' << "Beta" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 5);
			fout.addRow(NumberFormat(), TextAlignment::centered, 20);
			fout.addRow(NumberFormat(), TextAlignment::centered, 10);
			fout.addRow(NumberFormat(), TextAlignment::centered, 20);
			fout.addRow(NumberFormat(), TextAlignment::centered, 10);

			fout << '|' << '-' << '_' << "MO" << '|' << '-' << '_' << "Energy" << '|' << '-' << '_' << "Occupation" << '|' << '-' << '_' << "Energy" << '|' << '-' << '_' << "Occupation" << '|';

			fout.resetRows();
			fout.addRow(NumberFormat::crudeFormatPositive(5, 5), TextAlignment::right, 5);
			fout.addRow(NumberFormat::crudeFormatPositive(16, 12), TextAlignment::right, 20);
			fout.addRow(NumberFormat::crudeFormatPositive(8, 4), TextAlignment::right, 10);
			fout.addRow(NumberFormat::crudeFormatPositive(16, 12), TextAlignment::right, 20);
			fout.addRow(NumberFormat::crudeFormatPositive(8, 4), TextAlignment::right, 10);

			for (int i = 0; i < mo_levels[0].size(); ++i) {
				char underline_char = (i == electron_count[0] - 1 || i == electron_count[1] - 1) ? '_' : '\0';
				fout << '|' << underline_char << (i64)i << '|' << underline_char << mo_levels[0][i] << " Ha" << '|' << underline_char << (double)(i < electron_count[0]);
				fout << '|' << underline_char << mo_levels[1][i] << " Ha" << '|' << underline_char << (double)(i < electron_count[1]) << '|' << '\n';
			}
			fout << '-' << ',' << '-' << ',' << '-' << ',' << '-' << ',' << '-' << '\n';

			fout.resetRows();
			fout << '\n';
		}
		else {
			fout.resetRows();
			fout.offsetRight(2);

			fout.addRow(NumberFormat(), TextAlignment::centered, 20 + 25);
			fout << '|' << '-' << '_' << "MO levels" << '|';
			fout.resetRows();
			fout.addRow(NumberFormat(), TextAlignment::centered, 5);
			fout.addRow(NumberFormat(), TextAlignment::centered, 20);
			fout.addRow(NumberFormat(), TextAlignment::centered, 10);

			fout << '|' << '-' << '_' << "MO" << '|' << '-' << '_' << "Energy" << '|' << '-' << '_' << "Occupation" << '|';

			fout.resetRows();
			fout.addRow(NumberFormat::crudeFormatPositive(5, 5), TextAlignment::right, 5);
			fout.addRow(NumberFormat::crudeFormatPositive(16, 12), TextAlignment::right, 20);
			fout.addRow(NumberFormat::crudeFormatPositive(8, 4), TextAlignment::right, 10);

			for (int i = 0; i < mo_levels[0].size(); ++i) {
				char underline_char = (i == electron_count[0] - 1) ? '_' : '\0';
				fout << '|' << underline_char << (i64)i << '|' << underline_char << mo_levels[0][i] << " Ha" << '|' << underline_char << 2.0 * (double)(i < electron_count[0]) << '|' << '\n';
			}
			fout << '-' << ',' << '-' << ',' << '-' << '\n';

			fout.resetRows();
			fout << '\n';
		}
	}

	void printEnergyContributions() {
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::left, 20);
		fout.addRow(NumberFormat::crudeFormatPositive(16, 12), TextAlignment::centered, 20);
		fout.offsetRight(2);

		double kinetic_energy = getKineticEnergy();
		double potential_energy = getNuclearAttractionEnergy();
		double repulsion_energy = getElectronRepulsionEnergy();
		double exchange_energy = getExchangeEnergy();
		double correlation_energy = getCorrelationEnergy();

		fout << '|' << '-' << '_' << "Energy contribution" << '|' << '-' << '_' << "Value" << '|' << '\n';
		fout << '|' << "Kinetic" << '|' << fout.formatRow(TextAlignment::right) << kinetic_energy << " Ha" << '|' << '\n';
		fout << '|' << "Nuclear attraction" << '|' << potential_energy << " Ha" << '|' << '\r';
		fout << '|' << "Total 1e" << '|' << (kinetic_energy + potential_energy) << " Ha" << '|' << '\r';
		fout << '|' << "Electron repulsion" << '|' << repulsion_energy << " Ha" << '|' << '\n';
		fout << '|' << "Exchange" << '|' << exchange_energy << " Ha" << '|' << '\n';
		fout << '|' << "Correlation" << '|' << correlation_energy << " Ha" << '|' << '\r';
		fout << '|' << "Total 2e" << '|' << (repulsion_energy + exchange_energy + correlation_energy) << " Ha" << '|' << '\r';
		fout << '|' << "Total electronic" << '|' << (kinetic_energy + potential_energy + repulsion_energy + exchange_energy + correlation_energy) << " Ha" << '|' << '\r';
		fout << '|' << "Nuclear repulsion" << '|' << getNuclearRepulsionEnergy() << " Ha" << '|' << '\r';
		fout << '|' << '_' << "Total" << '|' << '_' << (kinetic_energy + potential_energy + repulsion_energy + exchange_energy + correlation_energy + getNuclearRepulsionEnergy()) << " Ha" << '|' << '\n' << '\n';

		fout.resetRows();
		fout << '\n';
	}

	void printMatrix(const MatrixBase<double>& matrix, const std::string title) {
		fout.resetRows();
		fout.offsetRight(2);
		fout.addRow(NumberFormat(), TextAlignment::centered, 87);
		fout << '|' << '_' << '-' << title << '|';
		fout.resetRows();
		fout.addRow(NumberFormat::crudeFormat(5, 5), TextAlignment::right, 5);
		fout.addRow(NumberFormat::scientificFormat(12, 5), TextAlignment::right, 12, TextPadding(2, 0, 0, 0));
		fout.addRow(NumberFormat::scientificFormat(12, 5), TextAlignment::right, 12, TextPadding(0, 0, 0, 0));
		fout.addRow(NumberFormat::scientificFormat(12, 5), TextAlignment::right, 12, TextPadding(0, 0, 0, 0));
		fout.addRow(NumberFormat::scientificFormat(12, 5), TextAlignment::right, 12, TextPadding(0, 0, 0, 0));
		fout.addRow(NumberFormat::scientificFormat(12, 5), TextAlignment::right, 12, TextPadding(0, 0, 0, 0));
		fout.addRow(NumberFormat::scientificFormat(12, 5), TextAlignment::right, 12, TextPadding(0, 2, 0, 0));
		fout << matrix;
	}
}