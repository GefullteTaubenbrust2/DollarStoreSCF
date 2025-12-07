#include "Thermochem.hpp"

#include "Constants.hpp"
#include "PointGroup.hpp"
#include "Rotation.hpp"

#include "../lalib/Lalib.hpp"

#include "../util/FormattedStream.hpp"

using namespace flo;

namespace scf {
	extern VectorNd frequencies;
	extern Molecule molecule;
	namespace freq {
		extern double equilibrium_total_energy;
	}
}

namespace flo {
	struct CanonicalPartitionFunction {
		double lnZ = 0.0;
		double dlnZ_dbeta = 0.0;
		double dlnZ_dV = 0.0;
		double temperature = 298.15;
		double volume = 1.0;

		CanonicalPartitionFunction() = default;

		CanonicalPartitionFunction(double temperature, double volume) : temperature(temperature), volume(volume) {};

		double getInternalEnergy() {
			return -dlnZ_dbeta;
		}

		double getEntropy() {
			return -1.0 / temperature * dlnZ_dbeta + kB_in_Ha_K * lnZ;
		}

		double getFreeEnergy() {
			return -kB_in_Ha_K * temperature * lnZ;
		}

		double getPressure() {
			return kB_in_Ha_K * temperature * dlnZ_dV;
		}

		double getEnthalpy() {
			return -dlnZ_dbeta + kB_in_Ha_K * temperature * volume * dlnZ_dV;
		}

		double getFreeEnthalpy() {
			return -kB_in_Ha_K * temperature * lnZ + kB_in_Ha_K * temperature * volume * dlnZ_dV;
		}

		CanonicalPartitionFunction operator+(const CanonicalPartitionFunction& other) {
			CanonicalPartitionFunction result = *this;
			result.lnZ += other.lnZ;
			result.dlnZ_dbeta += other.dlnZ_dbeta;
			result.dlnZ_dV += other.dlnZ_dV;
			return result;
		}
	};

	CanonicalPartitionFunction getElectronicPartitionFunction(double energy, double temperature, double volume) {
		CanonicalPartitionFunction result(temperature, volume);

		double beta = 1.0 / (kB_in_Ha_K * temperature);

		result.lnZ = -beta * energy;
		result.dlnZ_dbeta = -energy;

		return result;
	}

	CanonicalPartitionFunction getZeroPointFunction(VectorNd& frequencies, double temperature, double volume) {
		CanonicalPartitionFunction result(temperature, volume);

		double beta = 1.0 / (kB_in_Ha_K * temperature);

		for (int i = 0; i < frequencies.size(); ++i) {
			if (frequencies[i] < 1.0 / au_to_per_cm) continue;
			result.lnZ -= 0.5 * beta * frequencies[i];
			result.dlnZ_dbeta -= 0.5 * frequencies[i];
		}

		return result;
	}

	CanonicalPartitionFunction getVibrationalPartitionFunction(VectorNd& frequencies, double temperature, double volume) {
		CanonicalPartitionFunction result(temperature, volume);

		double beta = 1.0 / (kB_in_Ha_K * temperature);

		for (int i = 0; i < frequencies.size(); ++i) {
			if (frequencies[i] < 1.0 / au_to_per_cm) continue;
			double exponential_term = std::exp(-beta * frequencies[i]);
			result.lnZ -= std::log(1.0 - exponential_term);
			result.dlnZ_dbeta -= frequencies[i] * exponential_term / (1.0 - exponential_term);
		}

		return result;
	}

	CanonicalPartitionFunction getClassicRotatorPartitionFunction(vec3 moments_of_inertia, bool linear, double sigma, double temperature, double volume) {
		CanonicalPartitionFunction result(temperature, volume);

		double beta = 1.0 / (kB_in_Ha_K * temperature);

		if (linear) {
			result.lnZ = std::log(2.0 * moments_of_inertia[2] / (beta * sigma));
			result.dlnZ_dbeta = -kB_in_Ha_K * temperature;
		}
		else {
			result.lnZ = 0.5 * std::log(moments_of_inertia[0] * moments_of_inertia[1] * moments_of_inertia[2]) + 1.5 * std::log(2.0 * PI) - 1.5 * std::log(beta) - std::log(PI * sigma);
			result.dlnZ_dbeta = -1.5 * kB_in_Ha_K * temperature;
		}

		return result;
	}

	CanonicalPartitionFunction getTranslationPartitionFunction(double mass, double temperature, double volume) {
		CanonicalPartitionFunction result(temperature, volume);

		double beta = 1.0 / (kB_in_Ha_K * temperature);

		result.lnZ = std::log(volume) - 1.5 * std::log(2.0 * PI) + 1.5 * std::log(mass / beta) + 1.0;
		result.dlnZ_dbeta = -1.5 * kB_in_Ha_K * temperature;
		result.dlnZ_dV = 1.0 / volume;

		return result;
	}

#define PRINT_THERMO(TITLE, QUANTITY, FACTOR)\
	fout.resetRows();\
	fout.addRow(NumberFormat::crudeFormatPositive(10, 8), TextAlignment::centered, 62);\
	fout << '|' << '_' << TITLE << '|' << '\n';\
	fout.resetRows();\
	fout.addRow(NumberFormat(), TextAlignment::centered, 20);\
	fout.addRow(NumberFormat(), TextAlignment::centered, 16);\
	fout.addRow(NumberFormat(), TextAlignment::centered, 16);\
	fout << '|' << '_' << "Term" << '|' << '_' << "Value [Ha]" << '|' << '_' << "Value [kJ/mol]" << '|' << '\n';\
	fout.resetRows();\
	fout.addRow(NumberFormat(), TextAlignment::left, 20);\
	fout.addRow(NumberFormat::crudeFormat(16, 12), TextAlignment::right, 16);\
	fout.addRow(NumberFormat::crudeFormat(16, 12), TextAlignment::right, 16);\
	fout << '|' << "Translation" << '|' << (translation_Z.QUANTITY() * FACTOR) << '|' << (translation_Z.QUANTITY() * FACTOR * Ha_to_kJmol) << '|' << '\n';\
	fout << '|' << "Rotation" << '|' << (rotation_Z.QUANTITY() * FACTOR) << '|' << (rotation_Z.QUANTITY() * FACTOR * Ha_to_kJmol) << '|' << '\n';\
	fout << '|' << "Vibration" << '|' << (vibration_Z.QUANTITY() * FACTOR) << '|' << (vibration_Z.QUANTITY() * FACTOR * Ha_to_kJmol) << '|' << '\n';\
	fout << '|' << "Zero point energy" << '|' << (zero_point_Z.QUANTITY() * FACTOR) << '|' << (zero_point_Z.QUANTITY() * FACTOR * Ha_to_kJmol) << '|' << '\n';\
	fout << '|' << '_' << "Electronic" << '|' << '_' << (electronic_Z.QUANTITY() * FACTOR) << '|' << '_' << (electronic_Z.QUANTITY() * FACTOR * Ha_to_kJmol) << '|' << '\n';\
	fout << '|' << "Total correction" << '|' << (total_Z_minus_e.QUANTITY() * FACTOR) << '|' << (total_Z_minus_e.QUANTITY() * FACTOR * Ha_to_kJmol) << '|' << '\n';\
	fout << '|' << '_' << "Total" << '|' << '_' << (total_Z.QUANTITY() * FACTOR) << '|' << '_' << (total_Z.QUANTITY() * FACTOR * Ha_to_kJmol) << '|' << '\n';

	void doThermochem(double temperature, double pressure_Pa) {
		double pressure = pressure_Pa * Pa_to_au;
		double volume = kB_in_Ha_K * temperature / pressure;

		double mass = 0.0;
		for (int i = 0; i < scf::molecule.size(); ++i) {
			mass += scf::molecule[i].mass;
		}

		//scf::frequencies = VectorNd({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1770 / au_to_per_cm, 4150 / au_to_per_cm, 4270 / au_to_per_cm});
		//scf::freq::equilibrium_total_energy = -76.0231255;

		MatrixNd inertia_axes(3, 3);
		PointGroup point_group = findPointGroup(scf::molecule, 0.01);
		vec3 moments_of_inertia = calculateMomentOfInertia(scf::molecule, inertia_axes);

		CanonicalPartitionFunction translation_Z = getTranslationPartitionFunction(mass, temperature, volume);

		CanonicalPartitionFunction rotation_Z = getClassicRotatorPartitionFunction(moments_of_inertia, !point_group.rotation, point_group.rotation ? point_group.rotation : 1 + point_group.inversion, temperature, volume);

		CanonicalPartitionFunction zero_point_Z = getZeroPointFunction(scf::frequencies, temperature, volume);

		CanonicalPartitionFunction vibration_Z = getVibrationalPartitionFunction(scf::frequencies, temperature, volume);

		CanonicalPartitionFunction electronic_Z = getElectronicPartitionFunction(scf::freq::equilibrium_total_energy, temperature, volume);

		CanonicalPartitionFunction total_Z_minus_e = translation_Z + rotation_Z + zero_point_Z + vibration_Z;

		CanonicalPartitionFunction total_Z = translation_Z + rotation_Z + zero_point_Z + vibration_Z + electronic_Z;

		fout.resetRows();
		fout.offsetRight(2);
		fout.addRow(NumberFormat::crudeFormatPositive(10, -2), TextAlignment::centered, 62);
		fout << '|' << '-' << '_' << "Thermochemistry" << '|' << '\n';
		fout << '|' << '_' << "Temperature: " << temperature << " K Pressure: " << (pressure_Pa / 100000.0) << " bar" << '|' << '\n';
		fout.resetRows();
		fout.addRow(NumberFormat::crudeFormatPositive(10, 8), TextAlignment::block, 62);
		fout << '|' << '_' << "The calculation makes the following assumptions:\n\n";
		fout << "1.\tNo excited electronic states, including possible degenerate ground states, are accessible.\n";
		fout << "2.\tAll vibrations are strictly harmonic. This includes modes with very small frequencies that should likely be treated as hindered rotations instead.\n";
		fout << "3.\tThe rotations are treated within the semiclassical model of the rigid rotator, neglecting quantum effects relevant at low temperature.\n";
		fout << "4.\tThe molecule is treated as an ideal gas.\n";
		fout << "5.\tThe point group of the molecule is " << point_group.getSchoenfliesSymbol() << "." << '|' << '\n';

		PRINT_THERMO("Inner Energy", getInternalEnergy, 1.0);
		PRINT_THERMO("Enthalpy", getEnthalpy, 1.0);
		PRINT_THERMO("Entropy, TS", getEntropy, temperature);
		PRINT_THERMO("Gibbs Free Energy", getFreeEnthalpy, 1.0);
		PRINT_THERMO("Helmholtz Free Energy", getFreeEnergy, 1.0);
		fout.resetRows();
		fout << '\n';
	}
}