#pragma once

namespace scf {
	void calculateEnergies();

	double getKineticEnergy();

	double getNuclearAttractionEnergy();

	double getElectronRepulsionEnergy();

	double getExchangeEnergy();

	double getCorrelationEnergy();

	double getNuclearRepulsionEnergy();

	double getTotalElectronicEnergy();

	double getTotalEnergy();
}