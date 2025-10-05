#pragma once
#include "../../lalib/Lalib.hpp"

namespace scf {
	extern flo::DiagonalMatrixNd mulliken_populations[2];
	extern flo::DiagonalMatrixNd lowdin_populations[2];

	void calculateMullikenPopulations();

	void printMullikenPopulations();

	void printMullikenBondOrders();

	void calculateLowdinPopulations();

	void printLowdinPopulations();

	void printWibergBondOrders();

	void printMayerValence();

	void printMayerBondOrders();
}