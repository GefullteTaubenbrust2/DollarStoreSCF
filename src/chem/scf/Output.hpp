#pragma once
#include <string>
#include "../../lalib/Lalib.hpp"

namespace scf {
	void printTitle();

	void printMOLevels();

	void printEnergyContributions();

	void printMatrix(const flo::MatrixBase<double>& matrix, const std::string title);
}