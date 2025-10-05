#pragma once
#include "../../util/Types.hpp"
#include "SCFSolver.hpp"

namespace scf {
	double getExactCoulombMatrix(uint mu, uint nu);

	double getExactExchangeMatrix(uint mu, uint nu, Spin spin);

	void assignExactRepulsionTensor();
}