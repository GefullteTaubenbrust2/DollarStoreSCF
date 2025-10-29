#pragma once
#include "../../util/Types.hpp"
#include "SCFSolver.hpp"
#include "../../lalib/Lalib.hpp"

using namespace flo;

namespace scf {
	double getExactCoulombMatrix(uint mu, uint nu);

	double getExactExchangeMatrix(uint mu, uint nu, Spin spin);

	void assignExactRepulsionTensor();

	void addExactTwoElectronGradient(VectorNd& target);
}