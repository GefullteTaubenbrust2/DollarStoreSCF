#pragma once
#include "../../lalib/Lalib.hpp"
#include "SCFSolver.hpp"

namespace scf {
namespace diis {
	void setIterationCount(uint iterations);

	void addPulayErrorVector(const flo::SymmetricMatrixNd& fock_matrix, Spin spin);

	double getRMS();

	void solve();

	flo::SymmetricMatrixNd& getResult(Spin spin);
}
}