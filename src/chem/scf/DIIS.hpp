#pragma once
#include "../../lalib/Lalib.hpp"

namespace scf {
	void computeErrorVector(const flo::HermitianMatrixNd& fock_matrix);

	flo::HermitianMatrixNd& getDIISResult();
}