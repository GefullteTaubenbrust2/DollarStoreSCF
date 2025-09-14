#include "DIIS.hpp"

using namespace flo;

namespace scf {
	HermitianMatrixNd placeholder;

	void computeErrorVector(const HermitianMatrixNd& fock_matrix) {

	}

	HermitianMatrixNd& getDIISResult() {
		return placeholder;
	}
}