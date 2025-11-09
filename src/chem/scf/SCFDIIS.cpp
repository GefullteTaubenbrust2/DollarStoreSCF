#include "SCFDIIS.hpp"
#include "SCFCommon.hpp"
#include "../../lalib/Lapack.hpp"
#include "../../lalib/DIIS.ipp"

using namespace flo;

namespace scf {
namespace diis {
	typedef DIISSolver<double, MatrixNd, SymmetricMatrixNd> SCFDIIS;

	SCFDIIS solvers[2];

	void setIterationCount(uint iterations) {
		solvers[0].resize(iterations);
		if (spin_treatment == SpinTreatment::unrestricted) solvers[1].resize(iterations);
	}

	void computeErrorVector(flo::MatrixNd& error_vector, const flo::SymmetricMatrixNd& fock_matrix, Spin spin) {
		MatrixNd& spf = asymmetric_buffer[1];

		spf = (asymmetric_buffer[2] = overlap_matrix * (asymmetric_buffer[3] = density_matrix[(int)spin] * fock_matrix));

		error_vector = spf - (asymmetric_buffer[2] = T(spf));
	}

	void addPulayErrorVector(const flo::SymmetricMatrixNd& fock_matrix, Spin spin) {
		MatrixNd& error_vector = asymmetric_buffer[0];
		
		computeErrorVector(error_vector, fock_matrix, spin);

		solvers[(int)spin].addErrorVector(fock_matrix, error_vector);
	}

	double getRMS() {
		double rms = solvers[0].getRMS();
		if (spin_treatment == SpinTreatment::unrestricted) {
			rms = 0.5 * rms + solvers[1].getRMS();
		}
		return rms;
	}

	flo::SymmetricMatrixNd& getResult(Spin spin) {
		solvers[(int)spin].calculateResult(buffer[0]);
		return buffer[0];
	}

	void solve() {
		solvers[0].solve(&buffer[0]);
		if (spin_treatment == SpinTreatment::unrestricted) solvers[1].solve(&buffer[0]);
		buffer[0].resize(matrix_size);
	}
}
}