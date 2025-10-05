#include "Lapack.hpp"
#include "Lalib.hpp"

namespace flo {
	extern "C" {
		extern int dspev_(char*, char*, int*, double*, double*, double*, int*, double*, int*);
		extern int dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
		extern int dspsv_(char*, int*, int*, double*, int*, double*, int*, int*);
	}

	void computeEigenvectors(SymmetricMatrixNd& A, MatrixNd& Q, VectorNd& R, double* work_buffer) {
		Q.resize(A.getWidth(), A.getWidth());
		R.resize(A.getWidth());
		Array<double> buffer;

		if (!work_buffer) {
			buffer.resize(A.getWidth() * 3);
			work_buffer = buffer.getPtr();
		}

		char jobz = 'V';
		char uplo = 'U';
		int N = A.getWidth();
		double* ap = A.getRawData();
		double* w = R.getRawData();
		double* z = Q.getRawData();
		int ldz = N;
		int info;

		dspev_(&jobz, &uplo, &N, ap, w, z, &ldz, work_buffer, &info);

		if (info < 0) {
			std::cerr << "LAPACK error (symmetric eigenvalue QR iteration): illegal argument " << -info << '\n';
		}
		else if (info > 0) {
			std::cerr << "LAPACK error (symmetric eigenvalue QR iteration): failure to converge\n";
		}
	}

	void solveLinear(MatrixNd& A, VectorNd& B, VectorNd& x, int* ipiv) {
		Array<int> buffer;
		if (!ipiv) {
			buffer.resize(A.getHeight());
			ipiv = buffer.getPtr();
		}

		if (&B != &x) x = B;

		int N = A.getHeight();
		int nhrs = 1;
		double* ap = A.getRawData();
		int lda = A.getHeight();
		double* bp = x.getRawData();
		int ldb = A.getHeight();
		int info = 0;

		dgesv_(&N, &nhrs, ap, &lda, ipiv, bp, &ldb, &info);

		if (info < 0) {
			std::cerr << "LAPACK error (general LU decomposition): illegal argument " << -info << '\n';
		}
		else if (info > 0) {
			std::cerr << "LAPACK error (general LU decomposition): U(" << info << ',' << info << ") is exactly 0\n";
		}
	}

	void solveLinear(SymmetricMatrixNd& A, VectorNd& B, VectorNd& x, int* ipiv) {
		Array<int> buffer;
		if (!ipiv) {
			buffer.resize(A.getHeight());
			ipiv = buffer.getPtr();
		}

		if (&B != &x) x = B;

		char uplo = 'U';
		int N = A.getHeight();
		int nhrs = 1;
		double* ap = A.getRawData();
		double* bp = x.getRawData();
		int ldb = A.getHeight();
		int info = 0;

		dspsv_(&uplo, &N, &nhrs, ap, ipiv, bp, &ldb, &info);

		if (info < 0) {
			std::cerr << "LAPACK error (symmetric LDL decomposition): illegal argument " << -info << '\n';
		}
		else if (info > 0) {
			std::cerr << "LAPACK error (symmetric LDL decomposition): D(" << info << ',' << info << ") is exactly 0\n";
		}
	}
}