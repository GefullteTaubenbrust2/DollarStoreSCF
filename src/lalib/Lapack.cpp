#include "Lapack.hpp"
#include "Lalib.hpp"

namespace flo {
	extern "C" {
		extern int dspev_(char*, char*, int*, double*, double*, double*, int*, double*, int*);
		extern int dgesvd_(char*, char*, int*, int*, double*, int*, double*, double*, int*, double*, int*, double*, int*, int*);
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

	void computeEigenvalues(SymmetricMatrixNd& A, VectorNd& R, double* work_buffer) {
		R.resize(A.getWidth());
		Array<double> buffer;

		if (!work_buffer) {
			buffer.resize(A.getWidth() * 3);
			work_buffer = buffer.getPtr();
		}

		char jobz = 'N';
		char uplo = 'U';
		int N = A.getWidth();
		double* ap = A.getRawData();
		double* w = R.getRawData();
		double* z = nullptr;
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

	void singularValueDecomposition(MatrixNd& A, MatrixNd& U, VectorNd& sigma, MatrixNd& VT, double* work_buffer, int buffer_size) {
		Array<double> buffer;
		U.resize(A.getWidth(), A.getWidth());
		VT.resize(A.getHeight(), A.getHeight());

		if (!work_buffer) {
			char jobu = 'A';
			char jobvt = 'A';
			int M = 0;
			int N = 0;
			double ap = 0;
			int lda = 0;
			double sp = 0;
			double up = 0;
			int ldu = 0;
			double vtp = 0;
			int ldv = 0;
			double w = 0;
			int lwork = -1;
			int info = 0;
			dgesvd_(&jobu, &jobvt, &M, &N, &ap, &lda, &sp, &up, &ldu, &vtp, &ldv, &w, &lwork, &info);
			buffer_size = w;
			buffer.resize(buffer_size);
			work_buffer = buffer.getPtr();
		}

		char jobu = 'A';
		char jobvt = 'A';
		int M = A.getHeight();
		int N = A.getWidth();
		double* ap = A.getRawData();
		int lda = N;
		double* sp = sigma.getRawData();
		double* up = U.getRawData();
		int ldu = M;
		double* vtp = VT.getRawData();
		int ldv = N;
		int info = 0;

		dgesvd_(&jobu, &jobvt, &M, &N, ap, &lda, sp, up, &ldu, vtp, &ldv, work_buffer, &buffer_size, &info);

		if (info < 0) {
			std::cerr << "LAPACK error (general single value decomposition): illegal argument " << -info << '\n';
		}
		else if (info > 0) {
			std::cerr << "LAPACK error (general single value decomposition): off-diagonal element (" << info << ',' << info << ") is not 0\n";
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