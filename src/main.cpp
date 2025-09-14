#include <iostream>
#include "lalib/Lalib.hpp"
#include "lalib/Lapack.hpp"
#include "chem/GaussianIntegral.hpp"
#include "chem/BasisReader.hpp"
#include "chem/scf/SCFBasis.hpp"
#include "chem/scf/SCFSolver.hpp"
#include "chem/scf/ExactCoulomb.hpp"
#include "chem/XYZReader.hpp"
#include "chem/scf/Population.hpp"
#include "util/Time.hpp"
#include "util/FormattedStream.hpp"

using namespace flo;

int main() {
	/*SymmetricMatrixNd H(6, {-1.0, -0.5, -1.0, 0.0, -0.5, -1.0, 0.0, 0.0, -0.5, -1.0, 0.0, 0.0, 0.0, -0.5, -1.0, -0.5, 0.0, 0.0, 0.0, -0.5, -1.0});
	MatrixNd Q;
	VectorNd R;

	std::cout << H << '\n';

	computeEigenvectors(H, Q, R);

	std::cout << R << '\n';
	std::cout << Q << '\n';

	MatrixNd m(3, 3, { 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, 1.0, 2.0, 0.0 });
	VectorNd x({ 1.0, 0.0, 10.0 });
	VectorNd y({1.0, 1.0, 0.0});
	VectorNd z;
	MatrixNd n;
	MatrixNd o;

	std::cout << (z = (o = T(m)) * y) << '\n';*/

	//ContractedGaussian gto1({ GTOPrimitive(2.0, 1.0) }, 1, 0, vec3(0.0, 0.0, 0.0));
	//ContractedGaussian gto2({ GTOPrimitive(1.0, 1.0) }, 1, 0, vec3(1.0, 0.0, 0.0));

	//std::cout << flo::overlapIntegral(gto1, gto2) << '\n';
	//std::cout << flo::kineticEnergyIntegral(gto1, gto2) << '\n';
	//std::cout << flo::nuclearPotentialIntegral(gto1, gto2, flo::vec3(0.0, 0.0, 0.0)) << '\n';
	//std::cout << electronRepulsionIntegral(gto1, gto1, gto2, gto2) << '\n';
	//std::cout << electronRepulsionIntegral(gto1, gto2, gto2, gto1) << '\n';

	/*fout << '|' << '-' << "A1" << '|' << "B1" << '-' << '_' << ',' << "C1" << '|' << "D1" << '|' << '\n';
	fout << '|' << "A2" << '|' << "B2" << ',' << "C2" << ',' << "D2" << '|' << '\n';
	fout << '-' << "A3" << '_' << '|' << "B3" << ',' << "C3" << '-' << ',' << "D3" << '_' << '\n';
	fout << '|' << '_' << "A4" << '|' << "B4" << ',' << "C4" << '-' << ',' << "D4" << '\n' << '\n';*/

	/*for (int m = -3; m <= 3; ++m) {
		ContractedGaussian cgf({ GTOPrimitive(1.0, 1.0) }, 3, m);
		for (int i = 0; i < cgf.spherical_harmonic.size(); ++i) {
			std::cout << cgf.spherical_harmonic[i].x << ' ' << cgf.spherical_harmonic[i].y << ' ' << cgf.spherical_harmonic[i].z << "   ";
		}
		std::cout << '\n';
	}*/

	BasisSet basis = flo::Gaussian::readBasis("basis/6-31GPP");
	Molecule mol = readXYZFile("test/water.xyz");

	scf::assignMolecule(mol, 0, 1);
	scf::assignBasis(basis);
	scf::constructBasis();
	scf::solveMOs();
	scf::printMullikenPopulations();
	scf::printLowdinPopulations();

	system("pause");
	return 0;
}