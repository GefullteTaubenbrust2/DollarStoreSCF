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
#include "chem/scf/Title.hpp"
#include "chem/scf/WriteOrbitals.hpp"

using namespace flo;

int main() {
	/*ContractedGaussian a({GTOPrimitive(1.292, 1.0)}, 2, -2);
	ContractedGaussian b({ GTOPrimitive(0.75, 1.0) }, 1, 1, vec3(-1.4316565, 1.1092692, 0.0));

	std::cout << kineticEnergyIntegral(a, b) << '\n';*/

	scf::printTitle();

	BasisSet basis = flo::Gaussian::readBasis("basis/6-31GPP");
	Molecule mol = readXYZFile("test/water.xyz");

	scf::assignMolecule(mol, 0, 1);
	scf::assignBasis(basis);
	scf::constructBasis();
	scf::solveMOs();
	scf::printMullikenPopulations();
	scf::printLowdinPopulations();
	scf::printMullikenBondOrders();
	scf::printWibergBondOrders();
	scf::printMayerBondOrders();
	scf::printMayerValence();
	scf::writeOrbitalsMolden("test/water.molden");
	scf::writeOrbitalsWFX("test/water.wfx", true);

	system("pause");
	return 0;
}