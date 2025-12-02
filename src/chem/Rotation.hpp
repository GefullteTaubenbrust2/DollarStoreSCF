#pragma once
#include "../lalib/MathUtil.hpp"
#include "Molecule.hpp"

namespace flo {
	vec3 calculateMomentOfInertia(const Molecule& mol, MatrixNd& vectors);
}