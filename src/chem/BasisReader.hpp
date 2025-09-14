#pragma once
#include "GTO.hpp"
#include <string>

namespace flo {
namespace Gaussian {
	BasisSet readBasis(const std::string& path);
}
}