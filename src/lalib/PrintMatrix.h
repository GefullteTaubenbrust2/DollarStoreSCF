#pragma once
#include "../util/FormattedStream.hpp"
#include "Lalib.hpp"

namespace flo {
	FormattedStream& operator<<(FormattedStream& stream, const MatrixBase<double>& matrix);
}