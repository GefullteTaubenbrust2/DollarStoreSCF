#include "PrintMatrix.h"

namespace flo {
	FormattedStream& operator<<(FormattedStream& stream, const MatrixBase<double>& matrix) {
		if (stream.getColumnCount() < 2) return stream;
		uint entries = stream.getColumnCount() - 1;
		for (int xoff = 0; xoff < matrix.getWidth(); xoff += entries) {
			stream << '|' << '_' << '-' << '|';
			for (int x = xoff; x < xoff + entries; ++x) {
				char seperator = x == xoff + entries - 1 ? '|' : ',';
				if (x < matrix.getWidth()) stream << '_' << '-' << (i64)x << seperator;
				else stream << '_' << '-' << "" << seperator;
			}
			stream << '\n';
			for (int y = 0; y < matrix.getHeight(); ++y) {
				stream << '|' << (i64)y << '|';
				for (int x = xoff; x < xoff + entries; ++x) {
					char seperator = x == xoff + entries - 1 ? '|' : ',';
					if (x < matrix.getWidth()) stream << matrix(y, x) << seperator;
					else stream << "" << seperator;
				}
				stream << '\n';
			}
		}
		for (int i = 0; i < entries + 1; ++i) {
			char seperator = i == entries ? '\n' : ',';
			stream << '-' << " " << seperator;
		}
		return stream;
	}
}