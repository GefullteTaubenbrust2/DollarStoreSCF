#include "BasisReader.hpp"
#include "Element.hpp"
#include "../util/TextReader.hpp"
#include "../util/TextUtil.hpp"
#include <iostream>

namespace flo {
namespace Gaussian {
	void readAtomicBasis(AtomicBasis& basis) {
		TextReader::skipWhitespace();
		std::string first_entry = TextReader::readString();

		if (safeSubstr(first_entry, first_entry.size() - 3, first_entry.size()) == "ECP") return;
		basis.clear();
		TextReader::setCursorOffset(0);

		for (;;) {
			GTOShell shell[2];
			shell[1].l = 1;

			TextReader::skipWhitespace();
			std::string type = TextReader::readString();

			int angular_momentum = 0;
			if (safeGetChar(type, 0) == '*') return;
			else if (type == "SP") angular_momentum = -1;
			else if (type.size() != 1) {
				std::cerr << "ERROR: " << type[0] << " is not a recognized angular momentum abbreviation. Accepted values are S, P, D etc. as well as SP.\n";
				return;
			}
			else switch (type[0]) {
			case 'S':
				angular_momentum = 0;
				break;
			case 'P':
				angular_momentum = 1;
				break;
			case 'D':
				angular_momentum = 2;
				break;
			case 'F':
				angular_momentum = 3;
				break;
			case 'G':
				angular_momentum = 4;
				break;
			case 'H':
				angular_momentum = 5;
				break;
			case 'I':
				angular_momentum = 6;
				break;
			case 'J':
				angular_momentum = 7;
				break;
			default:
				std::cerr << "ERROR: " << type[0] << " is not a recognized angular momentum abbreviation. Accepted values are S, P, D etc. as well as SP.\n";
				return;
			}
			if (angular_momentum >= 0) shell[0].l = angular_momentum;
			else shell[0].l = 0;

			TextReader::skipWhitespace();
			int primitive_count = TextReader::readInt();
			if (primitive_count < 1 || primitive_count > 100) {
				std::cerr << "ERROR: Specified number of primitives (" << primitive_count << ") is outside of the accepted range of 1-100.\n";
				return;
			}

			for (int i = 0; i < primitive_count; ++i) {
				TextReader::nextLine();
				TextReader::skipWhitespace();
				double zeta = TextReader::readFloat();
				if (zeta <= 0.0) {
					std::cerr << "ERROR: Primitive GTOs must have positive exponents!\n";
					return;
				}

				TextReader::skipWhitespace();
				double weight = TextReader::readFloat();

				shell[0].primitives.push_back(GTOPrimitive(zeta, weight));

				if (angular_momentum < 0) {
					TextReader::skipWhitespace();
					double weight2 = TextReader::readFloat();
					shell[1].primitives.push_back(GTOPrimitive(zeta, weight2));
				}
			}

			TextReader::nextLine();

			basis.push_back(shell[0]);
			if (angular_momentum < 0) basis.push_back(shell[1]);
		}
	}

	BasisSet readBasis(const std::string& path) {
		BasisSet result;
		TextReader::openFile(path);
		TextReader::skipWhitespace();
		while (!TextReader::endOfFile()) {
			std::string element_name = TextReader::readString();
			int element = getElement(element_name);
			TextReader::skipWhitespace();
			if (element && TextReader::getChar() == '0') {
				TextReader::nextLine();
				readAtomicBasis(result.atomic_bases[element - 1]);
			}
			TextReader::nextLine();
		}
		return result;
	}
}
}