#include "XYZReader.hpp"
#include "../util/TextReader.hpp"
#include "Constants.hpp"

namespace flo {
	Molecule readXYZFile(const std::string& path) {
		Molecule result;
		TextReader::openFile(path);
		TextReader::skipWhitespace();
		int size = TextReader::readInt();
		result.atoms.reserve(size);
		TextReader::nextLine();
		TextReader::nextLine();
		for (int i = 0; i < size && !TextReader::endOfFile(); ++i) {
			TextReader::skipWhitespace();
			std::string element_symbol = TextReader::readString();
			int atomic_number = getElement(element_symbol);
			if (!atomic_number) {
				std::cerr << "ERROR: '" << element_symbol << "' is not a valid element symbol. Check capitalization (e.g. Rn instead of RN)!\n";
				continue;
			}
			TextReader::skipWhitespace();
			double x = TextReader::readFloat();
			TextReader::skipWhitespace();
			double y = TextReader::readFloat();
			TextReader::skipWhitespace();
			double z = TextReader::readFloat();
			TextReader::nextLine();
			result.atoms.push_back(Atom(vec3(x, y, z) * A_to_a0, (Element)atomic_number));
		}
		if (result.atoms.size() < size) {
			std::cerr << "ERROR: " << size << " atoms expected in XYZ file, but only " << result.atoms.size() << " were specified. Check for misplaced line breaks!\n";
		}
		return result;
	}
}