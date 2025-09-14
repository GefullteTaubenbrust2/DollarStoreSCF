#pragma once
#include <string>
#include <fstream>
#include "Types.hpp"

namespace flo {
	std::string safeSubstr(const std::string& str, uint start, uint end);

	char safeGetChar(const std::string& str, uint index);

	std::istream& safeGetline(std::istream& is, std::string& t);

	bool isWhitespace(char c);

	bool isLetter(char c);

	bool isDigit(char c);
}