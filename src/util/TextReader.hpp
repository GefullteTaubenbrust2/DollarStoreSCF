#pragma once

#include <string>
#include "Types.hpp"

namespace flo {
namespace TextReader {
	void openFile(const std::string& path);

	std::string getLine();

	void nextLine();

	void skipWhitespace();

	bool endOfLine();

	bool endOfFile();

	double readFloat();

	int readInt();

	std::string readString();

	char getChar();

	uint getLineNumber();

	uint getCursorOffset();

	void setLineNumber(uint line);

	void setCursorOffset(uint offset);
}
}