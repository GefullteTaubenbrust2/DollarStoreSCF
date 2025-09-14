#include "TextReader.hpp"
#include "TextUtil.hpp"
#include <iostream>
#include <vector>

namespace flo {
namespace TextReader {
	std::vector<std::string> file;
	std::string line;
	uint row = 0;
	uint column = 0;

	void openFile(const std::string& path) {
		file.clear();
		row = 0;
		column = 0;
		line = "";
		std::ifstream stream;
		stream.open(path);
		if (!stream) {
			std::cerr << "Warning: File '" << path << "' is empty or could not be found!\n";
			return;
		}
		std::string l;
		while (!stream.eof()) {
			safeGetline(stream, l);
			file.push_back(l);
		}
		if (file.size()) line = file[0];
	}

	std::string getLine() {
		return line;
	}

	void nextLine() {
		column = 0;
		++row;
		if (endOfFile()) line = "";
		else line = file[row];
	}

	void skipWhitespace() {
		for (; column < line.size(); ++column) {
			if (!isWhitespace(line[column])) return;
		}
	}

	bool endOfLine() {
		return column >= line.size();
	}

	bool endOfFile() {
		return row >= file.size();
	}

	double readFloat() {
		uint start_index = column;
		char c = safeGetChar(line, column);
		bool sign = c == '-' || c == '+' || c == '.';
		if (!sign && !isDigit(c)) return NAN;
		if (sign && !isDigit(safeGetChar(line, column + 1))) return NAN;
		++column;
		for (; column < line.size(); ++column) {
			c = line[column];
			bool exponent = c == 'E' || c == 'e' || c == 'D' || c == 'd';
			bool sign = c == '-' || c == '+';
			if (!isDigit(c) && c != '.' && !exponent && !sign) break;
			if (exponent) {
				char c1 = safeGetChar(line, column + 1);
				if (isDigit(c1)) continue;
				if (!(c1 == '-' || c1 == '+') || !isDigit(safeGetChar(line, column + 2))) break;
			}
			else if (sign && !isDigit(safeGetChar(line, column + 1))) break;
		}
		std::string substr = safeSubstr(line, start_index, column);
		for (int i = 0; i < substr.size(); ++i) {
			char c = substr[i];
			if (c == 'e' || c == 'D' || c == 'd') substr[i] = 'E';
		}
		return std::atof(substr.c_str());
	}

	int readInt() {
		uint start_index = column;
		char c = safeGetChar(line, column);
		bool sign = c == '-' || c == '+';
		if (!sign && !isDigit(c)) return NAN;
		if (sign && !isDigit(safeGetChar(line, column + 1))) return NAN;
		++column;
		for (; column < line.size(); ++column) {
			if (!isDigit(line[column])) break;
		}
		std::string substr = safeSubstr(line, start_index, column);
		return std::atoi(substr.c_str());
	}

	std::string readString() {
		uint start_index = column;
		uint end_index = line.size();
		char c = safeGetChar(line, column);
		if (c == '"') {
			++column;
			for (; column < line.size(); ++column) {
				if (line[column] == '"') {
					end_index = column;
					++column;
					break;
				}
			}
			return safeSubstr(line, start_index, end_index);
		}
		else if (!isWhitespace(c)) {
			for (; column < line.size(); ++column) {
				if (isWhitespace(line[column])) break;
			}
			return safeSubstr(line, start_index, column);

		}
		else return "";
	}

	char getChar() {
		return safeGetChar(line, column);
	}

	uint getLineNumber() {
		return row;
	}

	uint getCursorOffset() {
		return column;
	}

	void setLineNumber(uint line) {
		row = line;
	}

	void setCursorOffset(uint offset) {
		column = offset;
	}
}
}