#include "TextUtil.hpp"

namespace flo {
	std::string safeSubstr(const std::string& str, uint start, uint end) {
		if (start > str.size()) start = str.size();
		if (end >= str.size()) end = str.size();
		if (start >= end) return "";
		return str.substr(start, end - start);
	}

	char safeGetChar(const std::string& str, uint index) {
		return index >= str.size() ? 0 : str[index];
	}

	std::istream& safeGetline(std::istream& is, std::string& t) {
		t.clear();

		// The characters in the stream are read one-by-one using a std::streambuf.
		// That is faster than reading them one-by-one using the std::istream.
		// Code that uses streambuf this way must be guarded by a sentry object.
		// The sentry object performs various tasks,
		// such as thread synchronization and updating the stream state.

		std::istream::sentry se(is, true);
		std::streambuf* sb = is.rdbuf();

		for (;;) {
			int c = sb->sbumpc();
			switch (c) {
			case '\n':
				return is;
			case '\r':
				if (sb->sgetc() == '\n')
					sb->sbumpc();
				return is;
			case std::streambuf::traits_type::eof():
				// Also handle the case when the last line has no line ending
				if (t.empty())
					is.setstate(std::ios::eofbit);
				return is;
			default:
				t += (char)c;
			}
		}
	}

	bool isWhitespace(char c) {
		return c == ' ' || c == '\t' || c == '\n' || c == '\r';
	}

	bool isLetter(char c) {
		return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
	}

	bool isDigit(char c) {
		return c >= '0' && c <= '9';
	}
}