#pragma once
#include "Types.hpp"
#include <iostream>
#include <ostream>
#include <vector>

namespace flo {
	struct FormattedStream;

	extern FormattedStream fout;

	enum class TextAlignment {
		left = 0,
		right = 1,
		centered = 2,
		block = 3,
	};

	struct NumberFormat {
		uint max_characters = 16;
		int significant_digits = 8;
		bool adjust_positives = true;
		bool truncate_zeros = false;
		bool scientific_notation = true;
		uint max_mantissa = 4;

		static NumberFormat scientificFormat(uint max_characters, int significant_digits, bool adjust_positives = true, bool truncate_zeros = false);

		static NumberFormat scientificFormatPositive(uint max_characters, int significant_digits, bool truncate_zeros = false);

		static NumberFormat crudeFormat(uint max_characters, int significant_digits, bool adjust_positives = true);

		static NumberFormat crudeFormatPositive(uint max_characters, int significant_digits);
	};

	std::string toFormattedString(double x, const NumberFormat& format);

	std::string toFormattedString(i64 x, const NumberFormat& format);

	struct FormattedStream {
	private:
		struct Cell {
			NumberFormat number_format;
			TextAlignment alignment = TextAlignment::left;
			uint width = 16;
			bool separator_before = false;
			bool hline_top = false;
			bool hline_bottom = false;

			Cell(NumberFormat number_format, TextAlignment alignment, uint width) :
			number_format(number_format), alignment(alignment), width(width) {}
		};

		std::vector<Cell> cells;
		std::vector<std::string> content;
		std::string hline_string;
		uint current_cell = 0;

		Cell& getCurrentCell();

		bool row_has_content = false;

		uint left_offset = 0;

		bool finish_hline = false;
		bool start_hline = false;

		uint column_count = 0;

	public:
		uint horizontal_padding = 2;
		uint vertical_padding = 0;

		char vertical_line = '|';
		char horizontal_line = '-';
		char line_intersection = '+';

		std::ostream* write_to = &std::cout;

		FormattedStream() = default;

		void flush();

		void resetRows();

		void addRow(const NumberFormat& number_format, TextAlignment alignment, uint width);

		char formatRow(const NumberFormat& number_format);

		char formatRow(TextAlignment alignment);

		uint getColumnCount();

		//char drawHline(char c);

		void offsetRight(uint offset);

		void centerTable(uint total_width);

		void nextEntry();

		FormattedStream& operator<<(const std::string str);

		FormattedStream& operator<<(char c);

		FormattedStream& operator<<(double number);

		FormattedStream& operator<<(i64 number);
	};
}