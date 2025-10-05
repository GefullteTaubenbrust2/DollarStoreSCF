#include "FormattedStream.hpp"
#include "../lalib/MathUtil.hpp"
#include "TextUtil.hpp"

namespace flo {
	FormattedStream fout;

	NumberFormat NumberFormat::scientificFormat(uint max_characters, int significant_digits, bool adjust_positives, bool truncate_zeros) {
		NumberFormat result;
		result.max_characters = max_characters;
		result.significant_digits = significant_digits;
		result.adjust_positives = adjust_positives;
		result.truncate_zeros = truncate_zeros;
		return result;
	}

	NumberFormat NumberFormat::scientificFormatPositive(uint max_characters, int significant_digits, bool truncate_zeros) {
		NumberFormat result;
		result.max_characters = max_characters;
		result.significant_digits = significant_digits;
		result.adjust_positives = false;
		result.truncate_zeros = truncate_zeros;
		return result;
	}

	NumberFormat NumberFormat::crudeFormat(uint max_characters, int significant_digits, bool adjust_positives) {
		NumberFormat result;
		result.max_characters = max_characters;
		result.significant_digits = significant_digits;
		result.adjust_positives = adjust_positives;
		result.truncate_zeros = false;
		result.max_mantissa = max_characters;
		result.scientific_notation = false;
		return result;
	}

	NumberFormat NumberFormat::crudeFormatPositive(uint max_characters, int significant_digits) {
		NumberFormat result;
		result.max_characters = max_characters;
		result.significant_digits = significant_digits;
		result.adjust_positives = false;
		result.truncate_zeros = false;
		result.max_mantissa = max_characters;
		result.scientific_notation = false;
		return result;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	std::string intToString(i64 value, uint digits = 0) {
		if (digits) {
			i64 digit = pow((i64)10, digits - 1);

			std::string output;
			output.resize(digits);

			for (int i = 0; digit > 0; digit /= 10, ++i) {
				output[i] = '0' + (value / digit) % 10;
			}
			return output;
		}
		else {
			if (!value) return "0";
			uint exponent = 0;
			i64 digit = 1;
			while (digit <= value) {
				digit *= 10;
				++exponent;
			}
			digit /= 10;

			std::string output;
			output.resize(exponent);

			for (int i = 0; digit > 0; digit /= 10, ++i) {
				output[i] = '0' + (value / digit) % 10;
			}
			return output;
		}
	}

	void truncateZeros(std::string& str) {
		uint index = str.size();
		for (int i = index - 1; i >= 0; --i) {
			if (str[i] != '0') break;
			index = i;
		}
		str = safeSubstr(str, 0, index);
	}

	std::string toFormattedString(double x, const NumberFormat& format) {
		if (isnan(x)) return "NaN";
		if (isinf(x)) return "inf";
		std::string output;
		if (x < 0) {
			output += '-';
			x *= -1.0;
		}
		else if (format.adjust_positives) output += ' ';
		double mantissa_upper = pow(10.0, format.max_mantissa);
		double mantissa_lower = pow(0.1, format.max_mantissa);
		if (x == 0.0) {
			return "0." + std::string(format.significant_digits > 0 ? (format.significant_digits - 1) : -format.significant_digits, '0');
		}
		if ((x >= mantissa_upper || x <= mantissa_lower) && format.scientific_notation) {
			std::string exponent_str = "E";
			int exponent = 0;
			while (x >= 10.0) {
				x /= 10.0;
				++exponent;
			}
			while (x < 1.0) {
				x *= 10.0;
				--exponent;
			}
			exponent_str += exponent < 0 ? '-' : '+';
			exponent_str += intToString(std::abs(exponent));

			uint decimals = min((int)format.max_characters - (int)output.size() - (int)exponent_str.size() - 2, format.significant_digits > 0 ? (format.significant_digits - 1) : -format.significant_digits);

			x += 0.5 * pow(0.1, decimals);

			output += intToString(x);

			if (x - (u64)x != 0.0 && decimals) {
				output += '.';

				std::string mantissa_string = intToString((x - (u64)x) * pow(10.0, decimals), decimals);

				if (format.truncate_zeros) truncateZeros(mantissa_string);

				if (decimals) output += mantissa_string;
			}

			output += exponent_str;
		}
		else {
			if (format.significant_digits < 0) {
				x += 0.5 * pow(0.1, -format.significant_digits);
			}

			int exponent = 0;
			while (x >= 10.0) {
				x /= 10.0;
				++exponent;
			}
			while (x < 1.0) {
				x *= 10.0;
				--exponent;
			}

			if (format.significant_digits > 0) {
				x += 5.0 * pow(0.1, format.significant_digits);
			}
			while (x >= 10.0) {
				x /= 10.0;
				++exponent;
			}

			std::string full_string = intToString(x * pow(10.0, format.significant_digits < 0 ? (-format.significant_digits + exponent) : format.significant_digits - 1));

			if (exponent < 0) output += '0';
			else output += safeSubstr(full_string, 0, exponent + 1);

			uint decimals = min((int)format.max_characters - (int)output.size() - 1, format.significant_digits > 0 ? (format.significant_digits - exponent - 1) : -format.significant_digits);

			if (x - (i64)x != 0.0) {
				output += '.';

				for (int i = 1; i < -exponent && output.size() < format.max_characters; ++i) {
					output += '0';
				}

				std::string decimal_string = safeSubstr(full_string, max(0, exponent + 1), max(0, exponent + (int)decimals + 1));

				if (format.truncate_zeros) truncateZeros(decimal_string);

				output += decimal_string;
			}
		}
		return output;
	}

	std::string toFormattedString(i64 x, const NumberFormat& format) {
		std::string output;
		if (x < 0) {
			output += '-';
			x *= -1.0;
		}
		else if (format.adjust_positives) output += ' ';

		uint exponent = 0;
		for (i64 x2 = x; x2 > 9; x2 /= 10) {
			++exponent;
		}

		if ((exponent + 1 > format.significant_digits || exponent + 1 + output.size() > format.max_characters) && format.scientific_notation) {
			std::string exponent_string = exponent > 0 ? "E+" : "E-";
			exponent_string += intToString(exponent);

			int decimals = min((int)format.max_characters - (int)exponent_string.size() - 2, format.significant_digits > 0 ? (format.significant_digits - 1) : -format.significant_digits);

			i64 aux = pow((i64)10, exponent);
			i64 leading = x / aux;
			i64 mantissa = x % aux;

			mantissa += 5 * pow((i64)10, exponent - decimals - 1);

			for (int i = 0; i < (int)exponent - (int)decimals; ++i) mantissa /= 10;

			output += intToString(leading) + '.' + intToString(mantissa) + exponent_string;
		}
		else {
			output += intToString(x);
		}
		return output;
	}

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	FormattedStream::Cell& FormattedStream::getCurrentCell() {
		while (current_cell >= cells.size()) cells.push_back(Cell(NumberFormat(), TextAlignment::left, 16));
		return cells[current_cell];
	}

	void FormattedStream::flush() {
		while (cells.size() <= content.size()) nextEntry();

		cells[cells.size() - 1].width = 2;

		bool has_content = content.size();

		if (finish_hline) {
			if (hline_string.size()) {
				for (int i = 0; i < vertical_padding; ++i) {
					(*write_to) << std::string(left_offset, ' ');
					for (int j = 0; j < hline_string.size(); ++j) {
						if (hline_string[j] == vertical_line || hline_string[j] == line_intersection) (*write_to) << vertical_line;
						else (*write_to) << ' ';
					}
					(*write_to) << '\n';
				}
			}

			int width = 0;
			int x = 0;
			for (int i = 0; i < cells.size(); ++i) {
				int cell_width = (i == cells.size() - 1) ? 1 : cells[i].width + horizontal_padding * 2 + 1;
				width += cell_width;

				if (hline_string.size() < width) hline_string.resize(width, ' ');

				if (cells[i].separator_before) {
					bool intersects = cells[i].hline_top || hline_string[x] == horizontal_line || hline_string[x] == line_intersection;
					if (i > 0) {
						intersects |= cells[i - 1].hline_top;
					}

					if (intersects) hline_string[x] = line_intersection;
					else hline_string[x] = vertical_line;
				}
				else if (cells[i].hline_top) hline_string[x] = (hline_string[x] == vertical_line) ? line_intersection : horizontal_line;
				else if (i > 0) {
					if (cells[i - 1].hline_top) hline_string[x] = (hline_string[x] == vertical_line) ? line_intersection : horizontal_line;
				}

				if (i == cells.size() - 1) break;

				if (cells[i].hline_top) {
					for (int i = x + 1; i < width; ++i) {
						hline_string[i] = (hline_string[i] == vertical_line || hline_string[i] == line_intersection) ? line_intersection : horizontal_line;
					}
				}

				x += cell_width;
			}

			(*write_to) << std::string(left_offset, ' ') << hline_string << '\n';

			hline_string = "";

			finish_hline = false;

			if (has_content) {
				for (int i = 0; i < vertical_padding; ++i) {
					(*write_to) << std::string(left_offset, ' ');
					for (int j = 0; j < cells.size(); ++j) {
						if (cells[j].separator_before) (*write_to) << vertical_line;
						else (*write_to) << " ";

						if (j == cells.size() - 1) break;

						(*write_to) << std::string(cells[j].width + horizontal_padding * 2, ' ');
					}
					(*write_to) << '\n';
				}
			}
		}

		hline_string = "";
		for (int i = 0; i < cells.size(); ++i) {
			if (cells[i].separator_before) {
				bool intersects = cells[i].hline_bottom;
				if (i > 0) {
					intersects |= cells[i - 1].hline_bottom;
				}

				if (intersects) hline_string += line_intersection;
				else hline_string += vertical_line;
			}
			else if (cells[i].hline_bottom) hline_string += horizontal_line;
			else if (i > 0) {
				if (cells[i - 1].hline_bottom) hline_string += horizontal_line;
				else hline_string += ' ';
			}
			else hline_string += ' ';

			if (i == cells.size() - 1) break;

			if (cells[i].hline_bottom) hline_string += std::string(cells[i].width + horizontal_padding * 2, horizontal_line);
			else hline_string += std::string(cells[i].width + horizontal_padding * 2, ' ');
		}
		if (start_hline) {
			start_hline = false;
			finish_hline = true;
		}

		while (has_content) {
			(*write_to) << std::string(left_offset, ' ');

			has_content = false;
			for (int i = 0; i < content.size(); ++i) {
				uint linebreak_index = 0;
				bool midline_break = false;
				bool break_on_space = false;
				uint space_number = 0;
				std::string& cell_content = content[i];
				uint max_width = cells[i].width;
				uint chars_since_space = 5;
				std::string line;

				if (cell_content.size() > max_width) {
					for (int j = 0; j < cell_content.size() && j < max_width; ++j) {
						char c = cell_content[j];
						if (c == ' ') {
							++space_number;
							chars_since_space = 0;
							midline_break = false;
							break_on_space = true;
							if (safeGetChar(cell_content, j - 1) != ' ') {
								linebreak_index = j;
							}
						}
						else if (c == '\n') {
							midline_break = false;
							break_on_space = false;
							space_number = 0;
							linebreak_index = j;
							break;
						}
						else if (space_number < 1) {
							linebreak_index = j + 1;
							midline_break = true;
							break_on_space = false;
						}
						else if ((max_width - j) / (max(2, (int)space_number) - 1) > 0 && chars_since_space > 2) {
							if (isLetter(c) && isLetter(cell_content[j + 1])) {
								linebreak_index = j + 1;
								midline_break = true;
								break_on_space = false;
							}
						}
						++chars_since_space;
					}

					if (break_on_space) --space_number;

					if (midline_break) {
						cell_content.insert(cell_content.begin() + linebreak_index - 1, ' ');
						cell_content.insert(cell_content.begin() + linebreak_index - 1, '-');
					}

					line = safeSubstr(cell_content, 0, linebreak_index);
					cell_content = safeSubstr(cell_content, linebreak_index + 1, cell_content.size());

					for (int i = 0; i < cell_content.size(); ++i) {
						if (cell_content[i] != ' ') has_content = true;
					}
				}
				else {
					line = cell_content;
					cell_content = "";
				}

				if (cells[i].separator_before) (*write_to) << vertical_line;
				else (*write_to) << " ";

				(*write_to) << std::string(horizontal_padding, ' ');

				switch (cells[i].alignment) {
				case TextAlignment::left:
					(*write_to) << line;
					for (int j = 0; j < max_width - line.size(); ++j) {
						(*write_to) << ' ';
					}
					break;
				case TextAlignment::right:
					for (int j = 0; j < max_width - line.size(); ++j) {
						(*write_to) << ' ';
					}
					(*write_to) << line;
					break;
				case TextAlignment::centered:
					for (int j = 0; j < (max_width - line.size()) / 2; ++j) {
						(*write_to) << ' ';
					}
					(*write_to) << line;
					for (int j = 0; j < (max_width - line.size() + 1) / 2; ++j) {
						(*write_to) << ' ';
					}
					break;
				case TextAlignment::block:
					if (space_number) {
						uint space_index = 0;
						uint clearance = max_width - line.size();
						for (int j = 0; j < line.size(); ++j) {
							(*write_to) << line[j];
							if (line[j] == ' ') {
								uint count = (clearance + space_index) / space_number;
								for (int k = 0; k < count; ++k) {
									(*write_to) << ' ';
								}
								++space_index;
							}
						}
					}
					else {
						(*write_to) << line;
						for (int j = 0; j < max_width - line.size(); ++j) {
							(*write_to) << ' ';
						}
					}
					break;
				}

				(*write_to) << std::string(horizontal_padding, ' ');
			}
			if (cells.size() > content.size()) if (cells[cells.size() - 1].separator_before) (*write_to) << std::string(1, vertical_line);
			(*write_to) << '\n';
		}
		content.clear();
		current_cell = 0;
		row_has_content = false;

		for (int i = 0; i < cells.size(); ++i) {
			cells[i].hline_bottom = false;
			cells[i].hline_top = false;
			cells[i].separator_before = false;
		}
	}

	void FormattedStream::resetRows() {
		flush();
		current_cell = 0;
		column_count = 0;
		cells.clear();
	}

	void FormattedStream::addRow(const NumberFormat& number_format, TextAlignment alignment, uint width) {
		cells.push_back(Cell(number_format, alignment, width));
		++column_count;
	}

	char FormattedStream::formatRow(const NumberFormat& number_format) {
		getCurrentCell().number_format = number_format;
		return ' ';
	}

	char FormattedStream::formatRow(TextAlignment alignment) {
		getCurrentCell().alignment = alignment;
		return ' ';
	}

	uint FormattedStream::getColumnCount() {
		return column_count;
	}

	void FormattedStream::offsetRight(uint offset) {
		left_offset = offset;
	}

	void FormattedStream::centerTable(uint available) {
		uint total_width = 1;
		for (int i = 0; i < cells.size(); ++i) {
			total_width += cells[i].width + horizontal_padding * 2 + 1;
		}
		left_offset = max(0, ((int)available - (int)total_width) / 2);
	}

	void FormattedStream::nextEntry() {
		++current_cell;
		getCurrentCell().separator_before = false;
	}

	FormattedStream& FormattedStream::operator<<(const std::string str) {
		row_has_content = true;
		if (current_cell >= content.size()) content.resize(current_cell + 1);
		getCurrentCell();
		content[current_cell] += str;
		return *this;
	}

	FormattedStream& FormattedStream::operator<<(char c) {
		if (c == '-') {
			getCurrentCell();
			finish_hline = true;
			cells[current_cell].hline_top = true;
			return *this;
		}
		else if (c == '_') {
			getCurrentCell();
			start_hline = true;
			cells[current_cell].hline_bottom = true;
			return *this;
		}
		else if (c == '\n') {
			if (!cells.size()) (*write_to) << '\n';
			else flush();
			return *this;
		}
		else if (c == '\r') {
			start_hline = true;
			flush();
			return *this;
		}
		else if (c == ',') {
			nextEntry();
			return *this;
		}
		else if (c == '|') {
			getCurrentCell();
			if (row_has_content) nextEntry();
			row_has_content = true;
			cells[current_cell].separator_before = true;
			return *this;
		}
		else if (c == ' ') return *this;
		return operator<<(std::string(1, c));
	}

	FormattedStream& FormattedStream::operator<<(double number) {
		return operator<<(toFormattedString(number, getCurrentCell().number_format));
	}

	FormattedStream& FormattedStream::operator<<(i64 number) {
		return operator<<(toFormattedString(number, getCurrentCell().number_format));
	}
}