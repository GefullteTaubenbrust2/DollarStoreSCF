#pragma once
#include "Molecule.hpp"

namespace flo {
	struct PointGroup {
		uint rotation = 1;
		bool inversion = false;
		bool twist = false;
		bool spherical = false;
		bool vertical_mirror = false;
		bool horizontal_c2 = false;
		bool horizontal_mirror = false;

		PointGroup() = default;

		constexpr PointGroup(uint rotation, bool inversion, bool twist, bool spherical, bool vertical_mirror, bool horizontal_c2, bool horizontal_mirror);

		static constexpr PointGroup Cn(uint n);
		static constexpr PointGroup Cnh(uint n);
		static constexpr PointGroup Cnv(uint n);
		static constexpr PointGroup Dn(uint n);
		static constexpr PointGroup Dnh(uint n);
		static constexpr PointGroup Dnd(uint n);
		static constexpr PointGroup Sn(uint n);
		static constexpr PointGroup C1();
		static constexpr PointGroup Ci();
		static constexpr PointGroup Cs();
		static constexpr PointGroup Cinfv();
		static constexpr PointGroup Dinfh();
		static constexpr PointGroup Td();
		static constexpr PointGroup Th();
		static constexpr PointGroup Oh();
		static constexpr PointGroup I();
		static constexpr PointGroup Ih();

		std::string getSchoenfliesSymbol();
	};

	PointGroup findPointGroup(Molecule& molecule, double tolerance);
}