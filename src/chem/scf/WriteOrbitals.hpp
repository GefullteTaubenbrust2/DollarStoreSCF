#pragma once
#include <string>

namespace scf {
	void writeOrbitalsMolden(const std::string& path, bool write_virtual = true);

	void writeOrbitalsWFX(const std::string& path, bool write_virtual = false);
}