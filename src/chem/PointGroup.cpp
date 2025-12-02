#include "PointGroup.hpp"

namespace flo {
	constexpr PointGroup::PointGroup(uint rotation, bool inversion, bool twist, bool spherical, bool vertical_mirror, bool horizontal_c2, bool horizontal_mirror) :
	rotation(rotation), inversion(inversion), twist(twist), spherical(spherical), vertical_mirror(vertical_mirror), horizontal_c2(horizontal_mirror), horizontal_mirror(horizontal_mirror) {
	}

	constexpr PointGroup PointGroup::Cn(uint n) {
		return PointGroup(n, false, false, false, false, false, false);
	}
	constexpr PointGroup PointGroup::Cnh(uint n) {
		return PointGroup(n, n % 2 == 0, true, false, false, false, true);
	}
	constexpr PointGroup PointGroup::Cnv(uint n) {
		return PointGroup(n, false, false, false, true, false, false);
	}
	constexpr PointGroup PointGroup::Dn(uint n) {
		return PointGroup(n, false, false, false, false, true, false);
	}
	constexpr PointGroup PointGroup::Dnh(uint n) {
		return PointGroup(n, n % 2 == 0, true, false, true, true, true);
	}
	constexpr PointGroup PointGroup::Dnd(uint n) {
		return PointGroup(n, true, true, false, true, true, false);
	}
	constexpr PointGroup PointGroup::Sn(uint n) {
		return PointGroup(n, false, true, false, false, false, false);
	}
	constexpr PointGroup PointGroup::C1() {
		return PointGroup(1, false, false, false, false, false, false);
	}
	constexpr PointGroup PointGroup::Ci() {
		return PointGroup(1, true, false, false, false, false, false);
	}
	constexpr PointGroup PointGroup::Cs() {
		return PointGroup(1, false, false, false, false, false, true);
	}
	constexpr PointGroup PointGroup::Cinfv() {
		return PointGroup(0, false, true, false, true, false, false);
	}
	constexpr PointGroup PointGroup::Dinfh() {
		return PointGroup(0, true, true, false, true, true, true);
	}
	constexpr PointGroup PointGroup::Td() {
		return PointGroup(3, false, false, true, true, false, false);
	}
	constexpr PointGroup PointGroup::Th() {
		return PointGroup(3, true, true, true, false, false, false);
	}
	constexpr PointGroup PointGroup::Oh() {
		return PointGroup(4, true, true, true, true, true, true);
	}
	constexpr PointGroup PointGroup::I() {
		return PointGroup(5, false, false, true, false, true, false);
	}
	constexpr PointGroup PointGroup::Ih() {
		return PointGroup(5, true, true, true, false, true, false);
	}

	std::string PointGroup::getSchoenfliesSymbol() {
		if (!rotation) {
			if (inversion) return "Dinfh";
			return "Cinfv";
		}
		else if (spherical) {
			switch (rotation) {
			case 3:
				if (inversion) return "Th";
				return "Td";
			case 4:
				return "Oh";
			case 5:
				if (inversion) return "Ih";
				return "I";
			default:
				return "invalid";
				break;
			}
		}
		else if (rotation == 1) {
			if (horizontal_mirror) return "Cs";
			else if (inversion) return "Ci";
			else return "C1";
		}
		else if (horizontal_c2) {
			if (horizontal_mirror) return "D" + std::to_string(rotation) + "h";
			else if (vertical_mirror) return "D" + std::to_string(rotation) + "d";
			else return "D" + std::to_string(rotation);
		}
		else {
			if (horizontal_mirror) return "C" + std::to_string(rotation) + "h";
			else if (vertical_mirror) return "C" + std::to_string(rotation) + "v";
			else if (twist) "S" + std::to_string(2 * rotation);
			else return "C" + std::to_string(rotation);
		}
	}

	struct RotationAxis {
		uint order = 2;
		vec3 direction;
	};

	void addUnique(const RotationAxis& axis, std::vector<RotationAxis>& list, double tolerance) {
		bool added = false;
		for (int i = 0; i < list.size(); ++i) {
			if (length(list[i].direction - axis.direction) < tolerance || length(list[i].direction + axis.direction) < tolerance) {
				if (list[i].order < axis.order) list[i].order = axis.order;
				added = true;
			}
		}
		if (!added) list.push_back(axis);
	}

	void addUnique(const vec3& dir, int order, std::vector<RotationAxis>& list, double tolerance) {
		RotationAxis axis;
		axis.direction = dir;
		axis.order = order;
		addUnique(axis, list, tolerance);
	}

	void addUnique(const vec3& axis, std::vector<vec3>& list, double tolerance) {
		bool added = false;
		for (int i = 0; i < list.size(); ++i) {
			if (length(list[i] - axis) < tolerance || length(list[i] + axis) < tolerance) {
				added = true;
			}
		}
		if (!added) list.push_back(axis);
	}

	PointGroup findPointGroup(Molecule& molecule, double tolerance) {
		if (molecule.size() < 2) return PointGroup(0, true, true, true, true, true, true);

		PointGroup result;
		result.rotation = 0;
		result.inversion = true;

		vec3 base_dir = normalize(molecule[1].position - molecule[0].position);

		// Check for linearity
		for (int i = 2; i < molecule.size(); ++i) {
			vec3 offset = molecule[i].position - molecule[0].position;
			double deviation = length(offset - base_dir * dot(offset, base_dir));
			if (deviation > tolerance) {
				result.rotation = 1;
				break;
			}
		}

		vec3 center_of_mass;
		double total_mass = 0.0;
		for (int i = 0; i < molecule.size(); ++i) {
			center_of_mass += molecule[i].mass * molecule[i].position;
			total_mass += molecule[i].mass;
		}
		center_of_mass /= total_mass;

		std::vector<std::vector<uint>> equidistant_groups;

		for (int i = 0; i < molecule.size(); ++i) {
			double dist = length(molecule[i].position - center_of_mass);
			bool sorted = false;
			for (int j = 0; j < equidistant_groups.size(); ++j) {
				uint index = equidistant_groups[j][0];
				if (molecule[i].element == molecule[index].element && std::abs(length(molecule[index].position - center_of_mass) - dist) < tolerance) {
					equidistant_groups[j].push_back(i);
					sorted = true;
					break;
				}
			}
			if (!sorted) {
				equidistant_groups.push_back(std::vector<uint>());
				equidistant_groups[equidistant_groups.size() - 1].push_back(i);
			}
		}

		// Check for inversion
		for (int i = 0; i < equidistant_groups.size() && result.inversion; ++i) {
			for (int j = 0; j < equidistant_groups[i].size(); ++j) {
				vec3 pos = molecule[equidistant_groups[i][j]].position - center_of_mass;

				bool invertible = false;
				for (int k = 0; k < equidistant_groups[i].size(); ++k) {
					if (length(pos + molecule[equidistant_groups[i][k]].position - center_of_mass) < tolerance) {
						invertible = true;
						break;
					}
				}

				if (!invertible) {
					i = equidistant_groups.size();
					result.inversion = false;
					break;
				}
			}
		}

		// Linear point groups
		if (!result.rotation) {
			if (result.inversion) return PointGroup::Dinfh();
			else return PointGroup::Cinfv();
		}

		std::vector<RotationAxis> axes;

		int biggest_group_index = 0;
		int biggest_size = 0;
		for (int i = 0; i < equidistant_groups.size(); ++i) {
			if (equidistant_groups[i].size() > biggest_size) {
				biggest_group_index = i;
				biggest_size = equidistant_groups[i].size();
			}
		}

		std::vector<uint>& group = equidistant_groups[biggest_group_index];

		// Search for axes of rotation
		std::vector<RotationAxis> candidates;
		// We first look for atom-aligned axes
		int involved_atoms = 0;
		for (int j = 0; j < group.size(); ++j) {
			vec3 axis = normalize(molecule[group[j]].position - center_of_mass);
			for (int k = 0; k < group.size(); ++k) {
				if (j == k) continue;

				vec3 vec = molecule[group[k]].position - center_of_mass;
				double overlap = dot(vec, axis);

				if (length(vec - axis * overlap) > tolerance) {
					++involved_atoms;
				}
			}
			if (involved_atoms > 1) {
				addUnique(axis, involved_atoms, candidates, tolerance);
			}
		}
		// Then we check for C2-related tuples
		for (int j = 0; j < group.size(); ++j) {
			for (int k = j + 1; k < group.size(); ++k) {
				vec3 midpoint = 0.5 * (molecule[group[j]].position + molecule[group[k]].position) - center_of_mass;
				vec3 dir = normalize(molecule[group[j]].position - molecule[group[k]].position);
				if (std::abs(dot(midpoint, dir)) > tolerance) continue;
				if (length(midpoint) < tolerance) {
					for (int l = 0; l < molecule.size(); ++l) {
						if (l == group[k] || l == group[k]) continue;
						vec3 pos = molecule[j].position - center_of_mass;
						if (std::abs(dot(pos, dir)) < tolerance) {
							addUnique(normalize(pos), 2, candidates, tolerance);
						}
						if (length(pos - dir * dot(dir, pos)) > tolerance) {
							addUnique(normalize(cross(pos, dir)), 2, candidates, tolerance);
						}
					}
				}
				else {
					addUnique(normalize(midpoint), 2, candidates, tolerance);
				}
			}
		}
		// Finally we check all triplets of points
		for (int j = 0; j < group.size(); ++j) {
			for (int k = j + 1; k < group.size(); ++k) {
				for (int l = k + 1; l < group.size(); ++l) {
					vec3 axis = normalize(cross(molecule[group[l]].position - molecule[group[k]].position, molecule[group[j]].position - molecule[group[k]].position));
					int out_of_plane_atoms = 0;
					involved_atoms = 3;
					for (int m = 0; m < group.size(); ++m) {
						if (m == j || m == k || m == l) continue;

						vec3 vec1 = molecule[group[m]].position - molecule[group[k]].position;
						double overlap1 = dot(vec1, axis);
						vec3 vec2 = molecule[group[m]].position - center_of_mass;
						double overlap2 = dot(vec2, axis);

						if (std::abs(overlap1) < tolerance) {
							++involved_atoms;
						}
						else if (length(vec2 - axis * overlap2) > tolerance) {
							++out_of_plane_atoms;
						}
					}
					if (out_of_plane_atoms % involved_atoms == 0) {
						addUnique(axis, involved_atoms, candidates, tolerance);
					}
				}
			}
		}
			
		// Now check that the resulting axes are locally valid symmetry elements
		for (RotationAxis candidate : candidates) {
			int order = candidate.order;
			for (; order > 1; --order) {
				double angle = 2.0 * PI / order;
				bool found = false;
				for (int j = 0; j < group.size(); ++j) {
					vec3 pos = molecule[group[j]].position - center_of_mass;
					vec3 axial = candidate.direction * dot(candidate.direction, pos);
					vec3 radial = pos - axial;
					double radius = length(radial);
					if (radius) radial /= radius;
					else radial *= 0.0;
					vec3 tangential = cross(radial, candidate.direction);

					vec3 rotated_position = axial + radius * (std::cos(angle) * radial + std::sin(angle) * tangential);

					found = false;

					for (int k = 0; k < group.size(); ++k) {
						pos = molecule[group[k]].position - center_of_mass;
						if (length(pos - rotated_position) < tolerance) {
							found = true;
							break;
						}
					}

					if (!found) break;
				}
				if (found) break;
			}
			if (order > 1) {
				// If the check passed, check that the symmetry is valid globally
				order = candidate.order;
				for (; order > 1; --order) {
					double angle = 2.0 * PI / order;
					bool found = false;
					for (int j = 0; j < molecule.size(); ++j) {
						vec3 pos = molecule[j].position - center_of_mass;
						vec3 axial = candidate.direction * dot(candidate.direction, pos);
						vec3 radial = pos - axial;
						double radius = length(radial);
						if (radius) radial /= radius;
						else radial *= 0.0;
						vec3 tangential = cross(radial, candidate.direction);

						vec3 rotated_position = axial + radius * (std::cos(angle) * radial + std::sin(angle) * tangential);

						found = false;

						for (int k = 0; k < molecule.size(); ++k) {
							pos = molecule[k].position - center_of_mass;
							if (length(pos - rotated_position) < tolerance) {
								found = true;
								break;
							}
						}

						if (!found) break;
					}
					if (found) break;
				}
			}
			if (order > 1) {
				candidate.order = order;
				addUnique(candidate, axes, tolerance);
			}
		}

		// Search for mirror planes
		std::vector<vec3> mirrors_local;
		std::vector<vec3> mirror_planes;
		for (int i = 0; i < group.size(); ++i) {
			for (int j = i + 1; j < group.size(); ++j) {
				vec3 axis = normalize(molecule[group[i]].position - molecule[group[j]].position);
				addUnique(axis, mirrors_local, tolerance);
				axis = normalize(cross(molecule[group[i]].position - molecule[group[j]].position, molecule[group[i]].position - center_of_mass));
				addUnique(axis, mirrors_local, tolerance);
			}
		}
		for (int i = 0; i < mirrors_local.size(); ++i) {
			vec3 axis = mirrors_local[i];
			// Check locally first
			bool symmetric = false;
			for (int j = 0; j < group.size(); ++j) {
				symmetric = false;
				vec3 mirrored_pos = molecule[group[j]].position - center_of_mass;
				mirrored_pos -= 2.0 * axis * dot(axis, mirrored_pos);
				for (int k = 0; k < group.size(); ++k) {
					vec3 pos = molecule[group[k]].position - center_of_mass;
					if (length(pos - mirrored_pos) < tolerance) {
						symmetric = true;
						break;
					}
				}
				if (!symmetric) break;
			}
			// Check globally
			if (symmetric) {
				for (int j = 0; j < molecule.size(); ++j) {
					symmetric = false;
					vec3 mirrored_pos = molecule[j].position - center_of_mass;
					mirrored_pos -= 2.0 * axis * dot(axis, mirrored_pos);
					for (int k = 0; k < molecule.size(); ++k) {
						vec3 pos = molecule[k].position - center_of_mass;
						if (length(pos - mirrored_pos) < tolerance) {
							symmetric = true;
							break;
						}
					}
					if (!symmetric) break;
				}
			}
			if (symmetric) addUnique(axis, mirror_planes, tolerance);
		}
		int highest_rotation_index = 0;
		int highest_rotation = 1;
		for (int i = 0; i < axes.size(); ++i) {
			if (axes[i].order > highest_rotation) {
				highest_rotation = axes[i].order;
				highest_rotation_index = i;
			}
		}
		result.rotation = highest_rotation;

		vec3 principle_axis(0.0);
		if (highest_rotation > 1) principle_axis = axes[highest_rotation_index].direction;

		for (int i = 0; i < axes.size(); ++i) {
			if (i == highest_rotation_index) continue;
			if (axes[i].order == highest_rotation) {
				result.spherical = true;
			}
			if (axes[i].order == 2 && std::abs(dot(axes[i].direction, principle_axis)) < tolerance) {
				result.horizontal_c2 = true;
			}
		}

		for (vec3 mirror_plane : mirror_planes) {
			if (std::abs(dot(mirror_plane, principle_axis)) < tolerance) {
				result.vertical_mirror = true;
			}
			if (length(mirror_plane - principle_axis) < tolerance || length(mirror_plane + principle_axis) < tolerance) {
				result.horizontal_mirror = true;
			}
		}

		if (highest_rotation > 1) {
			double angle = PI / highest_rotation;
			bool symmetric = false;
			for (int i = 0; i < molecule.size(); ++i) {
				vec3 twisted_position = molecule[i].position - center_of_mass;
				vec3 axial = principle_axis * dot(principle_axis, twisted_position);
				vec3 radial = twisted_position - axial;
				double radius = length(radial);
				if (radius) radial /= radius;
				else radial = 0.0;
				vec3 tangential = cross(principle_axis, radial);

				twisted_position = -1.0 * axial + radius * (std::cos(angle) * radial + std::sin(angle) * tangential);

				symmetric = false;

				for (int j = 0; j < molecule.size(); ++j) {
					vec3 pos = molecule[j].position - center_of_mass;

					if (length(pos - twisted_position) < tolerance) {
						symmetric = true;
						break;
					}
				}

				if (!symmetric) break;
			}

			if (symmetric) {
				result.twist = true;
			}
		}

		return result;
	}
}
