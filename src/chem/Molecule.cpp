#include "Molecule.hpp"
#include "Constants.hpp"

#include "../util/FormattedStream.hpp"

#include <sstream>

namespace flo {
	const int affected_by_bond_length = 1;
	const int affected_by_dihedral_angle = 2;

	Atom::Atom(const vec3& position, Element element, double mass) :
	position(position), element(element), charge((uint)element), mass(mass) {
	}

	Atom::Atom(const vec3& position, Element element) :
	position(position), element(element), charge((uint)element), mass(common_atomic_masses[(int)element] * Da_to_me) {
	}

	vec3 Atom::getRepulsionGradient(const vec3& pos) {
		vec3 x = pos - position;

		double denominator = std::pow(length2(x), -1.5);

		return -((double)charge)* x * denominator;
	}

	size_t Molecule::size() const {
		return atoms.size();
	}
	Atom& Molecule::operator[](size_t index) {
		return atoms[index];
	}

	const Atom& Molecule::operator[](size_t index) const {
		return atoms[index];
	}

	std::string Molecule::getAtomName(size_t index) const {
		const Atom& a = atoms[index];
		std::ostringstream oss;
		oss << element_symbols[(int)a.element - 1] << (index + 1);
		return oss.str();
	}

	struct Node {
		std::vector<Node*> neighbours;
		uint index = 0;
		uint neighbour_index = 0;
		bool in_tree = false;

		void addNeighbour(Node* neighbour) {
			neighbour->neighbour_index = neighbours.size();
			neighbours.push_back(neighbour);
			neighbour->in_tree = true;
			neighbour->neighbours.push_back(this);
		}
	};

	size_t Molecule::getInternalCoordinateCount() const {
		size_t result = internal_coordinate_indices.size() / 4 * 3 - 2;
		if (internal_coordinate_indices.size() > 4) --result;
		return result;
	}

	void Molecule::generateInternalCoordinates(std::vector<uint>* bond_map_ptr) {
		base_position = atoms[0].position;
		if (size() > 1) base_direction = normalize(atoms[1].position - atoms[0].position);
		if (size() > 2) base_normal = normalize(atoms[2].position - atoms[0].position - base_direction * dot(base_direction, atoms[2].position - atoms[0].position));

		std::vector<uint> backup;
		if (!bond_map_ptr) bond_map_ptr = &backup;
		std::vector<uint>& bond_map = *bond_map_ptr;

		if (!bond_map.size()) {
			for (int i = 0; i < size(); ++i) {
				Atom atom0 = atoms[i];
				for (int j = i + 1; j < size(); ++j) {
					Atom atom1 = atoms[j];

					double bond_length = length(atom0.position - atom1.position);

					if (bond_length < 1.5 * A_to_a0 * (covalent_radii_A[(int)atom0.element - 1] + covalent_radii_A[(int)atom1.element - 1])) {
						bond_map.push_back(i);
						bond_map.push_back(j);
					}
				}
			}
		}

		std::vector<Node> tree;
		std::vector<Node*> ordered_tree;
		tree.resize(size());
		ordered_tree.reserve(size());
		for (int i = 0; i < size(); ++i) {
			tree[i].index = i;
		}
		tree[0].in_tree = true;
		ordered_tree.push_back(&tree[0]);

		std::vector<Node*> current_branches{&tree[0]}, next_branches;

		while (current_branches.size()) {
			for (int i = 0; i < current_branches.size(); ++i) {
				uint node_a = current_branches[i]->index;
				for (int j = 0; j < bond_map.size() / 2; ++j) {
					uint bond_a = bond_map[2 * j];
					uint bond_b = bond_map[2 * j + 1];
					if (node_a == bond_a || node_a == bond_b) {
						for (int k = 0; k < size(); ++k) {
							uint node_b = tree[k].index;
							if ((node_b == bond_a || node_b == bond_b) && !tree[k].in_tree) {
								current_branches[i]->addNeighbour(&tree[k]);
								ordered_tree.push_back(&tree[k]);
								next_branches.push_back(&tree[k]);
								break;
							}
						}
					}
				}
			}
			current_branches = next_branches;
			next_branches.clear();
		}

		affected_indices.resize(size() * (size() - 1));

		internal_coordinate_indices.clear();
		internal_coordinate_indices.reserve((size() - 1) * 4);

		internal_coordinate_indices.push_back(0);
		internal_coordinate_indices.push_back(0);
		internal_coordinate_indices.push_back(ordered_tree[0]->index);
		internal_coordinate_indices.push_back(ordered_tree[1]->index);

		if (size() > 2) {
			internal_coordinate_indices.reserve(5);
			
			if (ordered_tree[0]->neighbours.size() > 1) {
				internal_coordinate_indices.push_back(0);
				internal_coordinate_indices.push_back(ordered_tree[1]->index);
				internal_coordinate_indices.push_back(ordered_tree[0]->index);
				internal_coordinate_indices.push_back(ordered_tree[2]->index);
				
			}
			else {
				internal_coordinate_indices.push_back(0);
				internal_coordinate_indices.push_back(ordered_tree[0]->index);
				internal_coordinate_indices.push_back(ordered_tree[1]->index);
				internal_coordinate_indices.push_back(ordered_tree[2]->index);
			}
		}

		for (int i = 0; i < min((int)size(), 3); ++i) {
			ordered_tree[i]->in_tree = true;
		}

		for (int i = 3; i < size(); ++i) {
			ordered_tree[i]->in_tree = false;
		}

		for (int i = 3; i < size(); ++i) {
			Node* node1 = ordered_tree[i];
			Node* node2 = node1->neighbours[0];
			Node* node3 = node2->neighbours[0];
			if (node3->index == node1->index || !node3->in_tree) {
				for (Node* node : node2->neighbours) {
					if (node->in_tree && node != node1) {
						node3 = node;
						break;
					}
				}
			}

			Node* node4 = nullptr;
			for (Node* node : node3->neighbours) {
				if (node->in_tree && node != node2) {
					node4 = node;
					break;
				}
			}
			if (!node4) {
				for (Node* node : node2->neighbours) {
					if (node->in_tree && node != node1 && node != node3) {
						node4 = node;
						break;
					}
				}
			}

			internal_coordinate_indices.push_back(node4->index);
			internal_coordinate_indices.push_back(node3->index);
			internal_coordinate_indices.push_back(node2->index);
			internal_coordinate_indices.push_back(node1->index);

			node1->in_tree = true;
		}

		for (int i = 0; i < size() - 1; ++i) {
			uint active_index1 = internal_coordinate_indices[i * 4 + 3];
			uint active_index0 = internal_coordinate_indices[i * 4 + 2];

			for (Node& node : tree) {
				node.in_tree = false;
			}

			current_branches.clear();
			next_branches.clear();

			current_branches.push_back(&tree[active_index0]);
			affected_indices[i * size() + active_index0] = 0;
			current_branches.push_back(&tree[active_index1]);
			affected_indices[i * size() + active_index1] = affected_by_bond_length;

			for (int j = 0; j < size() - 1; ++j) {
				if (i == j) continue;
				uint index0 = internal_coordinate_indices[j * 4];
				uint index1 = internal_coordinate_indices[j * 4 + 1];
				uint index2 = internal_coordinate_indices[j * 4 + 2];
				uint index3 = internal_coordinate_indices[j * 4 + 3];

				if (index2 == active_index0 && index0 == active_index1) {
					current_branches.push_back(&tree[index3]);
					affected_indices[i * size() + index3] = affected_by_dihedral_angle;
				}
			}

			for (Node* root : current_branches) {
				root->in_tree = true;
			}

			while (current_branches.size()) {
				for (Node* root: current_branches) {
					for (Node* branch : root->neighbours) {
						if (!branch->in_tree) {
							next_branches.push_back(branch);
							branch->in_tree = true;
							affected_indices[i * size() + branch->index] = affected_indices[i * size() + root->index];
						}
					}
				}
				current_branches = next_branches;
				next_branches.clear();
			}
		}
	}

	void Molecule::calculateInternalCoordinates(VectorNd& result) const {
		int index = 0;
		result.resize(getInternalCoordinateCount());
		for (int i = 0; i < internal_coordinate_indices.size() / 4; ++i) {
			const vec3& atom0 = atoms[internal_coordinate_indices[i * 4]].position;
			const vec3& atom1 = atoms[internal_coordinate_indices[i * 4 + 1]].position;
			const vec3& atom2 = atoms[internal_coordinate_indices[i * 4 + 2]].position;
			const vec3& atom3 = atoms[internal_coordinate_indices[i * 4 + 3]].position;

			vec3 bond0 = atom0 - atom1;
			vec3 bond1 = atom1 - atom2;
			vec3 bond2 = atom3 - atom2;

			double bond_length = length(bond2);
			double bond_angle = std::acos(dot(bond2, bond1) / (length(bond2) * length(bond1)));
			vec3 in_plane_1 = bond2 - bond1 * dot(bond2, bond1) / length2(bond1);
			vec3 in_plane_0 = bond0 - bond1 * dot(bond0, bond1) / length2(bond1);
			double error0 = dot(in_plane_0, bond1);
			double error1 = dot(in_plane_1, bond1);
			double sign = dot(cross(in_plane_0, in_plane_1), bond1) > 0.0 ? -1.0 : 1.0;
			double dihedral_angle = sign * std::acos(dot(in_plane_1, in_plane_0) / (length(in_plane_1) * length(in_plane_0)));

			result[index] = bond_length;
			++index;
			if (i == 0) continue;
			result[index] = bond_angle;
			++index;
			if (i == 1) continue;
			result[index] = dihedral_angle;
			++index;
		}
	}

	void Molecule::assignInternalCoordinates(const VectorNd& internal_coordinates, double* step_size) {
		VectorNd old;
		if (step_size) {
			old.resize(size() * 3);
			for (int i = 0; i < size(); ++i) {
				old[i * 3]     = atoms[i].position.x;
				old[i * 3 + 1] = atoms[i].position.y;
				old[i * 3 + 2] = atoms[i].position.z;
			}
		}

		atoms[0].position = base_position;
		if (size() > 1) atoms[internal_coordinate_indices[3]].position = base_position + internal_coordinates[0] * base_direction;
		if (size() > 2) {
			double sign = internal_coordinate_indices[6] == 0 ? 1.0 : -1.0;
			atoms[internal_coordinate_indices[7]].position =
				atoms[internal_coordinate_indices[6]].position + internal_coordinates[1] * (std::sin(internal_coordinates[2]) * base_normal + sign * std::cos(internal_coordinates[2]) * base_direction);
		}

		for (int i = 2; i < size() - 1; ++i) {
			vec3 atom0 = atoms[internal_coordinate_indices[4 * i    ]].position;
			vec3 atom1 = atoms[internal_coordinate_indices[4 * i + 1]].position;
			vec3 atom2 = atoms[internal_coordinate_indices[4 * i + 2]].position;

			vec3 bond0 = atom0 - atom1;
			vec3 bond1 = atom1 - atom2;

			double bond_length    = internal_coordinates[(i - 1) * 3];
			double bond_angle     = internal_coordinates[(i - 1) * 3 + 1];
			double dihedral_angle = internal_coordinates[(i - 1) * 3 + 2];

			vec3 forward_vector  = normalize(bond1);
			vec3 ecliptic_vector = normalize(bond0 - bond1 * dot(bond0, bond1) / length2(bond1));
			vec3 clinal_vector   = cross(ecliptic_vector, forward_vector);

			atoms[internal_coordinate_indices[4 * i + 3]].position = 
				atom2 + bond_length * (std::cos(bond_angle) * forward_vector + std::sin(bond_angle) * (std::cos(dihedral_angle) * ecliptic_vector + std::sin(dihedral_angle) * clinal_vector));
		}

		vec3 center_of_mass;
		double mass = 0.0;
		for (int i = 0; i < size(); ++i) {
			center_of_mass += atoms[i].mass * atoms[i].position;
			mass += atoms[i].mass;
		}
		center_of_mass /= mass;

		for (Atom& atom : atoms) {
			atom.position -= center_of_mass;
		}

		if (step_size) {
			double rms = 0.0;
			for (int i = 0; i < size(); ++i) {
				double dx = old[i * 3    ] - atoms[i].position.x;
				double dy = old[i * 3 + 1] - atoms[i].position.y;
				double dz = old[i * 3 + 2] - atoms[i].position.z;
				rms += dx * dx + dy * dy + dz * dz;
			}
			*step_size = std::sqrt(rms);
		}
	}

	const MatrixNd& Molecule::calculateDisplacementMatrix() {
		displacement_matrix.resize(getInternalCoordinateCount(), size() * 3);

		int coordinate_index = 0;
		for (int i = 0; i < internal_coordinate_indices.size() / 4; ++i) {
			vec3& atom0 = atoms[internal_coordinate_indices[i * 4]].position;
			vec3& atom1 = atoms[internal_coordinate_indices[i * 4 + 1]].position;
			vec3& atom2 = atoms[internal_coordinate_indices[i * 4 + 2]].position;
			vec3& atom3 = atoms[internal_coordinate_indices[i * 4 + 3]].position;

			vec3 bond0 = atom0 - atom1;
			vec3 bond1 = atom1 - atom2;
			vec3 bond2 = atom3 - atom2;

			for (int j = 0; j < size(); ++j) {
				if (affected_indices[i * size() + j] & affected_by_bond_length) {
					vec3 dir = normalize(bond2);
					displacement_matrix.at(j * 3    , coordinate_index) = dir.x;
					displacement_matrix.at(j * 3 + 1, coordinate_index) = dir.y;
					displacement_matrix.at(j * 3 + 2, coordinate_index) = dir.z;
				}
				else {
					displacement_matrix.at(j * 3    , coordinate_index) = 0.0;
					displacement_matrix.at(j * 3 + 1, coordinate_index) = 0.0;
					displacement_matrix.at(j * 3 + 2, coordinate_index) = 0.0;
				}
			}
			++coordinate_index;
			if (i) {
				vec3 axis = normalize(cross(bond2, bond1));

				for (int j = 0; j < size(); ++j) {
					if (affected_indices[i * size() + j] & affected_by_bond_length) {
						vec3 dir = cross(atoms[j].position - atom2, axis);
						displacement_matrix.at(j * 3    , coordinate_index) = dir.x;
						displacement_matrix.at(j * 3 + 1, coordinate_index) = dir.y;
						displacement_matrix.at(j * 3 + 2, coordinate_index) = dir.z;
					}
					else {
						displacement_matrix.at(j * 3    , coordinate_index) = 0.0;
						displacement_matrix.at(j * 3 + 1, coordinate_index) = 0.0;
						displacement_matrix.at(j * 3 + 2, coordinate_index) = 0.0;
					}
				}
				++coordinate_index;
			}
			if (i > 1) {
				vec3 axis = normalize(bond1);

				for (int j = 0; j < size(); ++j) {
					if (affected_indices[i * size() + j]) {
						vec3 dir = cross(atoms[j].position - atom2, axis);
						displacement_matrix.at(j * 3    , coordinate_index) = dir.x;
						displacement_matrix.at(j * 3 + 1, coordinate_index) = dir.y;
						displacement_matrix.at(j * 3 + 2, coordinate_index) = dir.z;
					}
					else {
						displacement_matrix.at(j * 3    , coordinate_index) = 0.0;
						displacement_matrix.at(j * 3 + 1, coordinate_index) = 0.0;
						displacement_matrix.at(j * 3 + 2, coordinate_index) = 0.0;
					}
				}
				++coordinate_index;
			}
		}
		return displacement_matrix;
	}

	void Molecule::printXYZData(const std::string& table_title) const {
		fout.resetRows();
		fout.offsetRight(2);
		fout.addRow(NumberFormat(), TextAlignment::centered, 45);
		fout << '|' << '_' << '-' << table_title << '|';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 8);
		fout.addRow(NumberFormat(), TextAlignment::centered, 10, TextPadding(2, 0, 0, 0));
		fout.addRow(NumberFormat(), TextAlignment::centered, 10, TextPadding(0, 0, 0, 0));
		fout.addRow(NumberFormat(), TextAlignment::centered, 10, TextPadding(0, 2, 0, 0));
		fout << '|' << '-' << '_' << "Atom" << '|' << '-' << '_' << "x [A]" << ',' << '-' << '_' << "y [A]" << ',' << '-' << '_' << "z [A]" << '|';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::left, 8);
		fout.addRow(NumberFormat::crudeFormat(10, 6), TextAlignment::right, 10, TextPadding(2, 0, 0, 0));
		fout.addRow(NumberFormat::crudeFormat(10, 6), TextAlignment::right, 10, TextPadding(0, 0, 0, 0));
		fout.addRow(NumberFormat::crudeFormat(10, 6), TextAlignment::right, 10, TextPadding(0, 2, 0, 0));

		for (int i = 0; i < atoms.size(); ++i) {
			const Atom& a = atoms[i];
			fout << '|' << getAtomName(i) << '|' << a.position.x * a0_to_A << ',' << a.position.y * a0_to_A << ',' << a.position.z * a0_to_A << '|' << '\n';
		}

		fout << '-' << ',' << '-' << ',' << '-' << ',' << '-' << '\n';
		fout.resetRows();
		fout << '\n';
	}

	void Molecule::printInternalCoordinates(const std::string& table_title) const {
		fout.resetRows();
		fout.offsetRight(2);
		fout.addRow(NumberFormat(), TextAlignment::centered, 37);
		fout << '|' << '_' << '-' << table_title << '|';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::centered, 20);
		fout.addRow(NumberFormat(), TextAlignment::centered, 12);
		fout << '|' << '-' << '_' << "Coordinate" << '|' << '-' << '_' << "Value" << '|';
		fout.resetRows();
		fout.addRow(NumberFormat(), TextAlignment::left, 20);
		fout.addRow(NumberFormat::crudeFormat(10, 6), TextAlignment::right, 12);

		VectorNd coordinates;
		calculateInternalCoordinates(coordinates);

		int coordinate_index = 0;
		for (int i = 0; i < internal_coordinate_indices.size() / 4; ++i) {
			uint index0 = internal_coordinate_indices[i * 4];
			uint index1 = internal_coordinate_indices[i * 4 + 1];
			uint index2 = internal_coordinate_indices[i * 4 + 2];
			uint index3 = internal_coordinate_indices[i * 4 + 3];

			std::string label = getAtomName(index2) + '-' + getAtomName(index3);
			fout << '|' << label << '|' << coordinates[coordinate_index] * a0_to_A << " A" << '|' << '\n';
			++coordinate_index;

			if (i) {
				std::string label = getAtomName(index1) + '-' + getAtomName(index2) + '-' + getAtomName(index3);
				fout << '|' << label << '|' << coordinates[coordinate_index] * 180.0 / PI << " d" << '|' << '\n';
				++coordinate_index;
			}
			if (i > 1) {
				std::string label = getAtomName(index0) + '-' + getAtomName(index1) + '-' + getAtomName(index2) + '-' + getAtomName(index3);
				fout << '|' << label << '|' << coordinates[coordinate_index] * 180.0 / PI << " d" << '|' << '\n';
				++coordinate_index;
			}
		}

		fout << '-' << ',' << '-' << '\n';
		fout.resetRows();
		fout << '\n';
	}
}