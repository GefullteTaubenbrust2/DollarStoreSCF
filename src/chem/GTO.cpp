#include "GTO.hpp"
#include <sstream>

namespace flo {
	CartesianPolynomialTerm::CartesianPolynomialTerm(double weight, uint x, uint y, uint z) : 
	weight(weight), x(x), y(y), z(z) {}

	bool CartesianPolynomialTerm::operator&&(const CartesianPolynomialTerm& other) const {
		return x == other.x && y == other.y && z == other.z;
	}

	GTOPrimitive::GTOPrimitive(double zeta, double weight) :
	zeta(zeta), weight(weight) {}

	ContractedGaussian::ContractedGaussian(const std::vector<GTOPrimitive>& p, uint l, int m) : primitives(p), l (l) {
		int am = std::abs(m);

		spherical_harmonic.reserve((am - (am < 0 ? 1 : 0)) / 2 + 1);
		for (int i = 0, x = m <= 0 ? am : (am - 1); x >= 0; x -= 2, ++i) {
			spherical_harmonic.push_back(CartesianPolynomialTerm((i & 1 ? -1.0 : 1.0) * binomial((double)am, (double)x), x, am - x, 0));
		}

		double legendre_factor = doubleFactorial(2.0 * (double)am - 1.0);
		if (m < 0) {
			legendre_factor *= (double)factorial(l - am) / (double)factorial(l + am);
		}

		for (int i = 0; i < spherical_harmonic.size(); ++i) {
			spherical_harmonic[i].weight *= legendre_factor;
		}

		double y_normalization = std::sqrt((double)(2 * l + 1) * (double)factorial(l - m) / (4.0 * PI * (double)factorial(l + m)));
		if (m) y_normalization *= std::sqrt(2.0);

		if (l > am) {
			std::vector<CartesianPolynomialTerm> buffer1 = spherical_harmonic;
			std::vector<CartesianPolynomialTerm> buffer2;
			std::vector<CartesianPolynomialTerm>* current = &spherical_harmonic;
			std::vector<CartesianPolynomialTerm>* previous = &buffer1;
			std::vector<CartesianPolynomialTerm>* next = &buffer2;

			for (int i = 0; i < current->size(); ++i) {
				(*current)[i].weight *= (double)(2 * am + 1);
				(*current)[i].z += 1;
			}

			for (int i = am + 1; i < l; ++i) {
				next->resize(current->size() + 3 * previous->size());
				double scale_factor = (double)(i - am + 1);
				for (int j = 0; j < current->size(); ++j) {
					(*next)[j] = (*current)[j];
					(*next)[j].z += 1;
					(*next)[j].weight *= (double)(2 * i + 1) / scale_factor;
				}
				for (int j = 0, k = current->size(); j < previous->size(); ++j, k += 3) {
					(*next)[k] = (*previous)[j];
					(*next)[k + 1] = (*previous)[j];
					(*next)[k + 2] = (*previous)[j];
					(*next)[k].x += 2;
					(*next)[k + 1].y += 2;
					(*next)[k + 2].z += 2;
					(*next)[k].weight *= -(double)(i + am) / scale_factor;
					(*next)[k + 1].weight *= -(double)(i + am) / scale_factor;
					(*next)[k + 2].weight *= -(double)(i + am) / scale_factor;
				}

				for (int j = 0; j < next->size(); ++j) {
					for (int k = j + 1; k < next->size(); ++k) {
						if ((*next)[j] && (*next)[k]) {
							(*next)[j].weight += (*next)[k].weight;
							next->erase(next->begin() + k);
							--k;
						}
					}
				}

				std::vector<CartesianPolynomialTerm>* ptr = previous;
				previous = current;
				current = next;
				next = ptr;
			}

			spherical_harmonic = *current;
		}

		for (int i = 0; i < spherical_harmonic.size(); ++i) {
			spherical_harmonic[i].weight *= y_normalization;
		}

		for (int i = 0; i < primitives.size(); ++i) {
			double normalization_factor = (std::pow(2.0 * primitives[i].zeta, (double)l + 1.5) * std::pow(2.0, (double)l + 2)) / (std::sqrt(PI) * doubleFactorial(2.0 * (double)l + 1.0));
			primitives[i].weight *= std::sqrt(normalization_factor);
		}
	}

	ContractedGaussian::ContractedGaussian(const std::vector<GTOPrimitive>& primitives, uint l, int m, const vec3& c) : 
	ContractedGaussian(primitives, l, m) {
		center = c;
	}

	double ContractedGaussian::operator()(const vec3& position) const {
		return 0.0;
	}
	
	GTOShell::GTOShell(uint l, const std::vector<GTOPrimitive>& primitives) : l(l), primitives(primitives) {}

	std::string getOrbitalName(uint l, int m) {
		switch (l) {
		case 0:
			return "s";
		case 1: {
			const std::string p_strings[3] = { "px", "pz", "py" };
			return p_strings[m + 1];
		}
		case 2: {
			const std::string d_strings[5] = { "dxy", "dxz", "dz2", "dyz", "dx2-y2" };
			return d_strings[m + 2];
		}
		case 3: {
			const std::string f_strings[7] = { "fx3-3xy2", "fxyz", "fxz2", "fz3", "fyz2", "fx2z-y2z", "f3x2y-y3"};
			return f_strings[m + 3];
		}
		default:
			std::ostringstream oss;
			oss << "l=" << l << "m=" << m;
			return oss.str();
		}
	}
}