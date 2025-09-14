#include "GaussianIntegral.hpp"

#include <iostream>
#include <cmath>

namespace flo {
	const double precision = 0.000000001;

	std::vector<double> float_buffer[5];

	double obaraSaikaOverlap(double S00, double eta, double A, double B, double P, uint i, uint j) {
		if (j > i) {
			double X = A;
			A = B;
			B = X;
		}
		uint ij = max((int)i - (int)j, (int)j - (int)i);
		double current = S00;
		double previous = 0.0;
		for (int k = 0; k < ij; ++k) {
			double next = (P - A) * current + 0.5 * k * previous / eta;
			previous = current;
			current = next;
		}
		if (!min(i, j)) return current;
		std::vector<double> S;
		uint s = min(i, j) + 1;
		S.resize(s * 3);
		S[0] = current;
		for (int k = 1; k < s; ++k) {
			S[k] = (P - A) * current + 0.5 * (k + ij - 1) * previous / eta;
			previous = current;
			current = S[k];
		}
		for (int k = 1; k < s; ++k) {
			S[k + s] = (P - B) * S[k] + 0.5 * ((k + ij) * S[k - 1]) / eta;
		}
		for (int l = 1; l < s - 1; ++l) {
			for (int k = l + 1; k < s; ++k) {
				S[k + 2 * s] = (P - B) * S[k + s] + 0.5 * ((k + ij) * S[k + s - 1] + l * S[k]) / eta;
			}
			for (int k = l + 1; k < s; ++k) {
				S[k] = S[k + s];
				S[k + s] = S[k + 2 * s];
			}
		}
		return S[2 * s - 1];
	}

	double overlapIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu) {
		vec3 A = chi_mu.center;
		vec3 B = chi_nu.center;
		double integral = 0.0;
		for (int i = 0; i < chi_mu.primitives.size(); ++i) {
			for (int j = 0; j < chi_nu.primitives.size(); ++j) {
				double alpha = chi_mu.primitives[i].zeta;
				double beta = chi_nu.primitives[j].zeta;
				double eta = alpha + beta;
				vec3 P = (alpha * A + beta * B) / eta;
				vec3 K = P * P * eta - alpha * A * A - beta * B * B;
				double primitive_factor = chi_mu.primitives[i].weight * chi_nu.primitives[j].weight;
				double gaussian_factor = std::exp(-alpha * beta / eta * dot((A - B), (A - B)));
				double S00 = std::sqrt(PI / eta);

				for (int k = 0; k < chi_mu.spherical_harmonic.size(); ++k) {
					for (int l = 0; l < chi_nu.spherical_harmonic.size(); ++l) {
						double Sijx = obaraSaikaOverlap(S00, eta, A.x, B.x, P.x, chi_mu.spherical_harmonic[k].x, chi_nu.spherical_harmonic[l].x);
						double Sijy = obaraSaikaOverlap(S00, eta, A.y, B.y, P.y, chi_mu.spherical_harmonic[k].y, chi_nu.spherical_harmonic[l].y);
						double Sijz = obaraSaikaOverlap(S00, eta, A.z, B.z, P.z, chi_mu.spherical_harmonic[k].z, chi_nu.spherical_harmonic[l].z);
						integral += gaussian_factor * chi_mu.spherical_harmonic[k].weight * chi_nu.spherical_harmonic[l].weight * primitive_factor * Sijx * Sijy * Sijz;
					}
				}
			}
		}
		return integral;
	}

	// Dij=2aDi+1j-iDi-1j
	// Di + 1j = 2aDi + 2j - (i + 1)Dij
	// Di - 1j = 2aDij - (i - 1)Di - 2j
	// Dij = 2a(2aDi + 2j - (i + 1)Dij) - i(2aDij - (i - 1)Di - 2j)
	// Dij = 4a ^ 2Di + 2j - (2a(i + 1) + i2a)Dij + i(i - 1)Di - 2j
	vec2 obaraSaikaKineticEnergy(double S00, double alpha, double beta, double A, double B, double P, uint i, uint j) {
		double eta = alpha + beta;
		if (j > i) {
			double X = A;
			A = B;
			B = X;
			double x = alpha;
			alpha = beta;
			beta = x;
			uint k = i;
			i = j;
			j = k;
		}
		uint ij = i - j;
		double current = S00;
		double previous = 0.0;
		for (int k = 0; k < (int)ij - 2; ++k) {
			double next = (P - A) * current + 0.5 * k * previous / eta;
			previous = current;
			current = next;
		}
		std::vector<double> S;
		uint sj = j + 1;
		uint si = i + 3 - max(0, (int)ij - 2);
		S.resize(si * 3);
		for (int k = 0; k < si * 3; ++k) S[k] = 0.0;
		int index_offset = max(0, (int)ij - 2);
		S[0] = current;
		for (int k = 1; k < si; ++k) {
			S[k] = (P - A) * current + 0.5 * (k + index_offset - 1) * previous / eta;
			previous = current;
			current = S[k];
		}
		double S1 = 0.0, S2, S3;
		if (sj > 1) {
			for (int k = 1; k < si; ++k) {
				S[k + si] = (P - B) * S[k] + 0.5 * ((k + index_offset) * S[k - 1]) / eta;
			}
			for (int l = 1; l < (int)sj - 1; ++l) {
				for (int k = l - 1; k < si; ++k) {
					S[k + 2 * si] = (P - B) * S[k + si] + 0.5 * ((k + index_offset) * S[k + si - 1] + l * S[k]) / eta;
				}
				for (int k = l - 1; k < si; ++k) {
					S[k] = S[k + si];
					S[k + si] = S[k + 2 * si];
				}
			}
			if (si > 4) S1 = S[si * 2 - 5];
			S2 = S[si * 2 - 3];
			S3 = S[si * 2 - 1];
		}
		else {
			if (si > 4) S1 = S[si - 5];
			S2 = S[si - 3];
			S3 = S[si - 1];
		}
		double Dij = 4.0 * alpha * alpha * S3 - 2.0 * alpha * (2 * i + 1) * S2 + (i * i - i) * S1;
		return vec2(Dij, S2);
	}

	double kineticEnergyIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu) {
		vec3 A = chi_mu.center;
		vec3 B = chi_nu.center;
		double integral = 0.0;
		for (int i = 0; i < chi_mu.primitives.size(); ++i) {
			for (int j = 0; j < chi_nu.primitives.size(); ++j) {
				double alpha = chi_mu.primitives[i].zeta;
				double beta = chi_nu.primitives[j].zeta;
				double eta = alpha + beta;
				vec3 P = (alpha * A + beta * B) / eta;
				double gaussian_factor = std::exp(-alpha * beta / eta * dot((A - B), (A - B)));
				double primitive_factor = chi_mu.primitives[i].weight * chi_nu.primitives[j].weight;
				double S00 = std::sqrt(PI / eta);

				for (int k = 0; k < chi_mu.spherical_harmonic.size(); ++k) {
					for (int l = 0; l < chi_nu.spherical_harmonic.size(); ++l) {
						vec2 DSijx = obaraSaikaKineticEnergy(S00, chi_mu.primitives[i].zeta, chi_nu.primitives[j].zeta, A.x, B.x, P.x, chi_mu.spherical_harmonic[k].x, chi_nu.spherical_harmonic[l].x);
						double& Dijx = DSijx.x;
						double& Sijx = DSijx.y;

						vec2 DSijy = obaraSaikaKineticEnergy(S00, chi_mu.primitives[i].zeta, chi_nu.primitives[j].zeta, A.y, B.y, P.y, chi_mu.spherical_harmonic[k].y, chi_nu.spherical_harmonic[l].y);
						double& Dijy = DSijy.x;
						double& Sijy = DSijy.y;
						
						vec2 DSijz = obaraSaikaKineticEnergy(S00, chi_mu.primitives[i].zeta, chi_nu.primitives[j].zeta, A.z, B.z, P.z, chi_mu.spherical_harmonic[k].z, chi_nu.spherical_harmonic[l].z);
						double& Dijz = DSijz.x;
						double& Sijz = DSijz.y;
						
						integral += gaussian_factor * chi_mu.spherical_harmonic[k].weight * chi_nu.spherical_harmonic[l].weight * primitive_factor * (Dijx * Sijy * Sijz + Sijx * Dijy * Sijz + Sijx * Sijy * Dijz);
					}
				}
			}
		}
		return -0.5 * integral;
	}

	double boysFunctionTaylor(double x, uint m) {
		double xm = 1.0;
		double F = 0.0;
		double err = 1.0;
		for (int i = 0; err > 0.000000001; i += 2) {
			F += xm / ((2 * (m + i) + 1) * factorial(i));
			xm *= x;
			err = xm / ((2 * (m + i) + 3) * factorial(i + 1));
			F -= err;
			xm *= x;
		}
		return F;
	}

	double boysFunctionExact(double x, uint m) {
		// Crude error estimate
		if (x > (double)m * 1.1 + 25.0) {
			return doubleFactorial(2 * (int)m - 1) / (double)pow(2, m + 1) * std::sqrt(PI / pow(x, 2 * m + 1));
		}
		else {
			double F = 1.0;
			double a = 1.0;
			double F_large = doubleFactorial(2 * (int)m - 1) / pow(2, m + 1) * std::sqrt(PI / pow(x, 2 * m + 1));
			for (int i = 1; a > 0.000000001; ++i) {
				a = 2.0 * x / (2 * m + 2 * i + 1) * a;
				F += a;
			}
			return std::exp(-x) / (double)(2 * m + 1) * F;
		}
	}
	
	void obaraSaikaNuclearPotential(std::vector<double>& buffer, double eta, double A, double B, double C, double P, uint i, uint j, uint N) {
		if (j > i) {
			double X = A;
			A = B;
			B = X;
			uint x = i;
			i = j;
			j = x;
		}

		uint ij = i - j;
		uint sj = j + 1;
		uint sN = N + i + j + 1;
		std::vector<double> buffer1(sj * sN), buffer2(sj * sN), buffer3(sj * sN);
		std::vector<double>* previous = &buffer1, *current = &buffer2, *next = &buffer3;
		for (int M = 0; M < sN; ++M) {
			(*current)[M] = buffer[M];
		}
		for (int k = 0; k < ij; ++k) {
			for (int M = sN - k - 2; M >= 0; --M) {
				(*next)[M] = (P - A) * (*current)[M] + 0.5 / eta * k * (*previous)[M] - (P - C) * (*current)[M + 1] - 0.5 / eta * k * (*previous)[M + 1];
			}
			std::vector<double>* temp = previous;
			previous = current;
			current = next;
			next = temp;
		}
		if (sj > 1) {
			for (int M = sN - ij - 2; M >= 0; --M) {
				(*current)[M + sN] = (P - A) * (*current)[M] + 0.5 / eta * ij * (*previous)[M] - (P - C) * (*current)[M + 1] - 0.5 / eta * ij * (*previous)[M + 1];
			}
			for (int k = 1; k < sj - 1; ++k) {
				for (int M = sN - ij - k - 2; M >= 0; --M) {
					(*current)[(k + 1) * sN + M] =
						  (P - A) * (*current)[k * sN + M]     + 0.5 / eta * (k + ij) * (*current)[(k - 1) * sN + M]
						- (P - C) * (*current)[k * sN + M + 1] - 0.5 / eta * (k + ij) * (*current)[(k - 1) * sN + M + 1];
				}
			}
			for (int k = 1; k < sj; ++k) {
				for (int M = sN - i - 2; M >= 0; --M) {
					(*next)[k * sN + M] =
						  (P - B) * (*current)[k * sN + M]     + 0.5 / eta * (k + ij) * (*current)[(k - 1) * sN + M]
						- (P - C) * (*current)[k * sN + M + 1] - 0.5 / eta * (k + ij) * (*current)[(k - 1) * sN + M + 1];
				}
			}
			std::vector<double>* temp = previous;
			previous = current;
			current = next;
			next = temp;
			for (int l = 1; l < sj - 1; ++l) {
				for (int k = l + 1; k < sj; ++k) {
					for (int M = sN - i - l - 2; M >= 0; --M) {
						(*next)[k * sN + M] =
							  (P - B) * (*current)[k * sN + M]     + 0.5 / eta * ((k + ij) * (*current)[(k - 1) * sN + M]     + l * (*previous)[k * sN + M])
							- (P - C) * (*current)[k * sN + M + 1] - 0.5 / eta * ((k + ij) * (*current)[(k - 1) * sN + M + 1] + l * (*previous)[k * sN + M + 1]);
					}
				}
				std::vector<double>* temp = previous;
				previous = current;
				current = next;
				next = temp;
			}
		}
		for (int M = 0; M <= N; ++M) {
			buffer[M] = (*current)[(sj - 1) * sN + M];
		}
	}

	double nuclearPotentialIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu, const vec3& nucleus) {
		vec3 A = chi_mu.center;
		vec3 B = chi_nu.center;
		double integral = 0.0;
		for (int k = 0; k < chi_mu.spherical_harmonic.size(); ++k) {
			for (int l = 0; l < chi_nu.spherical_harmonic.size(); ++l) {
				uint exponents_mu[3] = { chi_mu.spherical_harmonic[k].x, chi_mu.spherical_harmonic[k].y, chi_mu.spherical_harmonic[k].z };
				uint exponents_nu[3] = { chi_nu.spherical_harmonic[l].x, chi_nu.spherical_harmonic[l].y, chi_nu.spherical_harmonic[l].z };

				uint N1 = exponents_mu[2] + exponents_nu[2];
				uint N2 = N1 + exponents_mu[1] + exponents_nu[1];
				uint N3 = N2 + exponents_mu[0] + exponents_nu[0];

				std::vector<double> buffer(N3 + 1);

				double angular_weight = chi_mu.spherical_harmonic[k].weight * chi_nu.spherical_harmonic[l].weight;

				for (int i = 0; i < chi_mu.primitives.size(); ++i) {
					for (int j = 0; j < chi_nu.primitives.size(); ++j) {
						double alpha = chi_mu.primitives[i].zeta;
						double beta = chi_nu.primitives[j].zeta;
						double eta = alpha + beta;
						vec3 P = (alpha * A + beta * B) / eta;
						double gaussian_factor = std::exp(-alpha * beta / eta * dot((A - B), (A - B)));
						double primitive_factor = chi_mu.primitives[i].weight * chi_nu.primitives[j].weight;
						double RPC2 = dot(P - nucleus, P - nucleus);

						for (int m = 0; m < buffer.size(); ++m) buffer[m] = boysFunctionExact(eta * RPC2, m);

						obaraSaikaNuclearPotential(buffer, eta, A.x, B.x, nucleus.x, P.x, exponents_mu[0], exponents_nu[0], N2);
						obaraSaikaNuclearPotential(buffer, eta, A.y, B.y, nucleus.y, P.y, exponents_mu[1], exponents_nu[1], N1);
						obaraSaikaNuclearPotential(buffer, eta, A.z, B.z, nucleus.z, P.z, exponents_mu[2], exponents_nu[2], 0);
						integral += 2.0 * PI / eta * gaussian_factor * primitive_factor * angular_weight * buffer[0];
					}
				}
			}
		}
		return integral;
	}

	void obaraSaikaRepulsion(std::vector<double>& buffer, double A, double B, double C, double D, double P, double Q, double eta, double zeta, uint it, uint jt, uint kt, uint lt, uint Nt) {
		if (jt > it) {
			double X = A;
			A = B;
			B = X;
			uint x = it;
			it = jt;
			jt = x;
		}
		if (lt > kt) {
			double X = C;
			C = D;
			D = X;
			uint x = kt;
			kt = lt;
			lt = x;
		}
		if (lt > jt) {
			double X = A;
			A = C;
			C = X;

			X = B;
			B = D;
			D = X;

			X = eta;
			eta = zeta;
			zeta = X;

			X = P;
			P = Q;
			Q = X;

			uint x = it;
			it = kt;
			kt = x;

			x = jt;
			jt = lt;
			lt = x;
		}
		double omega = eta * zeta / (eta + zeta);

		uint si = it + 1;
		uint sj = jt + 1;
		uint sk = kt + 1;
		uint sl = lt + 1;
		uint sN = Nt + it + jt + kt + lt + 1;

		uint oi = sN;
		uint oj = oi * si;
		uint ok = oj * sj;
		uint size_total = ok * sk;

		float_buffer[0].resize(size_total);
		float_buffer[1].resize(size_total);
		float_buffer[2].resize(size_total);
		std::vector<double>* previous = &float_buffer[0], *current = &float_buffer[1], *next = &float_buffer[2];

		for (int N = 0; N < sN; ++N) {
			(*current)[N] = buffer[N];
		}
		if (it) {
			for (int N = sN - 2; N >= 0; --N) {
				(*current)[sN + N] = (P - A) * (*current)[N] - omega / eta * (P - Q) * (*current)[N + 1];
			}
		}
		for (int i = 1; i < it; ++i) {
			for (int N = sN - i - 2; N >= 0; --N) {
				(*current)[(i + 1) * sN + N] =
					(P - A) * (*current)[i * sN + N] - omega / eta * (P - Q) * (*current)[i * sN + N + 1]
					+ 0.5 * i / eta * ((*current)[(i - 1) * sN + N] - omega / eta * (*current)[(i - 1) * sN + N + 1]);
			}
		}

		if (jt) {
			for (int N = sN - it - 2; N >= 0; --N) {
				(*current)[oj + N] = (P - B) * (*current)[N] - omega / eta * (P - Q) * (*current)[N + 1];
			}
			for (int i = 1; i < si; ++i) {
				for (int N = sN - it - 2; N >= 0; --N) {
					(*current)[(si + i) * sN + N] =
						(P - B) * (*current)[i * sN + N] - omega / eta * (P - Q) * (*current)[i * sN + N + 1] 
						+ 0.5 * i / eta * ((*current)[(i - 1) * sN + N] - omega / eta * (*current)[(i - 1) * sN + N + 1]);
				}
			}
			for (int j = 1; j < jt; ++j) {
				for (int N = sN - it - j - 2; N >= 0; --N) {
					(*current)[((j + 1) * si) * sN + N] =
						(P - B) * (*current)[j * si * sN + N] - omega / eta * (P - Q) * (*current)[j * si * sN + N + 1]
						+ 0.5 * j / eta * ((*current)[((j - 1) * si) * sN + N] - omega / eta * (*current)[((j - 1) * si) * sN + N + 1]);
				}
				for (int i = 1; i < si; ++i) {
					for (int N = sN - it - j - 2; N >= 0; --N) {
						(*current)[((j + 1) * si + i) * sN + N] =
							(P - B) * (*current)[(j * si + i) * sN + N] - omega / eta * (P - Q) * (*current)[(j * si + i) * sN + N + 1]
							+ 0.5 * i / eta * ((*current)[(j * si + i - 1) * sN + N] - omega / eta * (*current)[(j * si + i - 1) * sN + N + 1])
							+ 0.5 * j / eta * ((*current)[((j - 1) * si + i) * sN + N] - omega / eta * (*current)[((j - 1) * si + i) * sN + N + 1]);
					}
				}
			}
		}

		if (kt) {
			for (int N = sN - it - jt - 2; N >= 0; --N) {
				(*current)[ok + N] = (Q - C) * (*current)[N] - omega / zeta * (Q - P) * (*current)[N + 1];
			}
			for (int i = 1; i < si; ++i) {
				for (int N = sN - it - jt - 2; N >= 0; --N) {
					int index = i * oi + N;
					(*current)[index + ok] =
						(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * i / (eta + zeta) * ((*current)[index - oi + 1]);
				}
			}
			for (int j = 1; j < sj; ++j) {
				for (int N = sN - it - jt - 2; N >= 0; --N) {
					int index = j * oj + N;
					(*current)[index + ok] =
						(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * j / (eta + zeta) * ((*current)[index - oj + 1]);
				}
			}
			for (int j = 1; j < sj; ++j) {
				for (int i = 1; i < si; ++i) {
					for (int N = sN - it - jt - 2; N >= 0; --N) {
						int index = j * oj + i * oi + N;
						(*current)[index + ok] =
							(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
							+ 0.5 * i / (eta + zeta) * ((*current)[index - oi + 1])
							+ 0.5 * j / (eta + zeta) * ((*current)[index - oj + 1]);
					}
				}
			}
			for (int k = 1; k < kt; ++k) {
				for (int N = sN - it - jt - k - 2; N >= 0; --N) {
					int index = k * ok + N;
					(*current)[index + ok] =
						(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1]);
				}
				for (int i = 1; i < si; ++i) {
					for (int N = sN - it - jt - k - 2; N >= 0; --N) {
						int index = k * ok + i * oi + N;
						(*current)[index + ok] =
							(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
							+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
							+ 0.5 * i / (eta + zeta) * (*current)[index - oi + 1];
					}
				}
				for (int j = 1; j < sj; ++j) {
					for (int N = sN - it - jt - k - 2; N >= 0; --N) {
						int index = k * ok + j * oj + N;
						(*current)[index + ok] =
							(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
							+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
							+ 0.5 * j / (eta + zeta) * (*current)[index - oj + 1];
					}
				}
				for (int j = 1; j < sj; ++j) {
					for (int i = 1; i < si; ++i) {
						for (int N = sN - it - jt - k - 2; N >= 0; --N) {
							int index = k * ok + j * oj + i * oi + N;
							(*current)[index + ok] =
								(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
								+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
								+ 0.5 * i / (eta + zeta) * (*current)[index - oi + 1]
								+ 0.5 * j / (eta + zeta) * (*current)[index - oj + 1];
						}
					}
				}
			}
		}

		if (lt) {
			for (int ki = 0; ki < lt; ++ki) {
				int k = ki + kt - lt + 1;
				for (int ji = lt - ki - 1; ji < lt; ++ji) {
					int j = ji + jt - lt + 1;
					for (int ii = lt - ji - 1; ii < lt; ++ii) {
						int i = ii + it - lt + 1;
						for (int N = sN - it - jt - kt - 2; N >= 0; --N) {
							int index = k * ok + j * oj + i * oi + N;
							(*next)[index] =
								(Q - D) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
								+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
								+ 0.5 * i / (eta + zeta) * (*current)[index - oi + 1]
								+ 0.5 * j / (eta + zeta) * (*current)[index - oj + 1];
						}
					}
				}
			}
			std::vector<double>* temp = next;
			next = previous;
			previous = current;
			current = temp;
		}
		for (int l = 1; l < lt; ++l) {
			for (int ki = l; ki < lt; ++ki) {
				int k = ki + kt - lt + 1;
				for (int ji = lt - ki - 1; ji < lt; ++ji) {
					int j = ji + jt - lt + 1;
					for (int ii = lt - ji - 1; ii < lt; ++ii) {
						int i = ii + it - lt + 1;
						for (int N = sN - it - jt - kt - l - 2; N >= 0; --N) {
							int index = k * ok + j * oj + i * oi + N;
							(*next)[index] =
								(Q - D) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
								+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
								+ 0.5 * l / zeta * ((*previous)[index] - omega / zeta * (*previous)[index + 1])
								+ 0.5 * i / (eta + zeta) * (*current)[index - oi + 1]
								+ 0.5 * j / (eta + zeta) * (*current)[index - oj + 1];
						}
					}
				}
			}
			std::vector<double>* temp = next;
			next = previous;
			previous = current;
			current = temp;
		}
		int index = (si * sj * sk - 1) * sN;
		for (int N = 0; N <= Nt; ++N) {
			buffer[N] = (*current)[index + N];
		}
	}

	double electronRepulsionIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_sigma, const ContractedGaussian& chi_nu, const ContractedGaussian& chi_tau) {
		vec3 A = chi_mu.center;
		vec3 B = chi_nu.center;
		vec3 C = chi_sigma.center;
		vec3 D = chi_tau.center;
		double integral = 0.0;

		float_buffer[3].resize(chi_mu.l + chi_nu.l + chi_sigma.l + chi_tau.l + 1);
		float_buffer[4].resize(float_buffer[3].size());

		uint min_angular_size = min(min(chi_mu.spherical_harmonic.size(), chi_nu.spherical_harmonic.size()), min(chi_sigma.spherical_harmonic.size(), chi_tau.spherical_harmonic.size()));

		for (uint ip = 0; ip < chi_mu.primitives.size(); ++ip) {
		for (uint jp = 0; jp < chi_nu.primitives.size(); ++jp) {
		for (uint kp = 0; kp < chi_sigma.primitives.size(); ++kp) {
		for (uint lp = 0; lp < chi_tau.primitives.size(); ++lp) {
			double alpha = chi_mu.primitives[ip].zeta;
			double beta  = chi_nu.primitives[jp].zeta;
			double gamma = chi_sigma.primitives[kp].zeta;
			double delta = chi_tau.primitives[lp].zeta;
			double eta = alpha + beta;
			double zeta = gamma + delta;
			vec3 P = (alpha * A + beta * B) / eta;
			vec3 Q = (gamma * C + delta * D) / zeta;
			double gaussian_factor_1 = std::exp(-alpha * beta / eta * dot((A - B), (A - B)));
			double gaussian_factor_2 = std::exp(-gamma * delta / zeta * dot((C - D), (C - D)));
			double RPQ2 = dot(P - Q, P - Q);

			double primitive_factor = 
				chi_mu.primitives[ip].weight * 
				chi_nu.primitives[jp].weight * 
				chi_sigma.primitives[kp].weight * 
				chi_tau.primitives[lp].weight;

			for (int m = 0; m < float_buffer[3].size(); ++m) float_buffer[3][m] = boysFunctionExact(eta * zeta * RPQ2 / (eta + zeta), m);

			double y = 0.0;
			for (uint is = 0; is < chi_mu.spherical_harmonic.size(); ++is) {
			for (uint js = 0; js < chi_nu.spherical_harmonic.size(); ++js) {
			for (uint ks = 0; ks < chi_sigma.spherical_harmonic.size(); ++ks) {
			for (uint ls = 0; ls < chi_tau.spherical_harmonic.size(); ++ls) {
				double angular_multiplicity = 1.0;

				if (is == min_angular_size && js == min_angular_size && ks == min_angular_size && ls == min_angular_size) {
					if (is < js || ks < ls || is + js * min_angular_size < ks + ls * min_angular_size) continue;
					if (is != js) angular_multiplicity *= 2.0;
					if (ks != ls) angular_multiplicity *= 2.0;
					if (is != ks || js != ls) angular_multiplicity *= 2.0;
				}

				uint exp_mu[3]    = { chi_mu.spherical_harmonic[is].x,    chi_mu.spherical_harmonic[is].y,    chi_mu.spherical_harmonic[is].z    };
				uint exp_nu[3]    = { chi_nu.spherical_harmonic[js].x,    chi_nu.spherical_harmonic[js].y,    chi_nu.spherical_harmonic[js].z    };
				uint exp_sigma[3] = { chi_sigma.spherical_harmonic[ks].x, chi_sigma.spherical_harmonic[ks].y, chi_sigma.spherical_harmonic[ks].z };
				uint exp_tau[3]   = { chi_tau.spherical_harmonic[ls].x,   chi_tau.spherical_harmonic[ls].y,   chi_tau.spherical_harmonic[ls].z   };

				uint N1 = exp_mu[2] + exp_nu[2] + exp_sigma[2] + exp_tau[2];
				uint N2 = N1 + exp_mu[1] + exp_nu[1] + exp_sigma[1] + exp_tau[1];
				uint N3 = N2 + exp_mu[0] + exp_nu[0] + exp_sigma[0] + exp_tau[0];

				double angular_weight = 
					angular_multiplicity * 
					chi_mu.spherical_harmonic[is].weight * 
					chi_nu.spherical_harmonic[js].weight * 
					chi_sigma.spherical_harmonic[ks].weight * 
					chi_tau.spherical_harmonic[ls].weight;

				for (int m = 0; m < float_buffer[3].size(); ++m) float_buffer[4][m] = float_buffer[3][m];

				obaraSaikaRepulsion(float_buffer[4], A.x, B.x, C.x, D.x, P.x, Q.x, eta, zeta, exp_mu[0], exp_nu[0], exp_sigma[0], exp_tau[0], N2);
				obaraSaikaRepulsion(float_buffer[4], A.y, B.y, C.y, D.y, P.y, Q.y, eta, zeta, exp_mu[1], exp_nu[1], exp_sigma[1], exp_tau[1], N1);
				obaraSaikaRepulsion(float_buffer[4], A.z, B.z, C.z, D.z, P.z, Q.z, eta, zeta, exp_mu[2], exp_nu[2], exp_sigma[2], exp_tau[2], 0);
				integral += 2.0 * std::pow(PI, 2.5) / (eta * zeta * std::sqrt(eta + zeta)) * gaussian_factor_1 * gaussian_factor_2 * primitive_factor * angular_weight * float_buffer[4][0];
			}
			}
			}
			}
		}
		}
		}
		}
		return integral;
	}
}