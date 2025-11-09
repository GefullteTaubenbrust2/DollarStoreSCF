#include "GaussianIntegral.hpp"
#include "../util/Pointer.hpp"

#include <iostream>
#include <cmath>

namespace flo {
	const double precision = 0.000000001;

	Array<double> float_buffer[7];

	template<typename T>
	void swap(T& a, T& b) {
		T x = a;
		a = b;
		b = x;
	}

	// We use the Obara-Saika scheme [1] to compute relevant GTO integrals.
	// There may be better ways of doing this, especially since my implementation is probably less than optimal.
	// However, the OS scheme offers the benefit of being conceptually simple and easy to extend to the relevant
	// gradients and Hessians. The simplest implementation of the OS scheme would be a literal recursion scheme,
	// but it is easy to realize that for high angular momenta, the depth of recursion and number of redundant
	// function calls quickly becomes very large. For this reason, this implementation uses an iterative algorithm.
	// This is quite verbose and scary looking, but should be much more efficient.
	// 
	// [1] S. Obara and A. Saika, "Efficient recursive computation of molecular integrals over Cartesian Gaussian functions", J. Chem. Phys. 84, 3963 (1986).

	double OSOverlap(double S00, double eta, double A, double B, double P, uint i, uint j) {
		if (j > i) {
			swap(A, B);
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
		Array<double>& S = float_buffer[0];
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
						double Sijx = OSOverlap(S00, eta, A.x, B.x, P.x, chi_mu.spherical_harmonic[k].x, chi_nu.spherical_harmonic[l].x);
						double Sijy = OSOverlap(S00, eta, A.y, B.y, P.y, chi_mu.spherical_harmonic[k].y, chi_nu.spherical_harmonic[l].y);
						double Sijz = OSOverlap(S00, eta, A.z, B.z, P.z, chi_mu.spherical_harmonic[k].z, chi_nu.spherical_harmonic[l].z);
						integral += gaussian_factor * chi_mu.spherical_harmonic[k].weight * chi_nu.spherical_harmonic[l].weight * primitive_factor * Sijx * Sijy * Sijz;
					}
				}
			}
		}
		return integral;
	}

	// Derivation from recurrence relations:
	// Tij = 2aDi+1j - iDi-1j
	// Tij = 2a(2aSi+2j - (i+1)Sij) - i(2aSij - (i-1)Si-2j)
	// Tij = 4aaSi+2j - (4ai + 2a)Sij + i(i-1)Si-2j
	vec2 OSKineticEnergy(double S00, double alpha, double beta, double A, double B, double P, uint i, uint j) {
		double eta = alpha + beta;
		if (j > i) {
			swap(A, B);
			swap(alpha, beta);
			swap(i, j);
		}
		uint ij = i - j;
		double current = S00;
		double previous = 0.0;
		for (int k = 0; k < (int)ij - 2; ++k) {
			double next = (P - A) * current + 0.5 * k * previous / eta;
			previous = current;
			current = next;
		}
		Array<double>& S = float_buffer[0];

		int index_offset = (int)ij - 2;
		int j_deficit = -index_offset;
		index_offset = max(0, index_offset);

		uint sj = j + 1;
		uint si = i + 3 - index_offset;

		S.resize(si * 3);

		S[0] = current;
		for (int k = 1; k < si; ++k) {
			S[k] = (P - A) * current + 0.5 * (k + index_offset - 1) * previous / eta;
			previous = current;
			current = S[k];
		}
		double S1 = 0.0, S2, S3;
		if (sj > 1) {
			if (j_deficit >= 0) S[si] = (P - B) * S[0];
			for (int k = 1; k < si; ++k) {
				S[k + si] = (P - B) * S[k] + 0.5 * ((k + index_offset) * S[k - 1]) / eta;
			}
			for (int l = 1; l < (int)sj - 1; ++l) {
				if (j_deficit - l > 0) S[2 * si] = (P - B) * S[si] + 0.5 * l * S[0] / eta;
				for (int k = l; k < si; ++k) {
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
						vec2 DSijx = OSKineticEnergy(S00, alpha, beta, A.x, B.x, P.x, chi_mu.spherical_harmonic[k].x, chi_nu.spherical_harmonic[l].x);
						double& Dijx = DSijx.x;
						double& Sijx = DSijx.y;

						vec2 DSijy = OSKineticEnergy(S00, alpha, beta, A.y, B.y, P.y, chi_mu.spherical_harmonic[k].y, chi_nu.spherical_harmonic[l].y);
						double& Dijy = DSijy.x;
						double& Sijy = DSijy.y;
						
						vec2 DSijz = OSKineticEnergy(S00, alpha, beta, A.z, B.z, P.z, chi_mu.spherical_harmonic[k].z, chi_nu.spherical_harmonic[l].z);
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
			for (int i = 1; a > 0.000000001; ++i) {
				a = 2.0 * x / (2 * m + 2 * i + 1) * a;
				F += a;
			}
			return std::exp(-x) / (double)(2 * m + 1) * F;
		}
	}
	
	void OSNuclearPotential(Array<double>& buffer, double eta, double A, double B, double C, double P, uint i, uint j, uint N) {
		if (j > i) {
			swap(A, B);
			swap(i, j);
		}

		uint ij = i - j;
		uint sj = j + 1;
		uint sN = N + i + j + 1;

		float_buffer[0].resize(sj * sN);
		float_buffer[1].resize(sj * sN);
		float_buffer[2].resize(sj * sN);

		Array<double>* previous = &float_buffer[0], *current = &float_buffer[1], *next = &float_buffer[2];
		for (int M = 0; M < sN; ++M) {
			(*current)[M] = buffer[M];
		}
		for (int k = 0; k < ij; ++k) {
			for (int M = sN - k - 2; M >= 0; --M) {
				(*next)[M] = (P - A) * (*current)[M] + 0.5 / eta * k * (*previous)[M] - (P - C) * (*current)[M + 1] - 0.5 / eta * k * (*previous)[M + 1];
			}
			Array<double>* temp = previous;
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
			Array<double>* temp = previous;
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
				Array<double>* temp = previous;
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

				Array<double>& buffer = float_buffer[3];
				buffer.resize(N3 + 1);

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

						OSNuclearPotential(buffer, eta, A.x, B.x, nucleus.x, P.x, exponents_mu[0], exponents_nu[0], N2);
						OSNuclearPotential(buffer, eta, A.y, B.y, nucleus.y, P.y, exponents_mu[1], exponents_nu[1], N1);
						OSNuclearPotential(buffer, eta, A.z, B.z, nucleus.z, P.z, exponents_mu[2], exponents_nu[2], 0);
						integral += 2.0 * PI / eta * gaussian_factor * primitive_factor * angular_weight * buffer[0];
					}
				}
			}
		}
		return integral;
	}

	void OSRepulsion(Array<double>& buffer, double A, double B, double C, double D, double P, double Q, double eta, double zeta, int it, int jt, int kt, int lt, int Nt) {
		if (jt > it) {
			swap(A, B);
			swap(it, jt);
		}
		if (lt > kt) {
			swap(C, D);
			swap(kt, lt);
		}
		if (lt > jt) {
			swap(A, C);
			swap(B, D);
			swap(eta, zeta);
			swap(P, Q);
			swap(it, kt);
			swap(jt, lt);
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
		Array<double>* previous = &float_buffer[0], *current = &float_buffer[1], *next = &float_buffer[2];

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
			int i_start = max(0, it - jt - kt - lt + 1);

			if (!i_start) {
				for (int N = sN - it - 2; N >= 0; --N) {
					(*current)[oj + N] = (P - B) * (*current)[N] - omega / eta * (P - Q) * (*current)[N + 1];
				}
			}
			for (int i = max(1, i_start); i < si; ++i) {
				for (int N = sN - it - 2; N >= 0; --N) {
					(*current)[(si + i) * sN + N] =
						(P - B) * (*current)[i * sN + N] - omega / eta * (P - Q) * (*current)[i * sN + N + 1] 
						+ 0.5 * i / eta * ((*current)[(i - 1) * sN + N] - omega / eta * (*current)[(i - 1) * sN + N + 1]);
				}
			}
			for (int j = 1; j < jt; ++j) {
				i_start = max(0, it - jt - kt - lt + j + 1);

				if (!i_start) {
					for (int N = sN - it - j - 2; N >= 0; --N) {
						(*current)[((j + 1) * si) * sN + N] =
							(P - B) * (*current)[j * si * sN + N] - omega / eta * (P - Q) * (*current)[j * si * sN + N + 1]
							+ 0.5 * j / eta * ((*current)[((j - 1) * si) * sN + N] - omega / eta * (*current)[((j - 1) * si) * sN + N + 1]);
					}
				}
				for (int i = max(1, i_start); i < si; ++i) {
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
			int j_start = max(0, jt - kt - lt + 1);
			int i_start = max(0, it - kt - lt + jt);

			if (!j_start) {
				if (!i_start) {
					for (int N = sN - it - jt - 2; N >= 0; --N) {
						(*current)[ok + N] = (Q - C) * (*current)[N] - omega / zeta * (Q - P) * (*current)[N + 1];
					}
				}
				for (int i = max(1, i_start); i < si; ++i) {
					for (int N = sN - it - jt - 2; N >= 0; --N) {
						int index = i * oi + N;
						(*current)[index + ok] =
							(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
							+ 0.5 * i / (eta + zeta) * ((*current)[index - oi + 1]);
					}
				}
			}
			for (int j = 1; j < sj; ++j) {
				i_start = max(0, it - kt - lt + jt - j);

				if (!i_start) {
					for (int N = sN - it - jt - 2; N >= 0; --N) {
						int index = j * oj + N;
						(*current)[index + ok] =
							(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
							+ 0.5 * j / (eta + zeta) * ((*current)[index - oj + 1]);
					}
				}
				for (int i = max(1, i_start); i < si; ++i) {
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
				j_start = max(0, jt - kt - lt + k + 1);
				i_start = max(0, it - kt - lt + k + 1);

				if (!j_start) {
					for (int N = sN - it - jt - k - 2; N >= 0; --N) {
						int index = k * ok + N;
						(*current)[index + ok] =
							(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
							+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1]);
					}
					for (int i = max(1, i_start); i < si; ++i) {
						for (int N = sN - it - jt - k - 2; N >= 0; --N) {
							int index = k * ok + i * oi + N;
							(*current)[index + ok] =
								(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
								+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
								+ 0.5 * i / (eta + zeta) * (*current)[index - oi + 1];
						}
					}
				}
				for (int j = max(1, j_start); j < sj; ++j) {
					i_start = max(0, it - kt - lt + jt + k - j);

					if (!i_start) {
						for (int N = sN - it - jt - k - 2; N >= 0; --N) {
							int index = k * ok + j * oj + N;
							(*current)[index + ok] =
								(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
								+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
								+ 0.5 * j / (eta + zeta) * (*current)[index - oj + 1];
						}
					}
					for (int i = max(1, i_start); i < si; ++i) {
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
			Array<double>* temp = next;
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
			Array<double>* temp = next;
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

				OSRepulsion(float_buffer[4], A.x, B.x, C.x, D.x, P.x, Q.x, eta, zeta, exp_mu[0], exp_nu[0], exp_sigma[0], exp_tau[0], N2);
				OSRepulsion(float_buffer[4], A.y, B.y, C.y, D.y, P.y, Q.y, eta, zeta, exp_mu[1], exp_nu[1], exp_sigma[1], exp_tau[1], N1);
				OSRepulsion(float_buffer[4], A.z, B.z, C.z, D.z, P.z, Q.z, eta, zeta, exp_mu[2], exp_nu[2], exp_sigma[2], exp_tau[2], 0);
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

	// ------------------------------------------------------------------------------------------------------------------------------------------- //

	vec2 OSOverlapGradient(double S00, double alpha, double beta, double A, double B, double P, uint i, uint j) {
		double eta = alpha + beta;
		double gradient_sign = 1.0;
		if (j > i) {
			swap(A, B);
			swap(alpha, beta);
			swap(i, j);

			gradient_sign = -1.0;
		}
		uint ij = i - j;
		double current = S00;
		double previous = 0.0;
		for (int k = 0; k < (int)ij - 1; ++k) {
			double next = (P - A) * current + 0.5 * k * previous / eta;
			previous = current;
			current = next;
		}
		Array<double>& S = float_buffer[0];

		int index_offset = (int)ij - 1;
		int j_deficit = -index_offset;
		index_offset = max(0, index_offset);

		uint sj = j + 1;
		uint si = i + 2 - index_offset;
		S.resize(si * 3);
		S[0] = current;
		for (int k = 1; k < si; ++k) {
			S[k] = (P - A) * current + 0.5 * (k + index_offset - 1) * previous / eta;
			previous = current;
			current = S[k];
		}
		double S1 = 0.0, S2, S3;
		if (sj > 1) {
			if (j_deficit > 0) S[si] = (P - B) * S[0];
			for (int k = 1; k < si; ++k) {
				S[k + si] = (P - B) * S[k] + 0.5 * ((k + index_offset) * S[k - 1]) / eta;
			}
			for (int l = 1; l < (int)sj - 1; ++l) {
				for (int k = l; k < si; ++k) {
					S[k + 2 * si] = (P - B) * S[k + si] + 0.5 * ((k + index_offset) * S[k + si - 1] + l * S[k]) / eta;
				}
				for (int k = l - 1; k < si; ++k) {
					S[k] = S[k + si];
					S[k + si] = S[k + 2 * si];
				}
			}
			if (si > 2) S1 = S[si * 2 - 3];
			S2 = S[si * 2 - 2];
			S3 = S[si * 2 - 1];
		}
		else {
			if (si > 2) S1 = S[si - 3];
			S2 = S[si - 2];
			S3 = S[si - 1];
		}
		double Dij = 2.0 * alpha * S3 - i * S1;
		return vec2(gradient_sign * Dij, S2);
	}

	vec3 overlapGradient(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu) {
		vec3 A = chi_mu.center;
		vec3 B = chi_nu.center;
		vec3 integral = 0.0;
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
						vec2 DSijx = OSOverlapGradient(S00, alpha, beta, A.x, B.x, P.x, chi_mu.spherical_harmonic[k].x, chi_nu.spherical_harmonic[l].x);
						double& Dijx = DSijx.x;
						double& Sijx = DSijx.y;
						vec2 DSijy = OSOverlapGradient(S00, alpha, beta, A.y, B.y, P.y, chi_mu.spherical_harmonic[k].y, chi_nu.spherical_harmonic[l].y);
						double& Dijy = DSijy.x;
						double& Sijy = DSijy.y;
						vec2 DSijz = OSOverlapGradient(S00, alpha, beta, A.z, B.z, P.z, chi_mu.spherical_harmonic[k].z, chi_nu.spherical_harmonic[l].z);
						double& Dijz = DSijz.x;
						double& Sijz = DSijz.y;

						integral += gaussian_factor * chi_mu.spherical_harmonic[k].weight * chi_nu.spherical_harmonic[l].weight * primitive_factor * vec3(Dijx * Sijy * Sijz, Sijx * Dijy * Sijz, Sijx * Sijy * Dijz);
					}
				}
			}
		}
		return integral;
	}

	// Derivation from recurrence relations:
	// Gij = 2aTi+1j - iTi-1j
	// Gij = 2a(2aDi+2j - (i+1)Dij) - i(2aDij - (i-1)Di-2j)
	// Gij = 4aaDi+2j - (4ai + 2a)Dij + i(i-1)Di-2j
	// Gij = 4aa(2aSi+3j - (i+2)Si+1j) - (4ai + 2a)(2aSi+1j - iSi-1j) + i(i-1)(2aSi-1j - (i-2)Si-3j)
	// Gij = 8aaaSi+3j - 12aa(1+i)Si+1j + 2a(ii+2a)Si-1j - i(i-1)(i-2)Si-3j
	std::array<double, 4> OSKineticGradient(double S00, double alpha, double beta, double A, double B, double P, uint i, uint j) {
		double eta = alpha + beta;
		double gradient_sign = 1.0;
		if (j > i) {
			swap(A, B);
			swap(alpha, beta);
			swap(i, j);

			gradient_sign = -1.0;
		}
		uint ij = i - j;
		double current = S00;
		double previous = 0.0;
		for (int k = 0; k < (int)ij - 3; ++k) {
			double next = (P - A) * current + 0.5 * k * previous / eta;
			previous = current;
			current = next;
		}
		Array<double>& S = float_buffer[0];

		int index_offset = (int)ij - 3;
		int j_deficit = -index_offset;
		index_offset = max(0, index_offset);

		uint sj = j + 1;
		uint si = i + 4 - index_offset;
		S.resize(si * 3);
		S[0] = current;
		for (int k = 1; k < si; ++k) {
			S[k] = (P - A) * current + 0.5 * (k + index_offset - 1) * previous / eta;
			previous = current;
			current = S[k];
		}
		double S1 = 0.0, S2 = 0.0, S3 = 0.0, S4, S5, S6, S7;
		if (sj > 1) {
			if (j_deficit > 0) S[si] = (P - B) * S[0];
			for (int k = 1; k < si; ++k) {
				S[k + si] = (P - B) * S[k] + 0.5 * ((k + index_offset) * S[k - 1]) / eta;
			}
			for (int l = 1; l < (int)sj - 1; ++l) {
				if (j_deficit - l > 0) S[2 * si] = (P - B) * S[si] + 0.5 * l * S[0] / eta;
				for (int k = l; k < si; ++k) {
					S[k + 2 * si] = (P - B) * S[k + si] + 0.5 * ((k + index_offset) * S[k + si - 1] + l * S[k]) / eta;
				}
				for (int k = l - 1; k < si; ++k) {
					S[k] = S[k + si];
					S[k + si] = S[k + 2 * si];
				}
			}
			if (si > 6) S1 = S[2 * si - 7];
			if (si > 5) S2 = S[2 * si - 6];
			if (si > 4) S3 = S[2 * si - 5];
			S4 = S[2 * si - 4];
			S5 = S[2 * si - 3];
			S6 = S[2 * si - 2];
			S7 = S[2 * si - 1];
		}
		else {
			if (si > 6) S1 = S[si - 7];
			if (si > 5) S2 = S[si - 6];
			if (si > 4) S3 = S[si - 5];
			S4 = S[si - 4];
			S5 = S[si - 3];
			S6 = S[si - 2];
			S7 = S[si - 1];
		}
		double Gij = 2.0 * alpha * (alpha * (4.0 * alpha * S7 - 6 * (1 + i) * S5) + 3 * i * i * S3) - i * (i - 1) * (i - 2) * S1;
		double Tij = 4.0 * alpha * alpha * S6 - 2.0 * alpha * (2 * i + 1) * S4 + (i * i - i) * S2;
		double Dij = 2.0 * alpha * S5 - i * S3;
		return std::array<double, 4>{gradient_sign * Gij, Tij, gradient_sign * Dij, S4};
	}

	vec3 kineticEnergyGradient(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu) {
		vec3 A = chi_mu.center;
		vec3 B = chi_nu.center;
		vec3 integral = 0.0;
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
						auto GDSijx = OSKineticGradient(S00, alpha, beta, A.x, B.x, P.x, chi_mu.spherical_harmonic[k].x, chi_nu.spherical_harmonic[l].x);
						double& Gijx = GDSijx[0];
						double& Tijx = GDSijx[1];
						double& Dijx = GDSijx[2];
						double& Sijx = GDSijx[3];
						auto GDSijy = OSKineticGradient(S00, alpha, beta, A.y, B.y, P.y, chi_mu.spherical_harmonic[k].y, chi_nu.spherical_harmonic[l].y);
						double& Gijy = GDSijy[0];
						double& Tijy = GDSijy[1];
						double& Dijy = GDSijy[2];
						double& Sijy = GDSijy[3];
						auto GDSijz = OSKineticGradient(S00, alpha, beta, A.z, B.z, P.z, chi_mu.spherical_harmonic[k].z, chi_nu.spherical_harmonic[l].z);
						double& Gijz = GDSijz[0];
						double& Tijz = GDSijz[1];
						double& Dijz = GDSijz[2];
						double& Sijz = GDSijz[3];

						integral += -0.5 * gaussian_factor * chi_mu.spherical_harmonic[k].weight * chi_nu.spherical_harmonic[l].weight * primitive_factor * vec3(
							Gijx * Sijy * Sijz + Dijx * Tijy * Sijz + Dijx * Sijy * Tijz,
							Tijx * Dijy * Sijz + Sijx * Gijy * Sijz + Sijx * Dijy * Tijz,
							Tijx * Sijy * Dijz + Sijx * Tijy * Dijz + Sijx * Sijy * Gijz
						);
					}
				}
			}
		}
		return integral;
	}

	vec2 OSNuclearGradient(Array<double>& buffer, double alpha, double beta, double A, double B, double C, double P, uint i, uint j) {
		double eta = alpha + beta;

		bool swapped = false;
		if (j > i) {
			swap(A, B);
			swap(alpha, beta);
			swap(i, j);

			swapped = true;
		}

		uint ij = i - j;
		uint sj = j + 2;
		uint sN = i + j + 2;

		int index_offset = ij;

		float_buffer[0].resize(sj * sN);
		float_buffer[1].resize(sj * sN);
		float_buffer[2].resize(sj * sN);

		double V01 = 0.0, V21 = 0.0, V10 = 0.0, V12 = 0.0;

		Array<double>* previous = &float_buffer[0], * current = &float_buffer[1], * next = &float_buffer[2];
		for (int M = 0; M < sN; ++M) {
			(*current)[M] = buffer[M];
		}
		// Increase i up to index_offset inclusive
		for (int k = 0; k < index_offset; ++k) {
			for (int M = sN - k - 2; M >= 0; --M) {
				(*next)[M] = (P - A) * (*current)[M] + 0.5 / eta * k * (*previous)[M] - (P - C) * (*current)[M + 1] - 0.5 / eta * k * (*previous)[M + 1];
			}
			Array<double>* temp = previous;
			previous = current;
			current = next;
			next = temp;
		}
		if (!j && i) V01 = (*previous)[0];

		// Increase j by one for i = index_offset
		for (int M = sN - i - 2; M >= 0; --M) {
			(*next)[M] = 
				  (P - B) * (*current)[M]     + 0.5 / eta * index_offset * (*previous)[M]
				- (P - C) * (*current)[M + 1] - 0.5 / eta * index_offset * (*previous)[M + 1];
		}
		// Increase i by one
		for (int M = sN - index_offset - 2; M >= 0; --M) {
			(*current)[M + sN] = 
				  (P - A) * (*current)[M]     + 0.5 / eta * index_offset * (*previous)[M] 
				- (P - C) * (*current)[M + 1] - 0.5 / eta * index_offset * (*previous)[M + 1];
		}
		// Increase i the rest of the way
		for (int k = 1; k < sj - 1; ++k) {
			for (int M = sN - index_offset - k - 2; M >= 0; --M) {
				(*current)[(k + 1) * sN + M] =
					  (P - A) * (*current)[k * sN + M]     + 0.5 / eta * (k + index_offset) * (*current)[(k - 1) * sN + M]
					- (P - C) * (*current)[k * sN + M + 1] - 0.5 / eta * (k + index_offset) * (*current)[(k - 1) * sN + M + 1];
			}
		}

		if (j == 1) {
			V10 = (*current)[(sj - 2) * sN];
		}

		// Increase j to one
		for (int k = 1; k < sj - 1; ++k) {
			for (int M = sN - i - 2; M >= 0; --M) {
				(*next)[k * sN + M] =
					  (P - B) * (*current)[k * sN + M]     + 0.5 / eta * (k + index_offset) * (*current)[(k - 1) * sN + M]
					- (P - C) * (*current)[k * sN + M + 1] - 0.5 / eta * (k + index_offset) * (*current)[(k - 1) * sN + M + 1];
			}
		}
		// Increase j to one for the final i value
		for (int M = sN - i - 3; M >= 0; --M) {
			(*next)[(sj - 1) * sN + M] =
				  (P - B) * (*current)[(sj - 1) * sN + M]     + 0.5 / eta * (i + 1) * (*current)[(sj - 2) * sN + M]
				- (P - C) * (*current)[(sj - 1) * sN + M + 1] - 0.5 / eta * (i + 1) * (*current)[(sj - 2) * sN + M + 1];
		}
		Array<double>* temp = previous;
		previous = current;
		current = next;
		next = temp;

		if (j == 2) {
			V10 = (*current)[(sj - 2) * sN];
		}

		// Increase j the rest of the way
		for (int l = 1; l < sj - 1; ++l) {
			for (int k = l; k < sj - 1; ++k) {
				for (int M = sN - i - l - 2; M >= 0; --M) {
					(*next)[k * sN + M] =
						  (P - B) * (*current)[k * sN + M]     + 0.5 / eta * ((k + index_offset) * (*current)[(k - 1) * sN + M]     + l * (*previous)[k * sN + M])
						- (P - C) * (*current)[k * sN + M + 1] - 0.5 / eta * ((k + index_offset) * (*current)[(k - 1) * sN + M + 1] + l * (*previous)[k * sN + M + 1]);
				}
			}
			if (l < sj - 2) {
				for (int M = sN - i - l - 3; M >= 0; --M) {
					(*next)[(sj - 1) * sN + M] =
						  (P - B) * (*current)[(sj - 1) * sN + M]     + 0.5 / eta * ((i + 1) * (*current)[(sj - 2) * sN + M]     + l * (*previous)[(sj - 1) * sN + M])
						- (P - C) * (*current)[(sj - 1) * sN + M + 1] - 0.5 / eta * ((i + 1) * (*current)[(sj - 2) * sN + M + 1] + l * (*previous)[(sj - 1) * sN + M + 1]);
				}
			}
			Array<double>* temp = previous;
			previous = current;
			current = next;
			next = temp;

			if (l == j - 2) {
				V10 = (*current)[(sj - 2) * sN];
			}
		}

		if (j) V01 = (*previous)[(sj - 3) * sN];
		V21 = (*previous)[(sj - 1) * sN];
		V12 = (*current)[(sj - 2) * sN];

		double V11 = (*previous)[(sj - 2) * sN];

		double Di = 2.0 * alpha * V21 - i * V01;
		double Dj = 2.0 * beta * V12 - j * V10;

		if (swapped) return vec2(Dj, Di);
		return vec2(Di, Dj);
	}

	std::array<vec3, 2> nuclearPotentialGradient(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu, const vec3& nucleus) {
		vec3 A = chi_mu.center;
		vec3 B = chi_nu.center;
		vec3 left_integral = 0.0;
		vec3 right_integral = 0.0;
		for (int k = 0; k < chi_mu.spherical_harmonic.size(); ++k) {
			for (int l = 0; l < chi_nu.spherical_harmonic.size(); ++l) {
				uint exponents_mu[3] = { chi_mu.spherical_harmonic[k].x, chi_mu.spherical_harmonic[k].y, chi_mu.spherical_harmonic[k].z };
				uint exponents_nu[3] = { chi_nu.spherical_harmonic[l].x, chi_nu.spherical_harmonic[l].y, chi_nu.spherical_harmonic[l].z };

				uint Nx = exponents_mu[0] + exponents_nu[0];
				uint Ny = exponents_mu[1] + exponents_nu[1];
				uint Nz = exponents_mu[2] + exponents_nu[2];

				uint N_total = Nx + Ny + Nz;

				Array<double>& dx_buffer = float_buffer[3];
				Array<double>& dy_buffer = float_buffer[4];
				Array<double>& dz_buffer = float_buffer[5];
				dx_buffer.resize(N_total + 2);
				dy_buffer.resize(N_total + 2);
				dz_buffer.resize(N_total + 2);

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

						for (int m = 0; m < dx_buffer.size(); ++m) dx_buffer[m] = boysFunctionExact(eta * RPC2, m);
						dy_buffer = dx_buffer;
						dz_buffer = dx_buffer;

						OSNuclearPotential(dx_buffer, eta, A.y, B.y, nucleus.y, P.y, exponents_mu[1], exponents_nu[1], Nz + Nx + 1);
						OSNuclearPotential(dx_buffer, eta, A.z, B.z, nucleus.z, P.z, exponents_mu[2], exponents_nu[2], Nx + 1);
						vec2 x_integral = OSNuclearGradient(dx_buffer, alpha, beta, A.x, B.x, nucleus.x, P.x, exponents_mu[0], exponents_nu[0]);

						OSNuclearPotential(dy_buffer, eta, A.z, B.z, nucleus.z, P.z, exponents_mu[2], exponents_nu[2], Nx + Ny + 1);
						OSNuclearPotential(dy_buffer, eta, A.x, B.x, nucleus.x, P.x, exponents_mu[0], exponents_nu[0], Ny + 1);
						vec2 y_integral = OSNuclearGradient(dy_buffer, alpha, beta, A.y, B.y, nucleus.y, P.y, exponents_mu[1], exponents_nu[1]);

						OSNuclearPotential(dz_buffer, eta, A.x, B.x, nucleus.x, P.x, exponents_mu[0], exponents_nu[0], Ny + Nz + 1);
						OSNuclearPotential(dz_buffer, eta, A.y, B.y, nucleus.y, P.y, exponents_mu[1], exponents_nu[1], Nz + 1);
						vec2 z_integral = OSNuclearGradient(dz_buffer, alpha, beta, A.z, B.z, nucleus.z, P.z, exponents_mu[2], exponents_nu[2]);

						left_integral  += 2.0 * PI / eta * gaussian_factor * primitive_factor * angular_weight * vec3(x_integral.x, y_integral.x, z_integral.x);
						right_integral += 2.0 * PI / eta * gaussian_factor * primitive_factor * angular_weight * vec3(x_integral.y, y_integral.y, z_integral.y);
					}
				}
			}
		}
		return std::array<vec3, 2>{left_integral, right_integral};
	}

	std::array<double, 4> OSrepulsionGradient(Array<double>& buffer, double A, double B, double C, double D, double P, double Q, double alpha, double beta, double gamma, double delta, int it, int jt, int kt, int lt) {
		int result_permutation[4] = { 0, 1, 2, 3 };

		if (jt > it) {
			swap(A, B);
			swap(alpha, beta);
			swap(it, jt);
			swap(result_permutation[0], result_permutation[1]);
		}
		if (lt > kt) {
			swap(C, D);
			swap(gamma, delta);
			swap(kt, lt);
			swap(result_permutation[2], result_permutation[3]);
		}
		if (lt > jt) {
			swap(A, C);
			swap(B, D);
			swap(alpha, gamma);
			swap(beta, delta);
			swap(it, kt);
			swap(jt, lt);
			swap(result_permutation[0], result_permutation[2]);
			swap(result_permutation[1], result_permutation[3]);
			swap(P, Q);
		}
		double eta = alpha + beta;
		double zeta = gamma + delta;

		double omega = eta * zeta / (eta + zeta);

		uint si = it + 2;
		uint sj = jt + 2;
		uint sk = kt + 2;
		uint sl = lt + 2;
		uint sN = it + jt + kt + lt + 2;

		uint oi = sN;
		uint oj = oi * si;
		uint ok = oj * sj;
		uint size_total = ok * sk;

		float_buffer[0].resize(size_total);
		float_buffer[1].resize(size_total);
		float_buffer[2].resize(size_total);
		Array<double>* previous = &float_buffer[0], * current = &float_buffer[1], * next = &float_buffer[2];

		double results[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		
		for (int N = 0; N < sN; ++N) {
			(*current)[N] = buffer[N];
		}
		// Increase i
		for (int N = sN - 2; N >= 0; --N) {
			(*current)[sN + N] = (P - A) * (*current)[N] - omega / eta * (P - Q) * (*current)[N + 1];
		}
		for (int i = 1; i <= it; ++i) {
			for (int N = sN - i - 2; N >= 0; --N) {
				(*current)[(i + 1) * sN + N] =
					(P - A) * (*current)[i * sN + N] - omega / eta * (P - Q) * (*current)[i * sN + N + 1]
					+ 0.5 * i / eta * ((*current)[(i - 1) * sN + N] - omega / eta * (*current)[(i - 1) * sN + N + 1]);
			}
		}

		// Increase j
		int i_start = max(0, it - jt - kt - lt);

		if (!i_start) {
			for (int N = sN - it - 2; N >= 0; --N) {
				(*current)[oj + N] = (P - B) * (*current)[N] - omega / eta * (P - Q) * (*current)[N + 1];
			}
		}
		for (int i = max(1, i_start); i < si; ++i) {
			for (int N = sN - it - 2; N >= 0; --N) {
				(*current)[(si + i) * sN + N] =
					(P - B) * (*current)[i * sN + N] - omega / eta * (P - Q) * (*current)[i * sN + N + 1]
					+ 0.5 * i / eta * ((*current)[(i - 1) * sN + N] - omega / eta * (*current)[(i - 1) * sN + N + 1]);
			}
		}
		for (int j = 1; j <= jt; ++j) {
			i_start = max(0, it - jt - kt - lt + j);

			if (!i_start) {
				for (int N = sN - it - j - 2; N >= 0; --N) {
					(*current)[((j + 1) * si) * sN + N] =
						(P - B) * (*current)[j * si * sN + N] - omega / eta * (P - Q) * (*current)[j * si * sN + N + 1]
						+ 0.5 * j / eta * ((*current)[((j - 1) * si) * sN + N] - omega / eta * (*current)[((j - 1) * si) * sN + N + 1]);
				}
			}
			for (int i = max(1, i_start); i < si; ++i) {
				for (int N = sN - it - j - 2; N >= 0; --N) {
					(*current)[((j + 1) * si + i) * sN + N] =
						(P - B) * (*current)[(j * si + i) * sN + N] - omega / eta * (P - Q) * (*current)[(j * si + i) * sN + N + 1]
						+ 0.5 * i / eta * ((*current)[(j * si + i - 1) * sN + N] - omega / eta * (*current)[(j * si + i - 1) * sN + N + 1])
						+ 0.5 * j / eta * ((*current)[((j - 1) * si + i) * sN + N] - omega / eta * (*current)[((j - 1) * si + i) * sN + N + 1]);
				}
			}
		}

		// Increase k
		int j_start = max(0, jt - kt - lt);
		i_start = max(0, it - kt - lt + jt - 1);

		if (!j_start) {
			if (!i_start) {
				for (int N = sN - it - jt - 2; N >= 0; --N) {
					(*current)[ok + N] = (Q - C) * (*current)[N] - omega / zeta * (Q - P) * (*current)[N + 1];
				}
			}
			for (int i = max(1, i_start); i < si; ++i) {
				for (int N = sN - it - jt - 2; N >= 0; --N) {
					int index = i * oi + N;
					(*current)[index + ok] =
						(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * i / (eta + zeta) * ((*current)[index - oi + 1]);
				}
			}
		}
		for (int j = 1; j < sj; ++j) {
			i_start = max(0, it - kt - lt + jt - j - 1);

			if (!i_start) {
				for (int N = sN - it - jt - 2; N >= 0; --N) {
					int index = j * oj + N;
					(*current)[index + ok] =
						(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * j / (eta + zeta) * ((*current)[index - oj + 1]);
				}
			}
			for (int i = max(1, i_start); i < si; ++i) {
				for (int N = sN - it - jt - 2; N >= 0; --N) {
					int index = j * oj + i * oi + N;
					(*current)[index + ok] =
						(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * i / (eta + zeta) * ((*current)[index - oi + 1])
						+ 0.5 * j / (eta + zeta) * ((*current)[index - oj + 1]);
				}
			}
		}
		for (int k = 1; k <= kt; ++k) {
			j_start = max(0, jt - kt - lt + k);
			i_start = max(0, it - kt - lt + k);

			if (!j_start) {
				for (int N = sN - it - jt - k - 2; N >= 0; --N) {
					int index = k * ok + N;
					(*current)[index + ok] =
						(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1]);
				}
				for (int i = max(1, i_start); i < si; ++i) {
					for (int N = sN - it - jt - k - 2; N >= 0; --N) {
						int index = k * ok + i * oi + N;
						(*current)[index + ok] =
							(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
							+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
							+ 0.5 * i / (eta + zeta) * (*current)[index - oi + 1];
					}
				}
			}
			for (int j = 1; j < sj; ++j) {
				i_start = max(0, it - kt - lt + jt + k - j - 1);

				if (!i_start) {
					for (int N = sN - it - jt - k - 2; N >= 0; --N) {
						int index = k * ok + j * oj + N;
						(*current)[index + ok] =
							(Q - C) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
							+ 0.5 * k / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
							+ 0.5 * j / (eta + zeta) * (*current)[index - oj + 1];
					}
				}
				for (int i = max(1, i_start); i < si; ++i) {
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

		// Increase l
		if (lt) {
			if (lt == it) {
				for (int N = sN - it - jt - kt - 2; N >= 0; --N) {
					int index = kt * ok + jt * oj + N;
					(*next)[index] =
						(Q - D) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * kt / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
						+ 0.5 * jt / (eta + zeta) * (*current)[index - oj + 1];
				}
			}
			if (lt == jt) {
				for (int N = sN - it - jt - kt - 2; N >= 0; --N) {
					int index = kt * ok + it * oi + N;
					(*next)[index] =
						(Q - D) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * kt / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1])
						+ 0.5 * it / (eta + zeta) * (*current)[index - oi + 1];
				}
			}
			if (lt == kt) {
				for (int N = sN - it - jt - kt - 2; N >= 0; --N) {
					int index = jt * oj + it * oi + N;
					(*next)[index] =
						(Q - D) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1]
						+ 0.5 * it / (eta + zeta) * (*current)[index - oi + 1]
						+ 0.5 * jt / (eta + zeta) * (*current)[index - oj + 1];
				}
			}
			for (int ki = 0; ki <= lt + 1; ++ki) {
				int k = ki + kt - lt;
				if (k <= 0) continue;
				for (int ji = lt - ki - 1; ji <= lt + 1; ++ji) {
					int j = ji + jt - lt;
					if (j <= 0) continue;
					for (int ii = lt - ji - 1; ii <= lt + 1; ++ii) {
						int i = ii + it - lt;
						if (i <= 0) continue;
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
			Array<double>* temp = next;
			next = previous;
			previous = current;
			current = temp;
		}

		for (int l = 1; l < lt; ++l) {
			for (int ki = l - 1; ki <= lt + 1; ++ki) {
				int k = ki + kt - lt;
				if (k <= 0) continue;
				for (int ji = lt - ki - 1; ji <= lt + 1; ++ji) {
					int j = ji + jt - lt;
					if (j <= 0) continue;
					for (int ii = lt - ji - 1; ii <= lt + 1; ++ii) {
						int i = ii + it - lt;
						if (i <= 0) continue;
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
			Array<double>* temp = next;
			next = previous;
			previous = current;
			current = temp;
		}
		int index = it * oi + jt * oj + kt * ok;

		results[1] = (*current)[index + oi];
		results[3] = (*current)[index + oj];
		results[5] = (*current)[index + ok];

		results[7] = (Q - D) * (*current)[index] - omega / zeta * (Q - P) * (*current)[index + 1];
		if (it) results[7] += 0.5 * it / (eta + zeta) * (*current)[index - oi + 1];
		if (jt) results[7] += 0.5 * jt / (eta + zeta) * (*current)[index - oj + 1];
		if (kt) results[7] += 0.5 * kt / zeta * ((*current)[index - ok] - omega / zeta * (*current)[index - ok + 1]);
		if (lt) results[7] += 0.5 * lt / zeta * ((*previous)[index] - omega / zeta * (*previous)[index + 1]);

		if (it) results[0] = (*current)[index - oi];
		if (jt) results[2] = (*current)[index - oj];
		if (kt) results[4] = (*current)[index - ok];

		if (lt) {
			results[6] = (*previous)[it * oi + jt * oj + kt * ok];
		}

		std::array<double, 4> return_value;

		return_value[result_permutation[0]] = 2.0 * alpha * results[1] - it * results[0];
		return_value[result_permutation[1]] = 2.0 * beta  * results[3] - jt * results[2];
		return_value[result_permutation[2]] = 2.0 * gamma * results[5] - kt * results[4];
		return_value[result_permutation[3]] = 2.0 * delta * results[7] - lt * results[6];

		return return_value;
	}

	std::array<vec3, 4> electronRepulsionGradient(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_sigma, const ContractedGaussian& chi_nu, const ContractedGaussian& chi_tau) {
		vec3 A = chi_mu.center;
		vec3 B = chi_nu.center;
		vec3 C = chi_sigma.center;
		vec3 D = chi_tau.center;
		vec3 mu_integral = 0.0;
		vec3 nu_integral = 0.0;
		vec3 sigma_integral = 0.0;
		vec3 tau_integral = 0.0;

		uint N_total = chi_mu.l + chi_nu.l + chi_sigma.l + chi_tau.l;

		Array<double>& base_buffer = float_buffer[3];
		Array<double>& dx_buffer   = float_buffer[4];
		Array<double>& dy_buffer   = float_buffer[5];
		Array<double>& dz_buffer   = float_buffer[6];
		base_buffer.resize(N_total + 2);
		dx_buffer.resize(N_total + 2);
		dy_buffer.resize(N_total + 2);
		dz_buffer.resize(N_total + 2);

		uint min_angular_size = min(min(chi_mu.spherical_harmonic.size(), chi_nu.spherical_harmonic.size()), min(chi_sigma.spherical_harmonic.size(), chi_tau.spherical_harmonic.size()));

		for (uint ip = 0; ip < chi_mu.primitives.size(); ++ip) {
			for (uint jp = 0; jp < chi_nu.primitives.size(); ++jp) {
				for (uint kp = 0; kp < chi_sigma.primitives.size(); ++kp) {
					for (uint lp = 0; lp < chi_tau.primitives.size(); ++lp) {
						double alpha = chi_mu.primitives[ip].zeta;
						double beta = chi_nu.primitives[jp].zeta;
						double gamma = chi_sigma.primitives[kp].zeta;
						double delta = chi_tau.primitives[lp].zeta;
						double eta = alpha + beta;
						double zeta = gamma + delta;
						vec3 P = (alpha * A + beta * B) / eta;
						vec3 Q = (gamma * C + delta * D) / zeta;
						double gaussian_factor_1 = std::exp(-alpha * beta / eta * dot((A - B), (A - B)));
						double gaussian_factor_2 = std::exp(-gamma * delta / zeta * dot((C - D), (C - D)));
						double gaussian_factor = gaussian_factor_1 * gaussian_factor_2;
						double RPQ2 = dot(P - Q, P - Q);
						double normalization_factor = 2.0 * std::pow(PI, 2.5) / (eta * zeta * std::sqrt(eta + zeta));

						double primitive_factor =
							chi_mu.primitives[ip].weight *
							chi_nu.primitives[jp].weight *
							chi_sigma.primitives[kp].weight *
							chi_tau.primitives[lp].weight;

						for (int m = 0; m < dx_buffer.size(); ++m) base_buffer[m] = boysFunctionExact(eta * zeta * RPQ2 / (eta + zeta), m);

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

										uint exp_mu[3] =    { chi_mu.spherical_harmonic[is].x,    chi_mu.spherical_harmonic[is].y,    chi_mu.spherical_harmonic[is].z };
										uint exp_nu[3] =    { chi_nu.spherical_harmonic[js].x,    chi_nu.spherical_harmonic[js].y,    chi_nu.spherical_harmonic[js].z };
										uint exp_sigma[3] = { chi_sigma.spherical_harmonic[ks].x, chi_sigma.spherical_harmonic[ks].y, chi_sigma.spherical_harmonic[ks].z };
										uint exp_tau[3] =   { chi_tau.spherical_harmonic[ls].x,   chi_tau.spherical_harmonic[ls].y,   chi_tau.spherical_harmonic[ls].z };

										uint Nx = exp_mu[0] + exp_nu[0] + exp_sigma[0] + exp_tau[0];
										uint Ny = exp_mu[1] + exp_nu[1] + exp_sigma[1] + exp_tau[1];
										uint Nz = exp_mu[2] + exp_nu[2] + exp_sigma[2] + exp_tau[2];

										double angular_weight =
											angular_multiplicity *
											chi_mu.spherical_harmonic[is].weight *
											chi_nu.spherical_harmonic[js].weight *
											chi_sigma.spherical_harmonic[ks].weight *
											chi_tau.spherical_harmonic[ls].weight;

										dx_buffer = base_buffer;
										dy_buffer = base_buffer;
										dz_buffer = base_buffer;

										double gross_factor = normalization_factor * gaussian_factor * primitive_factor * angular_weight;

										OSRepulsion(dx_buffer, A.y, B.y, C.y, D.y, P.y, Q.y, eta, zeta, exp_mu[1], exp_nu[1], exp_sigma[1], exp_tau[1], Nz + Nx + 1);
										OSRepulsion(dx_buffer, A.z, B.z, C.z, D.z, P.z, Q.z, eta, zeta, exp_mu[2], exp_nu[2], exp_sigma[2], exp_tau[2], Nx + 1);
										auto Dx = OSrepulsionGradient(dx_buffer, A.x, B.x, C.x, D.x, P.x, Q.x, alpha, beta, gamma, delta, exp_mu[0], exp_nu[0], exp_sigma[0], exp_tau[0]);

										OSRepulsion(dy_buffer, A.z, B.z, C.z, D.z, P.z, Q.z, eta, zeta, exp_mu[2], exp_nu[2], exp_sigma[2], exp_tau[2], Nx + Ny + 1);
										OSRepulsion(dy_buffer, A.x, B.x, C.x, D.x, P.x, Q.x, eta, zeta, exp_mu[0], exp_nu[0], exp_sigma[0], exp_tau[0], Ny + 1);
										auto Dy = OSrepulsionGradient(dy_buffer, A.y, B.y, C.y, D.y, P.y, Q.y, alpha, beta, gamma, delta, exp_mu[1], exp_nu[1], exp_sigma[1], exp_tau[1]);

										OSRepulsion(dz_buffer, A.x, B.x, C.x, D.x, P.x, Q.x, eta, zeta, exp_mu[0], exp_nu[0], exp_sigma[0], exp_tau[0], Ny + Nz + 1);
										OSRepulsion(dz_buffer, A.y, B.y, C.y, D.y, P.y, Q.y, eta, zeta, exp_mu[1], exp_nu[1], exp_sigma[1], exp_tau[1], Nz + 1);
										auto Dz = OSrepulsionGradient(dz_buffer, A.z, B.z, C.z, D.z, P.z, Q.z, alpha, beta, gamma, delta, exp_mu[2], exp_nu[2], exp_sigma[2], exp_tau[2]);

										mu_integral    += gross_factor * vec3(Dx[0], Dy[0], Dz[0]);
										nu_integral    += gross_factor * vec3(Dx[1], Dy[1], Dz[1]);
										sigma_integral += gross_factor * vec3(Dx[2], Dy[2], Dz[2]);
										tau_integral   += gross_factor * vec3(Dx[3], Dy[3], Dz[3]);
									}
								}
							}
						}
					}
				}
			}
		}
		return std::array<vec3, 4>{mu_integral, sigma_integral, nu_integral, tau_integral};
	}
}