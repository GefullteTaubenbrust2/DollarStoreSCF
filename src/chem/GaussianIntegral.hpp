#pragma once
#include <array>
#include "GTO.hpp"

namespace flo {
	double overlapIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu);

	double kineticEnergyIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu);

	double nuclearPotentialIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu, const vec3& nucleus);

	double electronRepulsionIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_sigma, const ContractedGaussian& chi_nu, const ContractedGaussian& chi_tau);

	vec3 overlapGradient(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu);

	vec3 kineticEnergyGradient(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu);

	std::array<vec3, 2> nuclearPotentialGradient(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu, const vec3& nucleus);

	std::array<vec3, 4> electronRepulsionGradient(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_sigma, const ContractedGaussian& chi_nu, const ContractedGaussian& chi_tau);
}