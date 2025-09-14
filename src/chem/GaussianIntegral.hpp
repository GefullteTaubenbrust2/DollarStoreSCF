#pragma once
#include "GTO.hpp"

namespace flo {
	double overlapIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu);

	double kineticEnergyIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu);

	double nuclearPotentialIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_nu, const vec3& nucleus);

	double electronRepulsionIntegral(const ContractedGaussian& chi_mu, const ContractedGaussian& chi_sigma, const ContractedGaussian& chi_nu, const ContractedGaussian& chi_tau);
}