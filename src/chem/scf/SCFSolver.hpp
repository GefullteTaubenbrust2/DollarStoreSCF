#pragma once

namespace scf {
	enum class SpinTreatment {
		restricted = 0,
		unrestricted = 1,
	};

	enum class Spin {
		alpha = 0,
		beta = 1,
		up = 0,
		down = 1,
	};

	void useCoreGuess();

	void solveMOs();

	void computeCoulombTerms(Spin spin);

	void computeExchangeTerms(Spin spin);
}