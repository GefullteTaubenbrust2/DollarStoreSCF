#pragma once
#include "Lalib.hpp"

namespace flo {
	// Limited memory Broyden-Fletcher-Goldfarb-Shanno algorithm
	//
	// C. G. Broyden, The Convergence of a Class of Double-rank Minimization Algorithms 1. General Considerations, IMA Journal of Applied Mathematics, Volume 6, Issue 1, March 1970, Pages 76-90.
	// R. Fletcher, A new approach to variable metric algorithms, The Computer Journal, Volume 13, Issue 3, 1970, Pages 317-322.
	// D. Goldfarb, A family of variable-metric methods derived by variational means, Mathematics of Computation, Volume 24, Issue 109, 1970, Pages 23-26.
	// D. F. Shanno, Optimal Conditioning of Quasi-Newton Methods, Mathematics of Computation, Volume 24, Issue 111, 1970, Pages 657-657.
	struct LBFGSSolver {
	private:
		uint iterations = 8;

		std::vector<VectorNd> previous_vectors;
		VectorNd previous_position;
		VectorNd previous_gradient;
		double previous_objective;

		bool line_search = false;
		bool append_vectors = false;

		VectorNd q;
		VectorNd descent_direction;

		double left = 0.0, right = 0.0;

		double step_length = 0.0;

		std::vector<double> alpha;

	public:
		double base_step = 0.0001;
		double line_search_step = 0.01;
		double sufficient_decrease = 0.0001;
		double curvature_condition = 0.9;

		bool use_line_search = false;

	public:
		LBFGSSolver(uint iterations, uint dimensionality);

		void computeInverseHessianProduct(VectorNd& vector);

		void doStep(VectorNd& gradient, VectorNd& position, VectorNd& step, double objective);
	};
}