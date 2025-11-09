#include "LBFGS.hpp"

namespace flo {
	LBFGSSolver::LBFGSSolver(uint iterations, uint dimensionality) : 
	iterations(iterations) {
		previous_vectors.reserve(iterations * 2);
		alpha.reserve(iterations);

		q.resize(dimensionality);
		previous_position.resize(dimensionality);
	}

	void LBFGSSolver::computeInverseHessianProduct(VectorNd& vector) {
		uint vector_count = previous_vectors.size() / 2;

		q = vector;

		if (vector_count) {
			alpha.resize(vector_count);

			for (int i = vector_count - 1; i >= 0; --i) {
				VectorNd& delta_x = previous_vectors[i * 2];
				VectorNd& delta_g = previous_vectors[i * 2 + 1];

				alpha[i] = dot(delta_g, q) / dot(delta_x, delta_g);

				q -= alpha[i] * delta_x;
			}

			VectorNd& delta_x = previous_vectors[0];
			VectorNd& delta_g = previous_vectors[1];

			vector = dot(delta_g, delta_x) / dot(delta_x, delta_x) * q;

			for (int i = 0; i < vector_count; ++i) {
				VectorNd& delta_x = previous_vectors[i * 2];
				VectorNd& delta_g = previous_vectors[i * 2 + 1];

				vector += (alpha[i] - dot(delta_x, vector) / dot(delta_x, delta_g)) * delta_g;
			}
		}
	}

	void LBFGSSolver::doStep(VectorNd& gradient, VectorNd& position, VectorNd& step, double objective) {
		if (append_vectors) {
			if (previous_vectors.size() / 2 < iterations) {
				previous_vectors.resize(previous_vectors.size() + 2);
			}
			else {
				for (int i = 0; i < previous_vectors.size() / 2 - 1; ++i) {
					previous_vectors[i * 2]     = previous_vectors[i * 2 + 2];
					previous_vectors[i * 2 + 1] = previous_vectors[i * 2 + 3];
				}
			}

			previous_vectors[previous_vectors.size() - 2] = position - previous_position;
			previous_vectors[previous_vectors.size() - 1] = gradient - previous_gradient;

			append_vectors = false;
		}

		if (line_search) {
			double middle = 0.5 * (left + right);

			double previous_dot = dot(previous_gradient, descent_direction);

			if (objective - previous_objective < sufficient_decrease * middle * previous_dot) {
				if (dot(gradient, descent_direction) >= curvature_condition * previous_dot) {
					line_search = false;
					append_vectors = true;
					position = previous_position + (step = -middle * descent_direction);
					return;
				}
				step_length = middle;
				left = middle + line_search_step;
			}
			else {
				right = middle - line_search_step;
			}

			if (left >= right) {
				line_search = false;
				append_vectors = true;
				position = previous_position + (step = -step_length * descent_direction);
				return;
			}

			middle = 0.5 * (left + right);

			position = previous_position - middle * descent_direction;
		}
		else {
			descent_direction = gradient;

			computeInverseHessianProduct(descent_direction);

			step_length = base_step;

			left = 0.0;
			right = 1.0;

			previous_position = position;
			previous_gradient = gradient;
			previous_objective = objective;

			if (use_line_search) {
				position += step = -0.5 * descent_direction;
				line_search = true;
			}
			else {
				position += step = -1.0 * descent_direction;
				append_vectors = true;
			}
		}
	}
}