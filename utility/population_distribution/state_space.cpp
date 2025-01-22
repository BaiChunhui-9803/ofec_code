#include "state_space.h"

namespace ofec::population_distribution {
	StateSpace::StateSpace(size_t num_divisions, Problem *problem) :
		m_number_divisions(num_divisions), 
		m_problem(problem) {}

	void StateSpace::addReferencePoint(const SolutionBase &sol) {
		m_reference_points.emplace_back(m_problem->createSolution(sol));
	}

	size_t StateSpace::numberReferencePoints() const {
		return m_reference_points.size();
	}

	size_t StateSpace::numberStates() const {
		return pow(m_number_divisions, m_reference_points.size());
	}

	size_t StateSpace::stateToID(const std::vector<size_t> &state) const {
		if (state.size() != m_reference_points.size() || m_reference_points.empty()) {
			return 0;
		}
		else {
			size_t id = state.back();
			size_t multiple = 1;
			for (int i = state.size() - 2; i >= 0; --i) {
				multiple *= m_number_divisions;
				id += state[i] * multiple;
			}
			return id;
		}
	}

	size_t StateSpace::solutionsToStateID(const std::vector<const SolutionBase*> &sols) const {
		std::vector<size_t> state(m_reference_points.size());
		Real var_epsilon = 1e-8;
		Real obj_epsilon = 1e-4;
		Real var_eta = -log10(var_epsilon);
		Real obj_eta = -log10(obj_epsilon);
		std::vector<std::vector<Real>> distance(m_reference_points.size(), std::vector<Real>(sols.size()));
		for (size_t i = 0; i < m_reference_points.size(); ++i) {
			for (size_t j = 0; j < sols.size(); j++) {
				distance[i][j] = m_problem->normalizedVariableDistance(
					m_reference_points[i]->variableBase(), sols[j]->variableBase());
			}
		}
		for (size_t i = 0; i < m_reference_points.size(); ++i) {
			size_t idx_nearest = 0;
			for (size_t j = 1; j < sols.size(); j++) {
				if (distance[i][idx_nearest] > distance[i][j]) {
					idx_nearest = j;
				}
			}
			Real mapped_var_dist = std::max(0.0, log10(distance[i][idx_nearest]) + var_eta) / var_eta;
			Real obj_diff = m_reference_points[i]->objectiveDistance(*sols[idx_nearest]);
			Real mapped_obj_diff = 1.0 - 1.0 / (1.0 + std::max(0.0, log10(obj_diff) + obj_eta) / obj_eta);
			state[i] = mapped_var_dist * m_number_divisions * 0.9 + mapped_obj_diff * m_number_divisions * 0.1;
			//min = std::max(0.0, log(min) + 8) / 8;
			//state[i] = min * m_number_divisions;
			//Real mean = 0;
			//for (size_t j = 0; j < points.size(); j++) {
			//	mean += distance[i][j];
			//}
			//mean /= points.size();
			//mean = std::max(0.0, log(mean) + 8) / 8;
			////mean = std::max(0.0, log(mean) - log(epsilon)) / (-log(epsilon));
			////mean = 1 - log(std::max(epsilon, mean)) / log(epsilon);
			//state[i] = mean * m_number_divisions;
		}
		return stateToID(state);
	}
}
