#include "metrics_mmop.h"
#include "../../../utility/functional.h"

namespace ofec {
	void MetricsMMOP::updateCandidates(const SolutionBase &s, std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		if (m_variable_accuracy != -1) {
			updateByVar(s, candidates);
		}
		else if (m_objective_accuracy != -1 && m_variable_niche_radius != -1) {
			updateByObj(s, candidates);
		}	
	}

	size_t MetricsMMOP::numOptimaFound(const std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		if (m_variable_accuracy != -1) {
			return numOptimaByVar(candidates);
		}
		else if (m_objective_accuracy != -1 && m_variable_niche_radius != -1) {
			return numOptimaByObj(candidates);
		}	
		else {
			return 0;
		}
	}

	bool MetricsMMOP::isSolved(const std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		if (m_optima) {
			return numOptimaFound(candidates) == numberOptimumObjectives();
		}
		else {
			return false;
		}
	}

	std::vector<bool> MetricsMMOP::optimaFound(const std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		if (m_variable_accuracy != -1) {
			return optimaFoundByVar(candidates);
		}
		else if (m_objective_accuracy != -1 && m_variable_niche_radius != -1) {
			return optimaFoundByObj(candidates);
		}
		else {
			return std::vector<bool>(numberOptimumObjectives(), false);
		}
	}

	void MetricsMMOP::updateByVar(const SolutionBase &s, std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		if (candidates.empty()) {
			for (size_t i = 0; i < m_optima->numberSolutions(); i++) {
				candidates.emplace_back(createSolution(s));
			}
		}
		else {
			auto it = candidates.begin();
			for (size_t i = 0; i < m_optima->numberSolutions(); ++i, ++it) {
				if (variableDistance(s.variableBase(), m_optima->solutionBase(i).variableBase())
					< variableDistance((*it)->variableBase(), m_optima->solutionBase(i).variableBase())
				) {
					(*it).reset(createSolution(s));
				}
			}
		}
	}

	size_t MetricsMMOP::numOptimaByVar(const std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		size_t count = 0;
		for (size_t i = 0; i < m_optima->numberSolutions(); ++i) {
			for (auto &c : candidates) {
				if (variableDistance(c->variableBase(), m_optima->solutionBase(i).variableBase()) < m_variable_accuracy) {
					count++;
					break;
				}
			}
		}
		return count;
	}

	std::vector<bool> MetricsMMOP::optimaFoundByVar(const std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		std::vector<bool> result(m_optima->numberSolutions(), false);
		for (size_t i = 0; i < m_optima->numberSolutions(); ++i) {
			for (auto &c : candidates) {
				if (variableDistance(c->variableBase(), m_optima->solutionBase(i).variableBase()) < m_variable_accuracy) {
					result[i] = true;
					break;
				}
			}
		}
		return result;
	}

	void MetricsMMOP::updateByObj(const SolutionBase &s, std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		if (s.objectiveDistance(m_optima->objective(0)) < m_objective_accuracy) {
			for (auto &c : candidates) {
				if (variableDistance(c->variableBase(), s.variableBase()) < m_variable_niche_radius) {
					if (dominate(s, *c, m_optimize_mode)) {
						c.reset(createSolution(s));
					}
					return;
				}
			}
			candidates.emplace_back(createSolution(s));
		}
	}

	size_t MetricsMMOP::numOptimaByObj(const std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		std::list<SolutionBase*> opts_fnd;
		for (auto &c : candidates) {
			if (c->objectiveDistance(m_optima->objective(0)) < m_objective_accuracy) {
				bool is_new_opt = true;
				for (auto of : opts_fnd) {
					if (variableDistance(of->variableBase(), c->variableBase()) < m_variable_niche_radius) {
						is_new_opt = false;
						break;
					}
				}
				if (is_new_opt) {
					opts_fnd.push_back(c.get());
				}
			}
		}
		return opts_fnd.size();
	}

	std::vector<bool> MetricsMMOP::optimaFoundByObj(const std::list<std::unique_ptr<SolutionBase>> &candidates) const {
		std::vector<bool> result(m_optima->numberObjectives(), false);
		for (size_t i = 0; i < m_optima->numberSolutions(); ++i) {
			for (auto &c : candidates) {
				if (c->objectiveDistance(m_optima->objective(0)) < m_objective_accuracy) {
					if (variableDistance(c->variableBase(), m_optima->solutionBase(i).variableBase()) < m_variable_niche_radius) {
						result[i] = true;
						break;
					}
				}
			}
		}
		return result;
	}
}