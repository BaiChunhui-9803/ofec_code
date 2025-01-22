#ifndef OFEC_MMOP_METRICS_H
#define OFEC_MMOP_METRICS_H

#include "../../../core/problem/problem.h"

namespace ofec {
	class MetricsMMOP : virtual public Problem {
	public:
		Real variableAccuracy() const { return m_variable_accuracy; }
		Real objectiveAccuracy() const { return m_objective_accuracy; }
		Real variableNicheRadius() const { return m_variable_niche_radius; }
		void updateCandidates(const SolutionBase &s, std::list<std::unique_ptr<SolutionBase>> &candidates) const;
		bool isSolved(const std::list<std::unique_ptr<SolutionBase>> &candidates) const;
		std::vector<bool> optimaFound(const std::list<std::unique_ptr<SolutionBase>> &candidates) const;
		size_t numOptimaFound(const std::list<std::unique_ptr<SolutionBase>> &candidates) const;

	protected:
		Real m_variable_accuracy = -1;
		Real m_objective_accuracy = -1;
		Real m_variable_niche_radius = -1;

		void updateByObj(const SolutionBase &s, std::list<std::unique_ptr<SolutionBase>> &candidates) const;
		size_t numOptimaByObj(const std::list<std::unique_ptr<SolutionBase>> &candidates) const;
		std::vector<bool> optimaFoundByObj(const std::list<std::unique_ptr<SolutionBase>> &candidates) const;

		void updateByVar(const SolutionBase &s, std::list<std::unique_ptr<SolutionBase>> &candidates) const;
		size_t numOptimaByVar(const std::list<std::unique_ptr<SolutionBase>> &candidates) const;
		std::vector<bool> optimaFoundByVar(const std::list<std::unique_ptr<SolutionBase>> &candidates) const;
	};
}

#endif // !OFEC_MMOP_METRICS_H
