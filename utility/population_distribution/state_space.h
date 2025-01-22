#ifndef OFEC_POPULATION_DISTRIBUTION_STATE_SPACE_H
#define OFEC_POPULATION_DISTRIBUTION_STATE_SPACE_H

#include "../../core/problem/problem.h"
#include "../../core/algorithm/population.h"

namespace ofec::population_distribution {
	class StateSpace {
	protected:
		size_t m_number_divisions;
		std::vector<std::unique_ptr<SolutionBase>> m_reference_points;
		Problem* m_problem;

		size_t stateToID(const std::vector<size_t> &state) const;

	public:
		StateSpace(size_t number_divisions, Problem *problem);
		void addReferencePoint(const SolutionBase &sol);
		size_t numberReferencePoints() const;
		size_t numberStates() const;
		size_t solutionsToStateID(const std::vector<const SolutionBase*> &sols) const;
	};
}

#endif // !OFEC_POPULATION_DISTRIBUTION_STATE_SPACE_H
