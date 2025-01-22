#include "init_pop_selected.h"

namespace ofec {
	void InitPopSelectedBase::initializeVariables(VariableBase &x, Random *rnd) const {
		auto selected = rnd->uniform.nextElem(m_sols_init_pop.begin(), m_sols_init_pop.end());
		x = (*selected)->variableBase();
	}

	void InitPopSelectedBase::addSolInitPop(const SolutionBase &s) {
		m_sols_init_pop.emplace_back(createSolution(s));
	}

	void InitPopSelectedBase::clearSolsInitPop() {
		m_sols_init_pop.clear();
	}
}