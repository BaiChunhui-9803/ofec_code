#include "con_hc.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../core/environment/environment.h"

namespace ofec {
	void ConHC::initialize_(Environment *env) {
		BaseHC::initialize_(env);
		m_step_size.resize(env->problem()->numberVariables());
		for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
			auto &range = CAST_CONOP(env->problem())->range(j);
			m_step_size[j] = (range.second - range.first) / 100;
		}
	}

	void ConHC::pickNeighbour(Environment *env) {
		auto &cur = dynamic_cast<Solution<TypeVar>&>(*m_cur);
		auto &neighbour = dynamic_cast<Solution<TypeVar>&>(*m_neighbour);
		do {
			for (size_t j = 0; j < neighbour.variable().size(); j++)
				neighbour.variable()[j] = m_random->normal.nextNonStd(cur.variable()[j], m_step_size[j]);
		} while (CAST_CONOP(env->problem())->boundaryViolated(neighbour));
	}
}