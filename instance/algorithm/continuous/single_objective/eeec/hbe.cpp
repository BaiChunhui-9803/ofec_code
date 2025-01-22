#include "hbe.h"
#include "../../../../../../utility/clustering/nbc.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../core/environment/environment.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void HBE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_num_increase_each_iter, 5, 1000, 20));
	}

	void HBE::run_(Environment *env) {
		std::list<std::unique_ptr<SolutionBase>> sols;
		std::vector<const SolutionBase*> sol_ptrs;
		std::vector<size_t> ids_centers;
		auto optima = CAST_CONOP(env->problem())->optima();
		for (size_t i = 0; i < optima->numberSolutions(); ++i) {
			sols.emplace_back(new Solution<>(optima->solution(i)));
			sols.back()->evaluate(env);
			sol_ptrs.emplace_back(sols.back().get());
			ids_centers.emplace_back(i);
		}
		NBC nbc;
		while (!terminating()) {
			for (size_t i = 0; i < m_num_increase_each_iter; ++i) {
				sols.emplace_back(CAST_CONOP(env->problem())->createSolution());
				sols.back()->initialize(env, m_random.get());
				sols.back()->evaluate(env);
				sol_ptrs.emplace_back(sols.back().get());
			}
			nbc.setData(sol_ptrs);
			nbc.updateNbDistByKDTree(env);
			nbc.cutEdgesInGraph(ids_centers);
			nbc.updateClusters();
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(nbc.clusters().size());
			for (size_t k = 0; k < nbc.clusters().size(); ++k) {
				for (size_t i : nbc.clusters()[k]) {
					g_multi_pop.pops[k].emplace_back(sol_ptrs[i]);
				}
			}
#endif
			datumUpdated(env);
		}
	}
}
