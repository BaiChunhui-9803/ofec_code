#include "nbde.h"
#include "../../../../../../utility/clustering/nbc.h"
#include "../../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../template/selection/multi_objective/nsgaii.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void NBDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void NBDE::run_(Environment *env) {
		m_pop.resize(m_pop_size, env);
		m_pop.initialize(env, m_random.get());
		m_pop.evaluate(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop.at(i));
		}
#endif
		datumUpdated(env);
		std::vector<IndividualDE> pop_cmbnd(2 * m_pop_size, IndividualDE(env));
		std::vector<OptimizeMode> opt_mode = { env->problem()->optimizeMode(0) == OptimizeMode::kMaximize ? OptimizeMode::kMinimize : OptimizeMode::kMaximize, OptimizeMode::kMaximize};
		NSGAII nsgaii(2, opt_mode);	
		std::vector<std::multimap<Real, size_t>> nearest(2 * m_pop_size);
		std::vector<std::vector<Real>> objs(2 * m_pop_size, std::vector<Real>(2));
		while (!terminating()) {
			m_pop.evolve(env, m_random.get());
			for (size_t i = 0; i < m_pop_size; ++i) {
				pop_cmbnd[i] = m_pop[i];
				pop_cmbnd[m_pop_size + i].solution() = m_pop[i].trial();
			}
			for (size_t i = 0; i < pop_cmbnd.size(); ++i) {
				objs[i][0] = pop_cmbnd[i].objective(0);
				nearest[i].clear();
				Real dis;
				for (size_t j = 0; j < pop_cmbnd.size(); j++) {
					if (j != i) {
						dis = pop_cmbnd[i].variableDistance(pop_cmbnd[j], env);
						nearest[i].emplace(dis, j);
					}
				}
				//objs[i][0] = 0;
				//auto iter = nearest[i].begin();
				//for (size_t j = 0; j < 2 * m_k; j++, iter++) {
				//	if (pop_cmbnd[iter->second].dominate(pop_cmbnd[i], env))
				//		objs[i][0]++;
				//}
				objs[i][1] = nearest[i].begin()->first;
				//objs[i][1] = 1;
			}
			nsgaii.survivorSelection(m_pop, pop_cmbnd, objs);
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < m_pop.size(); ++i) {
				g_multi_pop.pops[0].push_back(&m_pop.at(i));
			}
#endif
			datumUpdated(env);
		}
	}
}
