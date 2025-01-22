#include "knpso_rd.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void kNPSO_RD::addInputParameters() {
		m_input_parameters.add("maximum stagnant iterations", new RangedSizeT(m_max_stagant_iters, 5, 100, 10));
	}

	void kNPSO_RD::initialize_(Environment *env) {
		kNPSO::initialize_(env);

	}

	void kNPSO_RD::run_(Environment *env) {
		size_t swarm_size = m_swarm_size;
		m_swarm.setNumNeighbors(m_num_nbrs);
		m_swarm.resize(swarm_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_swarm.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_swarm.at(i).pbest());
		}
#endif
		m_swarm.initialize(env, m_random.get());
		m_swarm.initVelocity(env, m_random.get());
		m_swarm.initVelocityMax(env, m_random.get());
		m_swarm.evaluate(env);
		m_swarm.initPbest(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		size_t counter_stagnant = 0;
		while (!terminating()) {
			int rf = m_swarm.evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
			if (rf == kNormalEval) {
				bool flag_stagnant = true;
				for (size_t i = 0; i < m_swarm.size(); ++i) {
					if (m_swarm[i].isImproved()) {
						flag_stagnant = false;
						break;
					}
				}
				if (flag_stagnant) {
					counter_stagnant++;
				}
				else {
					counter_stagnant = 0;
				}
				if (counter_stagnant > m_max_stagant_iters) {
					swarm_size *= 2;
					m_swarm.resize(swarm_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
					g_multi_pop.pops.clear();
					g_multi_pop.pops.resize(1);
					for (size_t i = 0; i < m_swarm.size(); ++i) {
						g_multi_pop.pops[0].push_back(&m_swarm.at(i).pbest());
					}
#endif
					m_swarm.initialize(env, m_random.get());
					m_swarm.initVelocity(env, m_random.get());
					m_swarm.initVelocityMax(env, m_random.get());
					m_swarm.evaluate(env);
					m_swarm.initPbest(env);
#ifdef OFEC_DATUM_MULTI_POP_H
					datumUpdated(env, g_multi_pop);
#endif
				}
			}
		}
	}
}
