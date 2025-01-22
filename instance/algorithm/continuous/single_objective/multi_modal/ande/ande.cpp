#include "ande.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS
#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void ANDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void ANDE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
	}

	void ANDE::run_(Environment *env) {
		m_pop.reset(new PopANDE(m_pop_size, env, 0.9, 100, 30));
		m_pop->initialize(env, m_random.get());
		m_pop->evaluate(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop->at(i));
		}
		datumUpdated(env, g_multi_pop);
#endif
		auto &clusters = m_pop->clusters();
		while (!terminating()) {
			m_pop->clustering(env);
#ifdef OFEC_DATUM_SEEDS_H
			g_seeds.sols.clear();
			for (size_t k = 0; k < clusters.size(); ++k) {
				SolutionBase *new_seed = nullptr;
				for (size_t i = 0; i < clusters[k].size(); ++i) {
					if (!new_seed || dominate(m_pop->at(clusters[k][i]), *new_seed, env->problem()->optimizeMode())) {
						new_seed = &m_pop->at(clusters[k][i]);
					}
				}
				if (new_seed) {
					g_seeds.sols.push_back(new_seed);
				}
			}
			datumUpdated(env, g_seeds);
#endif // OFEC_DATUM_SEEDS_H
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(clusters.size());
			for (size_t k = 0; k < clusters.size(); ++k) {
				for (size_t i = 0; i < clusters[k].size(); ++i) {
					g_multi_pop.pops[k].push_back(&m_pop->at(clusters[k][i]));
				}
			}
			datumUpdated(env, g_multi_pop);
#endif
			m_pop->evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(clusters.size());
			for (size_t k = 0; k < clusters.size(); ++k) {
				for (size_t i = 0; i < clusters[k].size(); ++i) {
					g_multi_pop.pops[k].push_back(&m_pop->at(clusters[k][i]));
				}
			}
			datumUpdated(env, g_multi_pop);
#endif
		}
	}
}