#include "lrde.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void LRDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
		m_input_parameters.add("radius", new RangedReal(m_repulsion_radius, 0, 1, 0.35));
		m_input_parameters.add("repulsion probability", new RangedReal(m_repulsion_probabilty, 0, 1.0, 0.8));
		m_input_parameters.add("maximum stagnant iterations", new RangedSizeT(m_max_stagant_iters, 5, 100, 10));
	}

	void LRDE::run_(Environment *env) {
		m_pop.resize(m_pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop.at(i));
		}
#endif
		m_archive.clear();
		m_best.reset();
		initializePopulation(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		int rf = kNormalEval;
		bool best_improved ;
		size_t counter_stagnant = 0;
		while (!terminating()) {
			best_improved = false;
			for (size_t i = 0; i < m_pop.size(); ++i) {
				do {
					m_pop.mutate(i, m_random.get(), env);
					m_pop.recombine(i, m_random.get(), env);
				} while (reject(m_pop[i].trial(), env));
				rf = m_pop[i].select(env);
				if (rf != kNormalEval) {
					break;
				}
				if (m_pop[i].isImproved() && (!m_best || dominate(m_pop[i], *m_best, env->problem()->optimizeMode()))) {
					m_best.reset(new Solution<>(m_pop[i]));
					best_improved = true;
				}
			}
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
			if (rf == kNormalEval) {
				if (!best_improved) {
					counter_stagnant++;
				}
				else {
					counter_stagnant = 0;
				}
				if (counter_stagnant > m_max_stagant_iters) {
					updateArchive(env);
					initializePopulation(env);
#ifdef OFEC_DATUM_MULTI_POP_H
					datumUpdated(env, g_multi_pop);
#endif
				}
			}
		}
	}

	bool LRDE::reject(const Solution<> &sol, Environment *env) const {
		for (auto &s : m_archive) {
			if (s.variableDistance(sol, env) < 
				m_repulsion_radius * env->problem()->maximumVariableDistance()) 
			{
				if (m_random->uniform.next() < m_repulsion_probabilty)
					return true;
			}
		}
		return false;
	}

	void LRDE::updateArchive(Environment *env) {
		size_t id_best = 0, id_worst = 0;
		for (size_t i = 1; i < m_pop.size(); ++i) {
			if (dominate(m_pop[i], m_pop[id_best], env->problem()->optimizeMode())) {
				id_best = i;
			}
			if (dominate(m_pop[id_worst], m_pop[i], env->problem()->optimizeMode())) {
				id_worst = i;
			}
		}
		if (m_pop[id_best].objectiveDistance(m_pop[id_worst]) < 1e-3) {
			bool duplicate = false;
			for (auto &s : m_archive) {
				if (s.variableDistance(m_pop[id_best], env) < 1e-4) {
					if (dominate(m_pop[id_best], s, env->problem()->optimizeMode())) {
						s = m_pop[id_best];
					}
					duplicate = true;
					break;
				}
			}
			if (!duplicate) {
				m_archive.emplace_back(m_pop[id_best]);
			}
		}
	}

	void LRDE::initializePopulation(Environment *env) {
		for (size_t i = 0; i < m_pop.size(); ++i) {
			do {
				m_pop[i].initialize(env, m_random.get());
			} while (reject(m_pop[i], env));
		}
		m_pop.evaluate(env);
		m_best.reset();
	}
}
