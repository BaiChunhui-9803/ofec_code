#include "crowding_de.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS
#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
	PopCrowdingDE::PopCrowdingDE(size_t size_pop, Environment *env) : 
		PopulationDE(size_pop, env) {}

	int PopCrowdingDE::evolve(Environment *env, Random *rnd) {
		if (m_individuals.size() < 5) {
			throw Exception("the population size cannot be smaller than 5@PopCrowdingDE::evolve()");
		}
		int tag = kNormalEval;
		for (size_t i = 0; i < size(); ++i) {
			mutate(i, rnd, env);
			recombine(i, rnd, env);
			tag = m_individuals[i]->trial().evaluate(env);
			size_t nn = nearestNeighbour(i, env);
			if (dominate(m_individuals[i]->trial(), *m_individuals[nn], env->problem()->optimizeMode())) {
				m_individuals[nn]->solution() = m_individuals[i]->trial();
			}
			if (tag & kTerminate) {
				return tag;
			}
		}
		m_iteration++;
		return tag;
	}

	size_t PopCrowdingDE::nearestNeighbour(int i, Environment *env) {
		size_t nn = 0;
		Real min_dis = m_individuals[i]->trial().variableDistance(*m_individuals[nn], env);
		for (size_t j = 1; j < m_individuals.size(); ++j) {
			Real dis = m_individuals[i]->trial().variableDistance(*m_individuals[j], env);
			if (dis < min_dis) {
				nn = j;
				min_dis = dis;
			}
		}
		return nn;
	}

	void CrowdingDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void CrowdingDE::run_(Environment *env) {
		PopCrowdingDE pop(m_pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop_size; ++i) {
			g_multi_pop.pops[0].push_back(&pop[i]);
		}
#endif
		pop.setParameter(0.9, 0.5);
		pop.initialize(env, m_random.get());
		pop.evaluate(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		while (!terminating()) {
			pop.evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		}
	}
}

