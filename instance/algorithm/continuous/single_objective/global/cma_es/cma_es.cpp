#include "cma_es.h"
#include "cma_es_pop.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void CMA_ES::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 2, 1000, 5));
	}

	void CMA_ES::initialize_(Environment *env) {
		Algorithm::initialize_(env);
	}

	void CMA_ES::run_(Environment *env) {
		PopCMA_ES pop(m_pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop_size; ++i) {
			g_multi_pop.pops[0].push_back(&pop[i]);
	}
#endif // OFEC_DATUM_MULTI_POP_H
#ifdef OFEC_DATUM_OFFSPRING_H
		g_offspring.sols.clear();
		for (size_t i = 0; i < m_pop_size; ++i) {
			g_offspring.sols.push_back(&pop[i]);
		}
#endif // OFEC_DATUM_OFFSPRING_H
		pop.initialize(env, m_random.get());
		pop.reproduce(env, m_random.get());
		pop.evaluate(env);
		while (!terminating()) {
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
			pop.evolve(env, m_random.get());
#ifdef OFEC_DATUM_OFFSPRING_H
			datumUpdated(env, g_offspring);
#endif // OFEC_DATUM_OFFSPRING_H
			if (pop.conditionNumber() > 1e14 || pop.equalFitness()) {
				terminate();
			}
		}
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
	}
}