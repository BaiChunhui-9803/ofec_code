#include "de_nrand_1.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
	void DE_nrand_1::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
		m_input_parameters.add("scaling factor", new RangedReal(m_scaling_factor, 0, 2, 0.5));
		m_input_parameters.add("crossover rate", new RangedReal(m_crossover_rate, 0, 1, 0.6));
	}

	void DE_nrand_1::initialize_(Environment *env) {
		Algorithm::initialize_(env);
	}

	void DE_nrand_1::run_(Environment *env) {
		initPop(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		while (!terminating()) {
			m_pop->evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		}
	}

	void DE_nrand_1::initPop(Environment *env) {
		m_pop.reset(new PopDE_nrand_1(m_pop_size, env));
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop->at(i));
		}
#endif
		m_pop->initialize(env, m_random.get());
		m_pop->evaluate(env);
		m_pop->crossoverRate() = m_crossover_rate;
		m_pop->scalingFactor() = m_scaling_factor;
	}
}
