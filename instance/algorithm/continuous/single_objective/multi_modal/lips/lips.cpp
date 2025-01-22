#include "lips.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void LIPS::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
		m_input_parameters.add("use history nearest", new Bool(m_use_history_nearest, false));
	}

	void LIPS::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		if (m_maximum_evaluations <= 0) {
			m_maximum_evaluations = 1e6;
		}
	}

	void LIPS::run_(Environment *env) {
		m_pop.reset(new SwarmLIP(m_pop_size, env, m_maximum_evaluations));
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop->at(i).pbest());
		}
#endif
		m_pop->setUseHistoryNearest(m_use_history_nearest);
		m_pop->initialize(env, m_random.get());
		m_pop->evaluate(env);
		m_pop->initPbest(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		while (!this->terminating()) {
			m_pop->evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
		}
	}
}