#include "jde.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
	void jDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void jDE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_pop.reset(new PopJDE(m_pop_size, env));
	}

	void jDE::run_(Environment *env) {
#ifdef OFEC_DATUM_JDE_POP_H
		g_jde_pop.value = m_pop.get();
#endif // OFEC_DATUM_JDE_POP_H
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop->at(i));
		}
#endif // OFEC_DATUM_MULTI_POP_H
		m_pop->initialize(env, m_random.get());
		m_pop->evaluate(env);
#ifdef OFEC_DATUM_JDE_POP_H
		datumUpdated(env, g_jde_pop);
#endif // OFEC_DATUM_JDE_POP_H
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		while (!terminating()) {
			m_pop->evolve(env, m_random.get());
#ifdef OFEC_DATUM_JDE_POP_H
			datumUpdated(env, g_jde_pop);
#endif // OFEC_DATUM_JDE_POP_H
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		}
	}
}


