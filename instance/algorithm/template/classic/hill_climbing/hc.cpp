#include "hc.h"
#include "../../../../../utility/functional.h"
#include "../../../../../core/environment/environment.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void BaseHC::initialize_(Environment *env) {
		Algorithm::initialize_(env);
	}

	void BaseHC::run_(Environment *env) {
		m_cur.reset(env->problem()->createSolution());
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		g_multi_pop.pops[0].push_back(m_cur.get());
#endif
		m_cur->initialize(env, m_random.get());
		m_cur->evaluate(env, this);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		m_neighbour.reset(env->problem()->createSolution());
		while (!terminating()) {
			pickNeighbour(env);
			m_neighbour->evaluate(env, this);
			if (dominate(*m_neighbour, *m_cur, env->problem()->optimizeMode())) {
				replaceCurrent();
			}
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
		}
	}
}