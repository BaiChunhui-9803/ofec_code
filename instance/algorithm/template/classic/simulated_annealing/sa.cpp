#include "sa.h"
#include "../../../../../utility/functional.h"
#include "../../../../../core/environment/environment.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void BaseSA::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		if (m_maximum_evaluations <= 0) {
			throw Exception("Maximum number of evaluations must be provided.");
		}
	}

	void BaseSA::run_(Environment *env) {
		m_cur.reset(env->problem()->createSolution());
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		g_multi_pop.pops[0].push_back(m_cur.get());
#endif
		m_cur->initialize(env, m_random.get());
		m_cur->evaluate(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		m_neighbour.reset(env->problem()->createSolution());
		Real t, p;
		while (!terminating()) {
			pickNeighbour(env);
			m_neighbour->evaluate(env);
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
			t = temperature((Real)m_evaluations / m_maximum_evaluations);
			p = criterion(t, env);
			if (p >= m_random->uniform.next())
				replaceCurrent();
		}
	}

	void BaseSA::record() {}

	Real BaseSA::temperature(Real r) {
		return -log10(r);
	}

	Real BaseSA::criterion(Real t, Environment *env) {
		if (dominate(*m_neighbour, *m_cur, env->problem()->optimizeMode())) {
			return 1.0;
		}
		else {
			return exp(-m_neighbour->objectiveDistance(*m_cur) / t);
		}
	}
}