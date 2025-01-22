#include "gl_cont.h"
#include "../../../../../../core/problem/continuous/continuous.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
	void ContGL::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
		m_input_parameters.add("alpha", new RangedReal(m_alpha, 0, 1, 0.5));
		m_input_parameters.add("beta", new RangedReal(m_beta, 0, 10, 3));
		m_input_parameters.add("gamma", new RangedReal(m_gamma, 0, 20, 6));
		m_input_parameters.add("update scheme", new Enumeration(m_update_scheme,
			{
				"best so far of each Solution",
				"all historical best solutions of each Solution",
				"improved Solutions in the best so far population",
				"all Solutions in the current population"
			}, 
			UpdateScheme::bsf));
	}

	void ContGL::initialize_(Environment *env) {
		Algorithm::initialize_(env);
	}

	void ContGL::run_(Environment *env) {
		m_pop.resize(m_pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop[i]);
		}
#endif
		m_pop.setAlpha(m_alpha);
		m_pop.setBeta(m_beta);
		m_pop.setGamma(m_gamma);
		m_pop.setUpdateScheme(m_update_scheme);
		m_pop.initialize(env, m_random.get());
		m_pop.initializeMemory(env);
		m_pop.evaluate(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		while (!terminating()) {
			m_pop.evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
		}
	}
}
