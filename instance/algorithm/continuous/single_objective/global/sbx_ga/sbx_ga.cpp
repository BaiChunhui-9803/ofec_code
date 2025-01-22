#include "sbx_ga.h"
#include "../../../../template/selection/single_objective/truncation.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
	void SBX_GA::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
		m_input_parameters.add("crossover rate", new RangedReal(m_cr, 0, 1, 0.9));
		m_input_parameters.add("mutation rate", new RangedReal(m_mr, 0, 1, 0.5));
		m_input_parameters.add("crossover eta", new RangedReal(m_ceta, 0, 50, 20));
		m_input_parameters.add("mutation eta", new RangedReal(m_meta, 0, 50, 20));
	}

	void SBX_GA::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		if (m_pop_size % 2) {
			throw Exception("Population size of NSGAII should be even.");
		}
		m_pop.resize(m_pop_size, env);
		m_pop_and_off.resize(2 * m_pop_size, env);
		m_pop.setRate(m_cr, m_mr);
		m_pop.setEta(m_ceta, m_meta);
	}

	void SBX_GA::run_(Environment *env) {
#ifdef OFEC_DATUM_OFFSPRING_H
		g_offspring.sols.clear();
#endif // OFEC_DATUM_OFFSPRING_H

		m_pop.initialize(env, m_random.get());
		m_pop.evaluate(env);
#ifdef OFEC_DATUM_OFFSPRING_H
		datumUpdated(env, g_offspring);
#endif // OFEC_DATUM_OFFSPRING_H
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop[i]);
		}
		datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		while (!terminating()) {
			m_pop.reproduction(m_pop_and_off, env, m_random.get());
			for (size_t i = 0; i < m_pop_size; ++i) {
				m_pop_and_off[i].evaluate(env);
			}
			selection::truncation(m_pop, m_pop_and_off, env);
			m_pop.increaseIteration();
#ifdef OFEC_DATUM_OFFSPRING_H
			datumUpdated(env, g_offspring);
#endif // OFEC_DATUM_OFFSPRING_H
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < m_pop.size(); ++i) {
				g_multi_pop.pops[0].push_back(&m_pop[i]);
			}
			datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		}
	}
}