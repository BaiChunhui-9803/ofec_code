#include "dblma.h"
#include "../../../../../../core/environment/environment.h"
#include "../../../../template/selection/diversity/multi.h"
#include "../../../../template/selection/diversity/multi_dynamic.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
	void DBLMA::addInputParameters() {
		m_input_parameters.add("survivor selection scheme", new Enumeration(
			m_sss, { "MULTI", "MULTI_DYNAMIC" }, SurvivorSelectionScheme::kMulti));
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
		m_input_parameters.add("crossover rate", new RangedReal(m_cr, 0, 1, 0.9));
		m_input_parameters.add("mutation rate", new RangedReal(m_mr, 0, 1, 0.5));
		m_input_parameters.add("crossover eta", new RangedReal(m_ceta, 0, 50, 20));
		m_input_parameters.add("mutation eta", new RangedReal(m_meta, 0, 50, 20));
	}

	void DBLMA::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		if (m_pop_size % 2) {
			throw Exception("Population size should be even.");
		}
		if (m_sss == SurvivorSelectionScheme::kMultiDynamic) {
			if (m_maximum_evaluations <= 0) {
				throw Exception("Maximum number of evaluations must be provided.");
			}
		}
		m_step_size.resize(env->problem()->numberVariables());
		for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
			auto &range = CAST_CONOP(env->problem())->range(j);
			m_step_size[j] = (range.second - range.first) / 10;
		}
		m_max_iter_stagnant = 10;
	}

	void DBLMA::run_(Environment *env) {
		m_pop.resize(m_pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop[i]);
		}
#endif
		m_pop_and_off.resize(2 * m_pop_size, env);
		m_pop.setRate(m_cr, m_mr);
		m_pop.setEta(m_ceta, m_meta);
		m_pop.initialize(env, m_random.get());
		m_pop.evaluate(env);
		m_neighbor.reset(dynamic_cast<Solution<>*>(env->problem()->createSolution()));
		localSearch(m_pop, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		while (!terminating()) {
			m_pop.reproduction(m_pop_and_off, env, m_random.get());
			for (size_t i = 0; i < m_pop_size; ++i) {
				m_pop_and_off[i].evaluate(env);
			}
			localSearch(m_pop_and_off, env);
			if (m_sss == SurvivorSelectionScheme::kMulti) {
				selection::diversity::multi(m_pop, m_pop_and_off, env, m_random.get());
			}
			else {
				Real elapsed = (Real)m_evaluations / m_maximum_evaluations;
				selection::diversity::multiDynamic(m_pop, m_pop_and_off, env, m_random.get(), elapsed);
			}
			m_pop.increaseIteration();
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
		}
	}

	void DBLMA::localSearch(Population<Solution<>> &pop, Environment *env) {
		for (size_t i = 0; i < m_pop_size; ++i) {
			size_t iter_stagnant = 0;
			while (iter_stagnant < m_max_iter_stagnant && !terminating()) {
				do {
					for (size_t j = 0; j < m_neighbor->variable().size(); j++)
						m_neighbor->variable()[j] = m_random->normal.nextNonStd(pop[i].variable()[j], m_step_size[j]);
				} while (CAST_CONOP(env->problem())->boundaryViolated(*m_neighbor));
				m_neighbor->evaluate(env);
				if (dominate(*m_neighbor, pop[i], env->problem()->optimizeMode())) {
					pop[i] = *m_neighbor;
					iter_stagnant = 0;
				}
				else {
					iter_stagnant++;
				}
			}
		}
	}
}