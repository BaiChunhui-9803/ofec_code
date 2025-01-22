#include "imbalance_de.h"
#include "../../../../../problem/continuous/free_peaks/free_peaks.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void ImbalanceDE::addInputParameters() {
		m_input_parameters.add("shrink rate", new RangedReal(m_shrink_rate, 0, 1, 1));
		m_input_parameters.add("quarter change", new Bool(m_quarter_change, false));
	}

	void ImbalanceDE::run_(Environment *env) {
		m_population->initialize(env, m_random.get());
		m_population->evaluate(env);
#ifdef OFEC_DATUM_OFFSPRING_H
		datum::g_offspring.clear();
#endif // OFEC_DATUM_OFFSPRING_H
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_population->size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_population->at(i));
		}
#endif // OFEC_DATUM_MULTI_POP_H
		datumUpdated(env);
		auto &kdtree = CAST_FPs(env->problem())->subspaceTree().tree;
		std::vector<Solution<>> offspring(m_population_size);
		bool terminate = false;
		while (!terminating()) {
			/* determine the number of offsprings in each BoA */
			std::vector<Real> num_offs_each_boa(kdtree->size(), 0);
			for (size_t i = 0; i < m_population_size; ++i) {
				size_t id_peak = kdtree->getRegionIdx(m_population->at(i).variable().vector());
				num_offs_each_boa[id_peak]++;
			}
			if (m_quarter_change) {
				if (m_evaluations <= m_maximum_evaluations * 0.25) {
					if (num_offs_each_boa.back() > 0) {
						Real sum_shrink = 0;
						for (size_t i = 0; i < kdtree->size() - 1; ++i) {
							num_offs_each_boa[i] *= m_shrink_rate;
							sum_shrink += num_offs_each_boa[i];
						}
						num_offs_each_boa.back() = m_population_size - sum_shrink;
					}
				}
				else if (m_evaluations <= m_maximum_evaluations * 0.5) {
					// remain unchanged
				}
				else if (m_evaluations <= m_maximum_evaluations * 0.75) {
					Real sum_increase = 0;
					for (size_t i = 0; i < kdtree->size() - 1; ++i) {
						num_offs_each_boa[i] /= m_shrink_rate;
						sum_increase += num_offs_each_boa[i];
					}
					num_offs_each_boa.back() = m_population_size - sum_increase;
				}
				else {
					// remain unchanged
				}
			}
			for (size_t i = 0; i < m_population_size; ++i) {
				size_t id_peak_to_generate = 0;
				Real rand_pos = m_population_size * m_random->uniform.next();
				Real accum = 0;
				for (auto prob : num_offs_each_boa) {
					accum += prob;
					if (rand_pos <= accum) {
						break;
					}
					id_peak_to_generate++;
				}
				size_t rand_ind;
				size_t num_attempts = 0;
				while (true) {
					rand_ind = m_random->uniform.next() * m_population_size;
					m_population->mutate(rand_ind, m_random.get(), env);
					m_population->recombine(rand_ind, m_random.get(), env);
					size_t id_peak_off = kdtree->getRegionIdx(
						m_population->at(rand_ind).trial().variable().vector());
					if (id_peak_off == id_peak_to_generate) {
						break;
					}
					num_attempts++;
					if (num_attempts > 1000) {
						break;
					}
				}
				offspring[i] = m_population->at(rand_ind).trial();
			}
			for (size_t i = 0; i < m_population_size; ++i) {
				m_population->at(i).trial() = offspring[i];
				int tag = m_population->at(i).select(env);
				if (!(tag & kNormalEval)) {
					break;
				}
			}
#ifdef OFEC_DATUM_OFFSPRING_H
			datum::g_offspring.clear();
			for (size_t i = 0; i < m_population->size(); i++) {
				datum::g_offspring.push_back(&m_population->at(i).trial());
			}
#endif // OFEC_DATUM_OFFSPRING_H
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < m_population->size(); ++i) {
				g_multi_pop.pops[0].push_back(&m_population->at(i));
			}
#endif // OFEC_DATUM_MULTI_POP_H
			datumUpdated(env);
		}
	}
}
