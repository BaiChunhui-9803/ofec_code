#include "imbalance_ga.h"
#include "../../../../template/selection/single_objective/truncation.h"
#include "../../../../../problem/continuous/free_peaks/free_peaks.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void ImbalanceGA::addInputParameters() {
		m_input_parameters.add("shrink rate", new RangedReal(m_shrink_rate, 0, 1, 1));
		m_input_parameters.add("quarter change", new Bool(m_quarter_change, false));
	}

	void ImbalanceGA::run_(Environment *env) {
		m_pop.initialize(env, m_random.get());
		m_pop.evaluate(env);
#ifdef OFEC_DATUM_OFFSPRING_H
		datum::g_offspring.clear();
#endif // OFEC_DATUM_OFFSPRING_H
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop[i]);
		}
#endif // OFEC_DATUM_MULTI_POP_H
		datumUpdated(env);
		auto &kdtree = CAST_FPs(env->problem())->subspaceTree().tree;
		while (!terminating()) {
			/* determine the number of offsprings in each BoA */
			std::vector<Real> num_offs_each_boa(kdtree->size(), 0);
			for (size_t i = 0; i < m_pop_size; ++i) {
				size_t id_peak = kdtree->getRegionIdx(m_pop[i].variable().vector());
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
						num_offs_each_boa.back() = m_pop_size - sum_shrink;
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
					num_offs_each_boa.back() = m_pop_size - sum_increase;
				}
				else {
					// remain unchanged
				}
			}
			/* generate offsprings */
			for (size_t i = 0; i < m_pop_size; i += 2) {
				/* determine the indexes of peaks to generate offsprings */
				size_t id_peak1_to_generate = 0, id_peak2_to_generate = 0;
				Real rand_pos = m_pop_size * m_random->uniform.next();
				Real accum = 0;
				for (auto prob : num_offs_each_boa) {
					accum += prob;
					if (rand_pos <= accum) {
						break;
					}
					id_peak1_to_generate++;
				}
				rand_pos = m_pop_size * m_random->uniform.next();
				accum = 0;
				for (auto prob : num_offs_each_boa) {
					accum += prob;
					if (rand_pos <= accum) {
						break;
					}
					id_peak2_to_generate++;
				}
				/* generate offsprings in the specified peaks */
				std::vector<size_t> p(2);
				size_t id_peak_off1, id_peak_off2;
				size_t num_attempts = 0;
				while (true) {
					p[0] = m_pop.tournamentSelection(env, m_random.get());
					do { p[1] = m_pop.tournamentSelection(env, m_random.get()); } while (p[1] == p[0]);
					m_pop.crossover(p[0], p[1], m_pop_and_off[i], m_pop_and_off[i + 1], env, m_random.get());
					m_pop.mutate(m_pop_and_off[i], env, m_random.get());
					m_pop.mutate(m_pop_and_off[i + 1], env, m_random.get());
					id_peak_off1 = kdtree->getRegionIdx(m_pop_and_off[i].variable().vector());
					id_peak_off2 = kdtree->getRegionIdx(m_pop_and_off[i + 1].variable().vector());
					if ((id_peak_off1 == id_peak1_to_generate && id_peak_off2 == id_peak2_to_generate) ||
						(id_peak_off1 == id_peak2_to_generate && id_peak_off2 == id_peak1_to_generate)
					) {
						break;
					}
					num_attempts++;
					if (num_attempts > 1000) {
						break;
					}
				}
			}
			for (size_t i = 0; i < m_pop_size; ++i) {
				m_pop_and_off[i + m_pop_size] = m_pop[i];
			}
			for (size_t i = 0; i < m_pop_size; ++i) {
				m_pop_and_off[i].evaluate(env);
			}
			selection::truncation(m_pop, m_pop_and_off, env);
			m_pop.increaseIteration();
#ifdef OFEC_DATUM_OFFSPRING_H
			datum::g_offspring.clear();
			for (size_t i = 0; i < m_pop_size; i++) {
				datum::g_offspring.push_back(&m_pop_and_off[i]);
			}
#endif // OFEC_DATUM_OFFSPRING_H
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < m_pop.size(); ++i) {
				g_multi_pop.pops[0].push_back(&m_pop[i]);
			}
#endif // OFEC_DATUM_MULTI_POP_H
			datumUpdated(env);
		}
	}
}