#include "guided_ga.h"
#include "../../../../../problem/continuous/free_peaks/free_peaks.h"
#include "../../../../template/selection/single_objective/truncation.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS


namespace ofec {
	void GuidedGA::addInputParameters() {}

	void GuidedGA::run_(Environment *env) {
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
			/* determine the ratio of number of offsprings in each BoA */
			std::vector<Real> objs(m_pop_size);
			for (size_t i = 0; i < m_pop_size; ++i) {
				objs[i] = m_pop.at(i).objective(0);
			}
			std::nth_element(objs.begin(), objs.begin() + objs.size() / 2, objs.end());
			Real threshold = *(objs.begin() + objs.size() / 2);
			std::vector<size_t> num_over_thrd_each_boa(kdtree->size(), 0);
			for (size_t i = 0; i < m_pop_size; ++i) {
				if (m_pop.at(i).objective(0) >= threshold) {
					size_t id_peak = kdtree->getRegionIdx(m_pop.at(i).variable().vector());
					num_over_thrd_each_boa[id_peak]++;
				}
			}
			std::list<size_t> empty_boas;
			for (size_t i = 0; i < kdtree->size(); ++i) {
				if (num_over_thrd_each_boa[i] == 0) {
					empty_boas.push_back(i);
				}
			}
			std::vector<Real> ratio_each_boa(kdtree->size());
			if (empty_boas.empty()) {
				for (size_t i = 0; i < kdtree->size(); ++i) {
					ratio_each_boa[i] = 1.0 / num_over_thrd_each_boa[i];
				}
			}
			else {
				ratio_each_boa.assign(kdtree->size(), 0);
				for (size_t i : empty_boas) {
					ratio_each_boa[i] = 1.0;
				}
			}
			Real sum_ratio = 0;
			for (size_t i = 0; i < kdtree->size(); ++i) {
				sum_ratio += ratio_each_boa[i];
			}
			for (size_t i = 0; i < m_pop_size; i += 2) {
				/* determine the indexes of peaks to generate offsprings */
				size_t id_peak1_to_generate = 0, id_peak2_to_generate = 0;
				Real rand_pos = sum_ratio * m_random->uniform.next();
				Real accum = 0;
				for (auto prob : ratio_each_boa) {
					accum += prob;
					if (rand_pos <= accum) {
						break;
					}
					id_peak1_to_generate++;
				}
				rand_pos = sum_ratio * m_random->uniform.next();
				accum = 0;
				for (auto prob : ratio_each_boa) {
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
