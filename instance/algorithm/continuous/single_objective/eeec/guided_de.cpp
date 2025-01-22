#include "guided_de.h"
#include "../../../../../problem/continuous/free_peaks/free_peaks.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS


namespace ofec {
	void GuidedDE::addInputParameters() {}

	void GuidedDE::run_(Environment *env) {
		auto &kdtree = CAST_FPs(env->problem())->subspaceTree().tree;
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
		std::vector<Solution<>> offspring(m_population_size);
		while (!terminating()) {
			/* determine the ratio of number of offsprings in each BoA */
			std::vector<Real> objs(m_population_size);
			for (size_t i = 0; i < m_population_size; ++i) {
				objs[i] = m_population->at(i).objective(0);
			}
			std::nth_element(objs.begin(), objs.begin() + objs.size() / 2, objs.end());
			Real threshold = *(objs.begin() + objs.size() / 2);
			std::vector<size_t> num_over_thrd_each_boa(kdtree->size(), 0);
			for (size_t i = 0; i < m_population_size; ++i) {
				if (m_population->at(i).objective(0) >= threshold) {
					size_t id_peak = kdtree->getRegionIdx(m_population->at(i).variable().vector());
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
			/*  */
			for (size_t i = 0; i < m_population_size; ++i) {
				size_t id_boa_to_generate = 0;
				Real rand_pos = sum_ratio * m_random->uniform.next();
				Real accum = 0;
				for (auto prob : ratio_each_boa) {
					accum += prob;
					if (rand_pos <= accum) {
						break;
					}
					id_boa_to_generate++;
				}
				size_t rand_ind;
				size_t num_attempts = 0;
				while (true) {
					rand_ind = m_random->uniform.next() * m_population_size;
					m_population->mutate(rand_ind, m_random.get(), env);
					m_population->recombine(rand_ind, m_random.get(), env);
					size_t id_boa_off = kdtree->getRegionIdx(
						m_population->at(rand_ind).trial().variable().vector());
					if (id_boa_off == id_boa_to_generate) {
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
