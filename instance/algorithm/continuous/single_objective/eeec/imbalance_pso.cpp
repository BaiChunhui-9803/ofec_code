#include "imbalance_pso.h"
#include "../../../../../problem/continuous/free_peaks/free_peaks.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void ImbalancePSO::addInputParameters() {
		m_input_parameters.add("shrink rate", new RangedReal(m_shrink_rate, 0, 1, 1));
		m_input_parameters.add("quarter change", new Bool(m_quarter_change, false));
	}

	void ImbalancePSO::run_(Environment *env) {
		m_pop->initialize(env, m_random.get());
		m_pop->initVelocity(env, m_random.get());
		m_pop->evaluate(env);
		m_pop->initPbest(env);
#ifdef OFEC_DATUM_OFFSPRING_H
		datum::g_offspring.clear();
#endif // OFEC_DATUM_OFFSPRING_H
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop->at(i));
		}
#endif
		datumUpdated(env);
		auto &kdtree = CAST_FPs(env->problem())->subspaceTree().tree;
		while (!terminating()) {
			/* determine the number of offsprings in each BoA */
			std::vector<Real> num_offs_each_boa(kdtree->size(), 0);
			for (size_t i = 0; i < m_pop_size; ++i) {
				size_t id_peak = kdtree->getRegionIdx(m_pop->at(i).pbest().variable().vector());
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
			//generate a permutation of particle index
			//std::vector<int> rindex(this->m_individuals.size());
			//std::iota(rindex.begin(), rindex.end(), 0);
			//rnd->uniform.shuffle(rindex.begin(), rindex.end());

			//this->setNeighborhood(rnd);

			//bool flag = false;
			//for (int i = 0; i < this->m_individuals.size(); i++) {
			//	auto &x = neighborhoodBest(rindex[i], env);

			//	this->m_individuals[rindex[i]]->nextVelocity(&x, m_weight, m_accelerator1, m_accelerator2, rnd);
			//	this->m_individuals[rindex[i]]->move();
			//	this->m_individuals[rindex[i]]->clampVelocity(env, rnd);

			//	rf = this->m_individuals[rindex[i]]->evaluate(env);

			//	if (dominate(*this->m_individuals[rindex[i]], this->m_individuals[rindex[i]]->pbest(), env->problem()->optimizeMode())) {
			//		this->m_individuals[rindex[i]]->pbest() = *(this->m_individuals[rindex[i]]);
			//		if (dominate(this->m_individuals[rindex[i]]->pbest(), this->m_best_particle->pbest(), env->problem()->optimizeMode())) {
			//			flag = true;
			//		}
			//	}

			//	if (rf != kNormalEval) break;
			//}
			//m_flag_best_impr = flag;
			//this->m_iteration++;
			//return rf;
#ifdef OFEC_DATUM_OFFSPRING_H
			datum::g_offspring.clear();
			for (size_t i = 0; i < m_pop->size(); i++) {
				datum::g_offspring.push_back(&m_pop->at(i));
			}
#endif // OFEC_DATUM_OFFSPRING_H
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < m_pop->size(); ++i) {
				g_multi_pop.pops[0].push_back(&m_pop->at(i).pbest());
			}
#endif
			datumUpdated(env);
		}
	}
}
