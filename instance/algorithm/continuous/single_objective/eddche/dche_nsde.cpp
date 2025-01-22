#include "dche_nsde.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS


namespace ofec {
	void DCHE_NSDE::addInputParameters() {
		m_input_parameters.add("maximum stagnant iterations", new RangedSizeT(m_max_stagant_iters, 5, 100, 10));
	}

	void DCHE_NSDE::initialize_(Environment *env) {
		DCHE::initialize_(env);
		NSDE::initialize_(env);
	}

	void DCHE_NSDE::run_(Environment *env) {
		m_pop.reset(new PopNSDE(m_pop_size, m_cluster_size, env));
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop->at(i));
		}
#endif
		m_pop->scalingFactor() = 0.9;
		m_pop->crossoverRate() = 0.1;
		Hill *cur_hill = nullptr;
		cur_hill = (*m_random->uniform.nextElem(m_hills.begin(), m_hills.end())).get();
		initPopInHill(cur_hill, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		size_t counter_stagnant = 0;
		while (!terminating()) {
			int rf = m_pop->evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
			if (rf == kNormalEval) {
				bool flag_stagnant = true;
				for (size_t i = 0; i < m_pop->size(); ++i) {
					if (m_pop->at(i).isImproved()) {
						flag_stagnant = false;
						break;
					}
				}
				if (flag_stagnant) {
					counter_stagnant++;
				}
				else {
					counter_stagnant = 0;
				}
				if (counter_stagnant > m_max_stagant_iters) {
					updateSeeds(env);
					updateHills(env);
					cur_hill = (*m_random->uniform.nextElem(m_hills.begin(), m_hills.end())).get();
					initPopInHill(cur_hill, env);
#ifdef OFEC_DATUM_MULTI_POP_H
					datumUpdated(env, g_multi_pop);
#endif
				}
			}
		}
	}

	void DCHE_NSDE::initPopInHill(const Hill *hill, Environment *env) {
		for (size_t i = 0; i < m_pop_size; ++i) {
			randomSolutionInHill(m_pop->at(i), hill);
		}
		m_pop->evaluate(env);
		m_pop->selectSubpop(env);
	}

	void DCHE_NSDE::identifyCandidates(std::list<const Solution<>*> &candidates, Environment *env) {
		if (m_omniscient_seeds) {
			const Real epsilon = 1e-2;
			for (size_t i = 0; i < env->problem()->numberOptimumSolutions(); ++i) {
				auto &optimum = CAST_CONOP(env->problem())->optima()->solution(i);
				size_t id_ssp = m_sp_tree->getRegionIdx(optimum.variable().vector());
				auto hill = m_subspaces.at(id_ssp)->assignedHills().front();
				Solution<> *nearest = nullptr;
				Real min_dis = std::numeric_limits<Real>::max();
				for (size_t i = 0; i < m_pop->size(); ++i) {
					if (isSolutionInHill(m_pop->at(i), hill)) {
						Real dis = optimum.variableDistance(m_pop->at(i), env);
						if (nearest == nullptr || dis < min_dis) {
							min_dis = dis;
							nearest = &m_pop->at(i);
						}
					}
				}
				if (nearest != nullptr && nearest->objectiveDistance(optimum) < epsilon) {
					candidates.push_back(nearest);
				}
			}
		}
		else {
			for (auto &c : m_pop->dist()) {
				int id_best = -1;
				for (auto &p : c) {
					if (id_best == -1 || dominate(m_pop->at(p.second), m_pop->at(id_best), env->problem()->optimizeMode())) {
						id_best = p.second;
					}
				}
				if (id_best != -1) {
					candidates.push_back(&m_pop->at(id_best));
				}
			}
		}
	}
}
