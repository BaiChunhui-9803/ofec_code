#include "dche_hillvallea.h"
#include "../../../../../core/environment/environment.h"
#include "../../../../../core/problem/continuous/continuous.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../datum/datum_inclusion.h"

namespace ofec {
	void DCHE_HillVallEA::addInputParameters() {

	}

	void DCHE_HillVallEA::initialize_(Environment *env) {
		DCHE::initialize_(env);
		HillVallEA::initialize_(env);
	}

	void DCHE_HillVallEA::run_(Environment *env) {
		Hill *cur_hill = nullptr;
		while (!terminating()) {
			cur_hill = (*m_random->uniform.nextElem(m_hills.begin(), m_hills.end())).get();
			auto P = uniformlySampleInHill(cur_hill, env);
#ifdef OFEC_DATUM_NUM_SAMPLES_H
			g_num_samples.value = m_N;
			datumUpdated(env, g_num_samples);
#endif // OFEC_DATUM_NUM_SAMPLES_H
			for (auto &e : m_E) { P.push_back(e); }
			auto S = truncationSelection(P, m_tau, env);
			auto K = clustering(S, env);
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(K.size());
			for (size_t k = 0; k < K.size(); k++) {
				for (size_t i = 0; i < K[k].size(); ++i) {
					g_multi_pop.pops[k].push_back(K[k][i].get());
				}
			}
			datumUpdated(env, g_multi_pop);
#endif
			// Second phase - niche optimization
			for (auto &C : K) {
				if (duplicate(*best(C, env), m_E, env)) {
					continue;
				}
				auto x = coreSearch(C, m_N_c, env, m_random.get());
				if (distinct(*x, m_E, env)) {
					if (betterThanAll(*x, m_E, m_TOL, env)) {
						m_E.clear();
					}
					m_E.push_back(x);
				}
			}
			if (!terminating()) {
				updateSeeds(env);
				updateHills(env);
			}
		}
	}

	std::vector<std::shared_ptr<Solution<>>> DCHE_HillVallEA::uniformlySampleInHill(const Hill *hill, Environment *env) {
		std::vector<std::shared_ptr<Solution<>>> sols;
		for (size_t i = 0; i < m_N && !terminating(); i++) {
			auto new_sol = dynamic_cast<Solution<>*>(env->problem()->createSolution());
			randomSolutionInHill(*new_sol, hill);
			new_sol->evaluate(env);
			sols.emplace_back(new_sol);
		}
		return sols;
	}

	void DCHE_HillVallEA::identifyCandidates(std::list<const Solution<>*> &candidates, Environment *env) {
		if (m_omniscient_seeds) {
			const Real epsilon = 1e-2;
			for (size_t k = 0; k < env->problem()->numberOptimumSolutions(); ++k) {
				auto &optimum = CAST_CONOP(env->problem())->optima()->solution(k);
				size_t id_ssp = m_sp_tree->getRegionIdx(optimum.variable().vector());
				auto hill = m_subspaces.at(id_ssp)->assignedHills().front();
				Solution<> *nearest = nullptr;
				Real min_dis = std::numeric_limits<Real>::max();
				for (auto &e : m_E) {
					if (isSolutionInHill(*e, hill)) {
						Real dis = optimum.variableDistance(*e, env);
						if (nearest == nullptr || dis < min_dis) {
							min_dis = dis;
							nearest = e.get();
						}
					}
				}
				if (nearest != nullptr && nearest->objectiveDistance(optimum) < epsilon) {
					candidates.push_back(nearest);
				}
			}
		}
		else {
			for (auto &e : m_E) {
				candidates.push_back(e.get());
			}
		}
	}
}
