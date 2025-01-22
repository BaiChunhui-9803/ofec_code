#include "dche_nea2.h"
#include "../../../../../utility/clustering/nbc.h"
#include "../../../../../core/problem/continuous/continuous.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../datum/datum_inclusion.h"

namespace ofec {
	void DCHE_NEA2::addInputParameters() {

	}

	void DCHE_NEA2::initialize_(Environment *env) {
		DCHE::initialize_(env);
		NEA2::initialize_(env);
	}

	void DCHE_NEA2::run_(Environment *env) {
		Hill *cur_hill = nullptr;
		while (!terminating()) {
			cur_hill = (*m_random->uniform.nextElem(m_hills.begin(), m_hills.end())).get();
			addSubpopsInHill(cur_hill, env);
			while (!terminating() && m_subpops.isActive()) {
				int rf = kNormalEval;
				for (size_t k = 0; k < m_subpops.size(); ++k) {
					if (!m_subpops[k].isActive()) continue;
					rf = m_subpops[k].evolve(env, m_random.get());
					if (rf & kTerminate) {
						break;
					}
					for (size_t k = 0; k < m_subpops.size(); ++k) {
						if (m_subpops[k].isActive() && stopTolFun(m_subpops[k])) {
							m_subpops[k].setActive(false);
						}
					}
#ifdef OFEC_DATUM_MULTI_POP_H
					datumUpdated(env, g_multi_pop);
#endif
				}
				if (rf & kTerminate) {
					break;
				}
			}
			if (!terminating()) {
				updateSeeds(env);
				updateHills(env);
			}
		}
	}

	void DCHE_NEA2::addSubpopsInHill(const Hill *hill, Environment *env) {
		size_t num_vars = env->problem()->numberVariables();
		Population<Solution<>> global_sample(m_pop_size, env, num_vars);
		for (size_t i = 0; i < m_pop_size; ++i) {
			randomSolutionInHill(global_sample[i], hill);
			global_sample[i].evaluate(env);
		}
		NBC nbc(2.0, NBC::kByKDTree, NBC::kByMean, true);
		nbc.setData(global_sample);
		nbc.clustering(env);
		auto &clusters = nbc.clusters();
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(clusters.size());
		for (size_t k = 0; k < clusters.size(); k++) {
			for (size_t i = 0; i < clusters[k].size(); ++i) {
				g_multi_pop.pops[k].push_back(&global_sample[clusters[k][i]]);
			}
		}
		datumUpdated(env, g_multi_pop);
#endif
		std::vector<size_t> filtered_start_Solutions;
		for (auto &cluster : clusters) {
			auto iter = cluster.begin();
			size_t center = *iter;
			while (++iter != cluster.end()) {
				if (dominate(global_sample[*iter], global_sample[center], env->problem()->optimizeMode())) {
					center = *iter;
				}
			}
			filtered_start_Solutions.push_back(center);
		}
		std::vector<Real> step_size(num_vars);
		std::vector<std::pair<Real, Real>> bnd(num_vars, { std::numeric_limits<Real>::max(),-std::numeric_limits<Real>::max() });
		Solution<> tmp_sol(1, 0, num_vars);
		for (size_t i = 0; i < 1000; ++i) {
			randomSolutionInHill(tmp_sol, hill);
			for (size_t j = 0; j < num_vars; ++j) {
				if (tmp_sol.variable()[j] < bnd[j].first) {
					bnd[j].first = tmp_sol.variable()[j];
				}
				if (tmp_sol.variable()[j] > bnd[j].second) {
					bnd[j].second = tmp_sol.variable()[j];
				}
			}
		}
		for (size_t j = 0; j < num_vars; ++j) {
			step_size[j] = 0.05 * (bnd[j].second - bnd[j].first);
		}
		size_t size_subpop = 4 + int(3 * log(num_vars));
		m_subpops.clear();
		for (size_t start_ind : filtered_start_Solutions) {
			auto pop = std::make_unique<PopCMA_ES>(size_subpop, env);
			pop->initializeByNonStd(env, m_random.get(), global_sample[start_ind].variable().vector(), step_size);
			pop->reproduce(env, m_random.get());
			pop->evaluate(env);
			m_subpops.append(pop);
		}
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(m_subpops.size());
		for (size_t k = 0; k < m_subpops.size(); k++) {
			for (size_t i = 0; i < m_subpops[k].size(); ++i) {
				g_multi_pop.pops[k].push_back(&m_subpops[k][i]);
			}
		}
		datumUpdated(env, g_multi_pop);
#endif
	}

	void DCHE_NEA2::identifyCandidates(std::list<const Solution<>*> &candidates, Environment *env) {
		if (m_omniscient_seeds) {
			const Real epsilon = 1e-2;
			for (size_t k = 0; k < env->problem()->numberOptimumSolutions(); ++k) {
				auto &optimum = CAST_CONOP(env->problem())->optima()->solution(k);
				size_t id_ssp = m_sp_tree->getRegionIdx(optimum.variable().vector());
				auto hill = m_subspaces.at(id_ssp)->assignedHills().front();
				Solution<> *nearest = nullptr;
				Real min_dis = std::numeric_limits<Real>::max();
				for (auto &subpop : m_subpops) {
					if (isSolutionInHill(*subpop->best(env), hill)) {
						Real dis = optimum.variableDistance(*subpop->best(), env);
						if (nearest == nullptr || dis < min_dis) {
							min_dis = dis;
							nearest = subpop->best();
						}
					}
				}
				if (nearest != nullptr && nearest->objectiveDistance(optimum) < epsilon) {
					candidates.push_back(nearest);
				}
			}
		}
		else {
			for (auto &subpop : m_subpops) {
				candidates.push_back(subpop->best(env));
			}
		}
	}
}
