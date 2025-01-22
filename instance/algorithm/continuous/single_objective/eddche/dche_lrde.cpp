#include "dche_lrde.h"

#ifdef OFEC_PLAYBACK
//#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

namespace ofec {
	void DCHE_LRDE::addInputParameters() {

	}

	void DCHE_LRDE::initialize_(Environment *env) {
		DCHE::initialize_(env);
		LRDE::initialize_(env);
	}

	void DCHE_LRDE::run_(Environment *env) {
		m_pop.resize(m_pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop.at(i));
		}
#endif
		m_archive.clear();
		m_best.reset();
		Hill *cur_hill = nullptr;
		cur_hill = (*m_random->uniform.nextElem(m_hills.begin(), m_hills.end())).get();
		initializePopulationInHill(cur_hill, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		int rf = kNormalEval;
		bool best_improved;
		size_t counter_stagnant = 0;
		while (!terminating()) { 
			best_improved = false;
			for (size_t i = 0; i < m_pop.size(); ++i) {
				do {
					m_pop.mutate(i, m_random.get(), env);
					m_pop.recombine(i, m_random.get(), env);
				} while (reject(m_pop[i].trial(), env));
				rf = m_pop[i].select(env);
				if (rf != kNormalEval) {
					break;
				}
				if (m_pop[i].isImproved() && (!m_best || dominate(m_pop[i], *m_best, env->problem()->optimizeMode()))) {
					m_best.reset(new Solution<>(m_pop[i]));
					best_improved = true;
				}
			}
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
			if (rf == kNormalEval) {
				if (!best_improved) {
					counter_stagnant++;
				}
				else {
					counter_stagnant = 0;
				}
				if (counter_stagnant > m_max_stagant_iters) {
					updateArchive(env);
					updateSeeds(env);
					updateHills(env);
					cur_hill = (*m_random->uniform.nextElem(m_hills.begin(), m_hills.end())).get();
					initializePopulationInHill(cur_hill, env);
#ifdef OFEC_DATUM_MULTI_POP_H
					datumUpdated(env, g_multi_pop);
#endif
				}
			}
		}
	}

	bool DCHE_LRDE::reject(const Solution<> &sol, Environment *env) const {
		for (auto &s : m_archive) {
			size_t id_ssp = m_sp_tree->getRegionIdx(s.variable().vector());
			auto &hills = m_subspaces.at(id_ssp)->assignedHills();
			Real ratio = pow(hills.front()->volume(), 1.0 / env->problem()->numberVariables());
			if (s.variableDistance(sol, env) <
				m_repulsion_radius * env->problem()->maximumVariableDistance() * ratio)
			{
				if (m_random->uniform.next() < m_repulsion_probabilty) {
					return true;
				}
			}
			std::string str;
		}
		return false;
	}

	void DCHE_LRDE::initializePopulationInHill(const Hill *hill, Environment *env) {
		for (size_t i = 0; i < m_pop.size(); ++i) {
			do {
				randomSolutionInHill(m_pop[i], hill);
			} while (reject(m_pop[i], env));
		}
		m_pop.evaluate(env);
		m_best.reset();
	}

	void DCHE_LRDE::identifyCandidates(std::list<const Solution<>*> &candidates, Environment *env) {
		if (m_omniscient_seeds) {
			const Real epsilon = 1e-2;
			for (size_t k = 0; k < env->problem()->numberOptimumSolutions(); ++k) {
				auto &optimum = CAST_CONOP(env->problem())->optima()->solution(k);
				size_t id_ssp = m_sp_tree->getRegionIdx(optimum.variable().vector());
				auto hill = m_subspaces.at(id_ssp)->assignedHills().front();
				Solution<> *nearest = nullptr;
				Real min_dis = std::numeric_limits<Real>::max();
				for (auto &a : m_archive) {
					if (isSolutionInHill(a, hill)) {
						Real dis = optimum.variableDistance(a, env);
						if (nearest == nullptr || dis < min_dis) {
							min_dis = dis;
							nearest = &a;
						}
					}
				}
				if (nearest != nullptr && nearest->objectiveDistance(optimum) < epsilon) {
					candidates.push_back(nearest);
				}
			}
		}
		else {
			for (auto &a : m_archive) {
				candidates.push_back(&a);
			}
		}
	}
}
