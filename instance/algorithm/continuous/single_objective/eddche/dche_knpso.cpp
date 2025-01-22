#include "dche_knpso.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS


namespace ofec {
	void DCHE_kNPSO::addInputParameters() {
		m_input_parameters.add("maximum stagnant iterations", new RangedSizeT(m_max_stagant_iters, 5, 100, 10));
	}

	void DCHE_kNPSO::initialize_(Environment *env) {
		DCHE::initialize_(env);
		kNPSO::initialize_(env);
		TracerPopNBD::initialize_(env);
	}

	void DCHE_kNPSO::run_(Environment *env) {
		m_swarm.setNumNeighbors(m_num_nbrs);
		m_swarm.resize(m_swarm_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_swarm.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_swarm.at(i).pbest());
	}
#endif
		Hill *cur_hill = nullptr;
		cur_hill = (*m_random->uniform.nextElem(m_hills.begin(), m_hills.end())).get();
		initializePopulationInHill(cur_hill, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		size_t counter_stagnant = 0;
		while (!terminating()) {
			int rf = m_swarm.evolve(env, m_random.get());
			addLatestPop(m_swarm, env);
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
			if (rf == kNormalEval) {
				bool flag_stagnant = true;
				for (size_t i = 0; i < m_swarm.size(); ++i) {
					if (m_swarm[i].isImproved()) {
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
#ifdef OFEC_DATUM_MULTI_POP_H
					datumUpdated(env, g_multi_pop);
#endif
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

	void DCHE_kNPSO::initializePopulationInHill(const Hill *hill, Environment *env) {
		VariableVector<> rnd_var_1(env->problem()->numberVariables()), rnd_var_2(rnd_var_1);
		for (size_t i = 0; i < m_swarm.size(); ++i) {
			randomSolutionInHill(m_swarm[i], hill);
			randomVariablesInHill(rnd_var_1, hill);
			randomVariablesInHill(rnd_var_2, hill);
			for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
				m_swarm[i].velocity()[j] = (rnd_var_1[j] - rnd_var_2[j]) / 4;
			}
			m_swarm[i].initVelocityMax(env, nullptr);
		}
		m_swarm.evaluate(env);
		m_swarm.initPbest(env);
		clearEvoPop();
		addLatestPop(m_swarm, env);
	}

	void DCHE_kNPSO::identifyCandidates(std::list<const Solution<>*> &candidates, Environment *env) {
		if (m_omniscient_seeds) {
			const Real epsilon = 1e-2;
			for (size_t k = 0; k < env->problem()->numberOptimumSolutions(); ++k) {
				auto &optimum = CAST_CONOP(env->problem())->optima()->solution(k);
				size_t id_ssp = m_sp_tree->getRegionIdx(optimum.variable().vector());
				auto hill = m_subspaces.at(id_ssp)->assignedHills().front();
				Solution<> *nearest = nullptr;
				Real min_dis = std::numeric_limits<Real>::max();
				for (size_t i = 0; i < m_swarm.size(); ++i) {
					if (isSolutionInHill(m_swarm[i].pbest(), hill)) {
						Real dis = optimum.variableDistance(m_swarm[i].pbest(), env);
						if (nearest == nullptr || dis < min_dis) {
							min_dis = dis;
							nearest = &m_swarm[i].pbest();
						}
					}
				}
				if (nearest != nullptr && nearest->objectiveDistance(optimum) < epsilon) {
					candidates.push_back(nearest);
				}
			}
		}
		else {
			identifyOutliers(candidates, env);
		}
	}
}
