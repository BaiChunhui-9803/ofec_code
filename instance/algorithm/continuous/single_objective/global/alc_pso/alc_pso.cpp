#include "alc_pso.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
	void ALC_PSO::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_swarm_size, 5, 1000, 20));
	}

	void ALC_PSO::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_weight = 0.4;
		m_accelerator1 = m_accelerator2 = 2.0;
		m_initial_lifespan = 10;
		m_T = 2;
		m_pro = 1.0 / env->problem()->numberVariables();
	}

	void ALC_PSO::run_(Environment *env) {
		m_swarm.resize(m_swarm_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_swarm.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_swarm[i].pbest());
		}
#endif
		m_challenger = m_swarm[0];
		m_old.resize(m_swarm_size, m_swarm[0]);
		/* Step 1 Initialization */
		m_swarm.initialize(env, m_random.get());
		m_swarm.evaluate(env);
		m_swarm.initPbest(env);
		m_swarm.initVelocityMax(env, m_random.get());
		m_leader = m_swarm.best(env)->pbest();
		m_age = 0;
		m_lifespan = m_initial_lifespan;
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
		while (!terminating()) {
			/* Step 2 Velocity and Position Updating */
			for (size_t i = 0; i < m_swarm_size; ++i) {
				m_swarm[i].nextVelocity(&m_leader, m_weight, m_accelerator1, m_accelerator2, m_random.get());
				m_swarm[i].move();
			}
			/* Step 3 Updating pBest and Leader */
			for (size_t i = 0; i < m_swarm_size; ++i) {
				m_swarm[i].evaluate(env);
				if (dominate(m_swarm[i], m_swarm[i].pbest(), env->problem()->optimizeMode())) {
					m_swarm[i].pbest() = m_swarm[i];
				}
				if (dominate(m_swarm[i], m_leader, env->problem()->optimizeMode())) {
					m_leader = m_swarm[i];
				}
			}
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
			/* Step 4 Lifespan Control */
			adjustLifespan();
			m_age++;
			if (m_age >= m_lifespan) {
				/* Step 5 Generating a Challenger */
				generateChallenger(env);
				/* Step 6 Evaluating the Challenger */
				evaluateChallenger(env);
			}
		}
	}

	void ALC_PSO::adjustLifespan() {
		std::array<Real, 3> cur_obj;
		cur_obj[0] = m_swarm.best()->pbest().objective(0);
		cur_obj[1] = 0;
		for (size_t i = 0; i < m_swarm_size; ++i) {
			cur_obj[1] += m_swarm[i].pbest().objective(0);
		}
		cur_obj[2] = m_leader.objective(0);
		if (m_age != 0) {
			if (cur_obj[0] != m_pre_obj[0]) {
				m_lifespan += 2;
			}
			else if (cur_obj[1] != m_pre_obj[1]) {
				m_lifespan++;
			}
			else if (cur_obj[2] == m_pre_obj[2]) {
				m_lifespan--;
			}
		}
		m_pre_obj[0] = cur_obj[0];
		m_pre_obj[1] = cur_obj[1];
		m_pre_obj[2] = cur_obj[2];
	}

	void ALC_PSO::generateChallenger(Environment *env) {
		size_t num_vars = CAST_CONOP(env->problem())->numberVariables();
		size_t count = 0;
		for (size_t j = 0; j < num_vars; ++j) {
			if (m_random->uniform.next() < m_pro) {
				Real L_j = CAST_CONOP(env->problem())->range(j).first;
				Real U_j = CAST_CONOP(env->problem())->range(j).second;
				m_challenger.variable()[j] = m_random->uniform.nextNonStd(L_j, U_j);
				count++;
			}
			else {
				m_challenger.variable()[j] = m_leader.variable()[j];
			}
		}
		if (count == 0) {
			size_t ran = m_random->uniform.nextNonStd<size_t>(0, num_vars);
			Real L_ran = CAST_CONOP(env->problem())->range(ran).first;
			Real U_ran = CAST_CONOP(env->problem())->range(ran).second;
			m_challenger.variable()[ran] = m_random->uniform.nextNonStd(L_ran, U_ran);
		}
	}

	void ALC_PSO::evaluateChallenger(Environment *env) {
		for (size_t i = 0; i < m_swarm_size; ++i) {
			m_old[i] = m_swarm[i];
		}
		for (size_t iters = 0; iters < m_T && !terminating(); ++iters) {
			for (size_t i = 0; i < m_swarm_size; ++i) {
				m_swarm[i].nextVelocity(&m_challenger, m_weight, m_accelerator1, m_accelerator2, m_random.get());
				m_swarm[i].move();
				m_swarm[i].evaluate(env);
			}
			bool improved = false;
			for (size_t i = 0; i < m_swarm_size; ++i) {
				if (dominate(m_swarm[i], m_swarm[i].pbest(), env->problem()->optimizeMode())) {
					m_swarm[i].pbest() = m_swarm[i];
					improved = true;
				}
				if (dominate(m_swarm[i], m_challenger, env->problem()->optimizeMode())) {
					m_challenger = m_swarm[i];
				}
			}
			if (improved) {
				m_leader = m_challenger;
				m_age = 0;
				m_lifespan = m_initial_lifespan;
				return;
			}
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif
		}
		for (size_t i = 0; i < m_swarm_size; ++i) {
			m_swarm[i] = m_old[i];
		}
		m_age = m_lifespan - 1;
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif
	}
}
