#include "lips_pop.h"
#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	LI_particle::LI_particle(size_t num_obj, size_t num_con, size_t size_var) :
		Particle(num_obj, num_con, size_var),
		m_pos(size_var),
		m_rdm(size_var)
	{
		m_improved = true;
	}

	void LI_particle::nextVelocityByWeight(Real w) {
		for (size_t j = 0; j < variable().size(); j++) {
			m_vel[j] = w * (m_vel[j] + m_rdm[j] * (m_pos[j] - variable()[j]));
		}
	}

	SwarmLIP::SwarmLIP(size_t size_pop, Environment *env, int max_evals) :
		Swarm(size_pop, env),
		m_M(2),
		m_maximum_evalutions(max_evals),
		m_dis(size_pop) {}

	void SwarmLIP::resize(size_t size_pop, Environment *env, int max_evals) {
		Swarm::resize(size_pop, env);
		m_M = 2;
		m_maximum_evalutions = max_evals;
		m_dis.resize(size_pop);
	}

	void SwarmLIP::clear() {
		Swarm::clear();
		m_dis.clear();
	}

	int SwarmLIP::evolve(Environment *env, Random *rnd) {
		int rf = kNormalEval;
		// Update the value of M
		if (m_iteration * m_individuals.size() > m_maximum_evalutions)
			m_M = 5;
		else
			m_M = 2 + 3. * static_cast<Real>(m_iteration * m_individuals.size()) / static_cast<Real>(m_maximum_evalutions);
		// Breeding process
		for (int i = 0; i < m_individuals.size(); i++) {
			sortDistance(i, env);
			setBestPos(i, env, rnd);
			m_individuals[i]->nextVelocityByWeight(m_weight);   // lbest and c1, c2 are actually not used in LIPS
			m_individuals[i]->move();
			m_individuals[i]->clampVelocity(env, rnd);
			rf = m_individuals[i]->evaluate(env);
			m_individuals[i]->updatePBest(env);
			handleEvaluationTag(rf);
			if (rf != kNormalEval) break;
		}
		this->m_iteration++;
		return rf;
	}

	void SwarmLIP::setBestPos(int idx_ind, Environment *env, Random *rnd) {
		for (auto& j : m_individuals[idx_ind]->m_pos) j = 0;
		for (auto& j : m_individuals[idx_ind]->m_rdm) j = 0;
		for (auto& i : m_individuals[idx_ind]->m_nbr) {
			for (auto j = 0; j < env->problem()->numberVariables(); j++) {
				double rdom = rnd->uniform.next() * (4.1 / m_M);
				m_individuals[idx_ind]->m_pos[j] += m_individuals[i]->pbest().variable()[j] * rdom;
				m_individuals[idx_ind]->m_rdm[j] += rdom;
			}
		}
		for (auto j = 0; j < env->problem()->numberVariables(); ++j) {
			m_individuals[idx_ind]->m_pos[j] = m_individuals[idx_ind]->m_pos[j] / m_individuals[idx_ind]->m_rdm[j];
		}
	}

	void SwarmLIP::sortDistance(int idx_ind, Environment *env) {
		if (!m_use_history_nearest)
			m_dis[idx_ind].clear();
		for (size_t j = 0; j < m_individuals.size(); j++) {
			if (idx_ind == j) continue;
			if (m_individuals[idx_ind]->isImproved() || m_individuals[j]->isImproved()) {
				std::pair<Real, size_t> dis = std::make_pair(
					m_individuals[idx_ind]->pbest().variableDistance(m_individuals[j]->pbest(), env), j);
				auto it = m_dis[idx_ind].begin();
				while (it != m_dis[idx_ind].end() && it->first < dis.first) {
					it++;
				}
				if (m_dis[idx_ind].size() >= m_M - 1) {
					if (it != m_dis[idx_ind].end()) {
						m_dis[idx_ind].insert(it, dis);
						m_dis[idx_ind].pop_back();
					}
				}
				else {
					m_dis[idx_ind].insert(it, dis);
				}
			}
		}
		m_individuals[idx_ind]->m_nbr.resize(0);
		for (auto& i : m_dis[idx_ind])
			m_individuals[idx_ind]->m_nbr.push_back(i.second);
		m_individuals[idx_ind]->m_nbr.push_back(idx_ind);
	}
}