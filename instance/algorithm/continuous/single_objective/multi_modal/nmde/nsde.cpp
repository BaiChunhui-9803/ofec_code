#include "nsde.h"
#include <algorithm>
#include <numeric>

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	PopNSDE::PopNSDE(size_t size_pop, size_t cluster_size, Environment *env) : 
		PopulationDE<>(size_pop, env),
		m_m(cluster_size - 1),
		m_dis((size_pop / (m_m + 1)) + 1), 
		m_seed((size_pop / (m_m + 1)) + 1),
		m_order_list(size_pop)
	{
		std::iota(m_order_list.begin(), m_order_list.end(), 0);
	}

	void PopNSDE::selectSubpop(Environment *env) {
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			std::sort(m_order_list.begin(), m_order_list.end(),
				[this](int i, int j) {
					return m_individuals[i]->objective(0) < m_individuals[j]->objective(0);
				}
			);
		}
		else {
			std::sort(m_order_list.begin(), m_order_list.end(),
				[this](int i, int j) {
					return m_individuals[i]->objective(0) > m_individuals[j]->objective(0);
				}
			);
		}
		int idx_species = 0;
		while (m_order_list.size() > 0) {
			m_seed[idx_species] = m_order_list[0];
			//m_dis[idx_species].emplace_back(std::make_pair(0, m_order_list[0]));
			for (size_t j = 0; j < m_order_list.size(); ++j) {
				std::pair<Real, int> dis = std::make_pair(m_individuals[m_seed[idx_species]]->variableDistance(*m_individuals[m_order_list[j]], env), m_order_list[j]);
				auto it = m_dis[idx_species].begin();
				while (it != m_dis[idx_species].end() && it->first < dis.first) {
					it++;
				}
				if (m_dis[idx_species].size() >= m_m + 1) {
					if (it != m_dis[idx_species].end()) {
						m_dis[idx_species].insert(it, dis);
						m_dis[idx_species].pop_back();
					}
				}
				else {
					m_dis[idx_species].insert(it, dis);
				}
			}
			for (auto it = m_dis[idx_species].begin(); it != m_dis[idx_species].end(); it++) {
				auto it_ = std::find(m_order_list.begin(), m_order_list.end(), it->second);
				m_order_list.erase(it_);
			}
			++idx_species;
		}
	}

	int PopNSDE::evolve(Environment *env, Random *rnd) {
		int tag = kNormalEval;
		Real u, l;
		for (size_t i = 0; i < m_dis.size(); ++i) {
			std::vector<size_t> cluster;
			for (auto &p : m_dis[i]) {
				cluster.push_back(p.second);
			}
			for (auto it = m_dis[i].begin(); it != m_dis[i].end(); it++) {
				std::vector<size_t> ridx(3);
				selectInCandidates(3, cluster, ridx, rnd);
				m_individuals[it->second]->mutate(
					m_scaling_factor, 
					m_individuals[ridx[0]].get(),
					m_individuals[ridx[1]].get(),
					m_individuals[ridx[2]].get(),
					env
				);
				m_individuals[it->second]->recombine(m_crossover_rate, m_recombine_strategy, rnd, env);
				tag = m_individuals[it->second]->select(env);
				//if (m_seed[i] != it->second) {
				//	if (m_individuals[m_seed[i]]->objective(0) == m_individuals[it->second]->objective(0)) {
				//		m_individuals[it->second]->initialize(env, rnd);
				//		m_individuals[it->second]->evaluate(env);
				//	}
				//	else if (dominate(*m_individuals[it->second], *m_individuals[m_seed[i]], env->problem()->optimizeMode())) {
				//		*m_individuals[m_seed[i]] = *m_individuals[it->second];
				//	}
				//}
				if (tag != kNormalEval) break;
			}
			if (tag != kNormalEval) break;
		}
		++m_iteration;
		return tag;
	}

	void NSDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 10, 1000, 20));
		m_input_parameters.add("cluster size", new RangedSizeT(m_cluster_size, 5, 1000, 5));

	}

	void NSDE::initPop(size_t pop_size, Environment *env) {
		m_pop.reset(new PopNSDE(pop_size, m_cluster_size, env));
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop->at(i));
		}
#endif
		m_pop->initialize(env, m_random.get());
		m_pop->evaluate(env);
		m_pop->scalingFactor() = 0.9;
		m_pop->crossoverRate() = 0.1;
		m_pop->selectSubpop(env);
	}

	void NSDE::run_(Environment *env) {
		initPop(m_pop_size, env);
#ifdef OFEC_DATUM_MULTI_POP_H
		datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		while (!terminating()) {
			m_pop->evolve(env, m_random.get());
#ifdef OFEC_DATUM_MULTI_POP_H
			datumUpdated(env, g_multi_pop);
#endif // OFEC_DATUM_MULTI_POP_H
		}
	}
}


