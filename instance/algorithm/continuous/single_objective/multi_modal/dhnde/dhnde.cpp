#include "dhnde.h"
#include "../../../../../../utility/functional.h"
#include <numeric>

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"


namespace ofec {
	void DHNDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
		m_input_parameters.add("cluster size", new RangedSizeT(m_m, 5, 1000, 5));
		m_input_parameters.add("without crowding", new Bool(m_without_crowding, false));
	}

	void DHNDE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_scaling_factor = 0.9;
		m_Cr = 0.3;
		m_d_min = 0.1;
		m_f_min = 1.0e-5;
		m_max_io = 2 * m_pop_size;
		m_max_s = 5 * m_pop_size;
		m_max_T = 20;
		m_ms = de::MutateStrategy::kCurrentToBest1;
	}

	void DHNDE::run_(Environment *env) {
		m_A_os.clear();
		m_A_io.clear();
		m_pop.resize(m_pop_size, env);
		m_pop.scalingFactor() = m_scaling_factor;
		m_pop.crossoverRate() = m_Cr;
		m_pop.mutationStrategy() = m_ms;
		
		m_pop.initialize(env, m_random.get());
		m_pop.evaluate(env);
#ifdef OFEC_DATUM_MULTI_POP_H
		g_multi_pop.pops.clear();
		g_multi_pop.pops.resize(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			g_multi_pop.pops[0].push_back(&m_pop[i]);
		}
		datumUpdated(env, g_multi_pop);
#endif
		bool flag;
		while (!terminating()) {
			CDE_Aio(env);
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < m_pop.size(); ++i) {
				g_multi_pop.pops[0].push_back(&m_pop[i]);
			}
			datumUpdated(env, g_multi_pop);
#endif
			flag = 1;
			while (m_A_io.size() >= m_max_io && !terminating()) {
				if (flag == 1) {
					for (auto &s : m_A_io) {
						m_pop.append(IndDHNDE(s));
					}
					flag = 0;
#ifdef OFEC_DATUM_MULTI_POP_H
					g_multi_pop.pops.clear();
					g_multi_pop.pops.resize(1);
					for (size_t i = 0; i < m_pop.size(); ++i) {
						g_multi_pop.pops[0].push_back(&m_pop[i]);
					}
					datumUpdated(env, g_multi_pop);
#endif
				}
				INSDE(env);
#ifdef OFEC_DATUM_MULTI_POP_H
				g_multi_pop.pops.clear();
				g_multi_pop.pops.resize(1);
				for (size_t i = 0; i < m_pop.size(); ++i) {
					g_multi_pop.pops[0].push_back(&m_pop[i]);
				}
				datumUpdated(env, g_multi_pop);
#endif
				if (m_pop.size() < m_pop_size) {
					m_A_io.clear();
					for (size_t i = m_pop.size(); i < m_pop_size; ++i) {
						m_pop.append(IndDHNDE(env));
						m_pop.back().initialize(env, m_random.get());
						m_pop.back().evaluate(env);
					}
#ifdef OFEC_DATUM_MULTI_POP_H
					g_multi_pop.pops.clear();
					g_multi_pop.pops.resize(1);
					for (size_t i = 0; i < m_pop.size(); ++i) {
						g_multi_pop.pops[0].push_back(&m_pop[i]);
					}
					datumUpdated(env, g_multi_pop);
#endif
				}
			}
		}
	}

	void DHNDE::archiveUpdate(const Solution<> &x, Real d, Environment *env) {
		if (m_A_os.empty()) {
			m_A_os.push_back(x);
		}
		else {
			auto x_c = m_A_os.begin(), iter = m_A_os.begin();
			while (++iter != m_A_os.end()) {
				if (x_c->variableDistance(x, env) > iter->variableDistance(x, env)) {
					x_c = iter;
				}
			}
			if (x_c->variableDistance(x, env) < d) {
				if (dominate(x, *x_c, env->problem()->optimizeMode())) {
					(*x_c) = x;
				}
			}
			else {
				if (m_A_os.size() < m_max_s) {
					m_A_os.push_back(x);
				}
				else {
					if (dominate(x, *x_c, env->problem()->optimizeMode())) {
						(*x_c) = x;
					}
				}
			}
		}
	}
	
	void DHNDE::CDE_Aio(Environment *env) {
		m_pop.updateBest(env);
		for (size_t i = 0; i < m_pop_size; ++i) {
			m_pop.mutate(i, m_random.get(), env);
			m_pop.recombine(i, m_random.get(), env);
			m_pop[i].trial().evaluate(env);
		}
		for (size_t i = 0; i < m_pop_size; ++i) {
			auto &u_i = m_pop[i].trial();
			decltype(m_pop.begin()) y;
			if (m_without_crowding) {
				y = m_pop.begin() + i;
			}
			else {
				y = m_pop.begin();
				auto it = m_pop.begin();
				while (++it != m_pop.end()) {
					if ((*y)->variableDistance(u_i, env) > (*it)->variableDistance(u_i, env)) {
						y = it;
					}
				}
			}
			if (dominate(u_i, **y, env->problem()->optimizeMode())) {
				**y = u_i;
			}
			else {
				m_A_io.push_back(u_i);
			}
		}
	}
	
	void DHNDE::INSDE(Environment *env) {
		std::vector<size_t> order(m_pop.size());
		std::iota(order.begin(), order.end(), 0);
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			std::sort(order.begin(), order.end(),
				[this](size_t i, size_t j) {
					return m_pop[i].objective(0) < m_pop[j].objective(0);
				}
			);
		}
		else {
			std::sort(order.begin(), order.end(),
				[this](size_t i, size_t j) {
					return m_pop[i].objective(0) > m_pop[j].objective(0);
				}
			);
		}
		std::list<size_t> sort;
		for (size_t i = 0; i < m_pop.size(); ++i) {
			sort.push_back(order[i]);
		}
		std::vector<bool> processed(m_pop.size(), false);
		std::vector<PopulationDE<IndDHNDE>> subpop;
		for (size_t i = 0; i < m_pop.size() / m_m; ++i) {
			subpop.resize(subpop.size() + 1);
			auto &new_pop = subpop.back();
			size_t first = sort.front();
			new_pop.append(m_pop[first]);
			processed[first] = true;
			sort.pop_front();
			std::multimap<Real, size_t> msd;// map sorted by distance
			for (size_t ind : sort) {
				msd.emplace(m_pop[first].variableDistance(m_pop[ind], env), ind);
			}
			for (auto it = msd.begin(); it != msd.end() && new_pop.size() < m_m; ++it) {
				new_pop.append(m_pop[it->second]);
				processed[it->second] = true;
			}
			for (auto it = sort.begin(); it != sort.end();) {
				if (processed[*it]) {
					it = sort.erase(it);
				}
				else {
					it++;
				}
			}
			new_pop.scalingFactor() = m_scaling_factor;
			new_pop.crossoverRate() = m_Cr;
			new_pop.mutationStrategy() = m_ms;
		}
		for (size_t i = 0; i < subpop.size(); ++i) {
			subpop[i].updateBest(env);
			for (size_t j = 0; j < subpop[i].size(); j++) {
				subpop[i].mutate(j, m_random.get(), env);
				subpop[i].recombine(j, m_random.get(), env);
				subpop[i][j].trial().evaluate(env);
			}
		}




		for (size_t i = 0; i < subpop.size(); ++i) {
			for (size_t j = 0; j < subpop[i].size(); j++) {
				if (dominate(subpop[i][j].trial(), subpop[i][j], env->problem()->optimizeMode())) {
					subpop[i][j].resetCounter();
					subpop[i][j] = subpop[i][j].trial();
				}
				else {
					subpop[i][j].increaseCounter();
				}
			}
		}
		m_pop.clear();
		for (size_t i = 0; i < subpop.size(); ++i) {
			for (size_t j = 0; j < subpop[i].size(); j++) {
				bool add = true;
				for (size_t k = 0; k < m_pop.size(); ++k) {
					Real vd = m_pop[k].variableDistance(subpop[i][j], env);
					bool df = dominate(m_pop[k], subpop[i][j], env->problem()->optimizeMode());
					Real od = m_pop[k].objectiveDistance(subpop[i][j]);
					if (vd < m_d_min && df && od < m_f_min) {
						add = false;
						break;
					}
				}
				if (add) {
					m_pop.append(subpop[i][j]);
				}
			}
		}
		for (auto it = m_pop.begin(); it != m_pop.end();) {
			if ((*it)->counter() >= m_max_T) {
				archiveUpdate((**it), 0.01, env);
				it = m_pop.remove(it);
			}
			else {
				it++;
			}
		}
	}
}