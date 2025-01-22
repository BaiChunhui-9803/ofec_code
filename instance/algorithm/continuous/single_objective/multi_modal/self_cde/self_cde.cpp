#include "self_cde.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void SelfCDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 10, 1000, 100));
		m_input_parameters.add("cluster size", new RangedSizeT(m_cluster_size, 2, 1000, 20));
		m_input_parameters.add("without crowding", new Bool(m_without_crowding, false));
	}

	void SelfCDE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		if (m_pop_size % m_cluster_size != 0) {
			throw Exception("The population size should be a multiple of the cluster size.");
		}
	}

	void SelfCDE::run_(Environment *env) {
		m_pop.resize(m_pop_size, env);
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
		Real C_rm = 0.5;
		while (!terminating()) {
			m_pop.updateBest(env);
			m_pop.updateWorst(env);
			clusteringPartition(env);
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(m_clusters.size());
			for (size_t k = 0; k < m_clusters.size(); ++k) {
				for (size_t i = 0; i < m_clusters[k].size(); ++i)
					g_multi_pop.pops[k].push_back(&m_pop[m_clusters[k][i]]);
			}
			datumUpdated(env, g_multi_pop);
#endif
			std::vector<size_t> r;
			std::vector<Real> CR(m_pop_size);
			for (auto &cluster : m_clusters) {
				for (size_t i : cluster) {
					m_pop.selectInCandidates(3, cluster, r, m_random.get());
					Real F = (m_pop[r[1]].objective(0) - m_pop[r[2]].objective(0))
						/ (m_pop.best()->objective(0) - m_pop.worst()->objective(0));
					m_pop[i].mutate(F, &m_pop[r[0]], &m_pop[r[1]], &m_pop[r[2]], env);
					CR[i] = m_random->normal.nextNonStd(C_rm, 0.1);
					CR[i] = CR[i] > 1.0 ? 1.0 : CR[i]; CR[i] = CR[i] < 0.0 ? 0.0 : CR[i];
					m_pop[i].recombine(CR[i], de::RecombineStrategy::kBinomial, m_random.get(), env);
				}
			}
			std::list<Real> S_Cr;
			for (size_t i = 0; i < m_pop_size; ++i) {
				m_pop[i].trial().evaluate(env);
				if (m_without_crowding) {
					if (dominate(m_pop[i].trial(), m_pop[i], env->problem()->optimizeMode())) {
						m_pop[i].solution() = m_pop[i].trial();
						S_Cr.push_back(CR[i]);
					}
				}
				else {
					size_t nearest = nearestInd(m_pop[i].trial(), env);
					if (dominate(m_pop[i].trial(), m_pop[nearest], env->problem()->optimizeMode())) {
						m_pop[nearest].solution() = m_pop[i].trial();
						S_Cr.push_back(CR[i]);
					}
				}
			}
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1);
			for (size_t i = 0; i < m_pop.size(); ++i) {
				g_multi_pop.pops[0].push_back(&m_pop[i]);
			}
			datumUpdated(env, g_multi_pop);
#endif
			if (S_Cr.empty())
				C_rm = 0.5;
			else {
				C_rm = 0.0;
				for (Real Cr : S_Cr) {
					C_rm += Cr;
				}
				C_rm /= S_Cr.size();
			}
		}
	}

	size_t SelfCDE::nearestInd(const Solution<> &s, Environment *env) {
		size_t nearest = 0;
		Real min_dis = s.variableDistance(m_pop[0], env);
		for (size_t j = 1; j < m_pop_size; ++j) {
			Real dis = s.variableDistance(m_pop[j], env);
			if (dis < min_dis) {
				min_dis = dis;
				nearest = j;
			}
		}
		return nearest;
	}

	void SelfCCDE::clusteringPartition(Environment *env) {
		VariableVector<> R(env->problem()->numberVariables());
		env->problem()->initializeVariables(R, m_random.get());
		std::set<size_t> P;
		for (size_t i = 0; i < m_pop_size; ++i) {
			P.insert(i);
		}
		m_clusters.clear();
		while (P.size() > m_cluster_size) {
			m_clusters.resize(m_clusters.size() + 1);
			size_t X_tilde = *P.begin();
			Real min_dis = env->problem()->variableDistance(R, m_pop[X_tilde].variableBase());
			for (size_t i : P) {
				Real dis = env->problem()->variableDistance(R, m_pop[i].variableBase());
				if (dis < min_dis) {
					X_tilde = i;
					min_dis = dis;
				}
			}
			P.erase(X_tilde);
			m_clusters.back().push_back(X_tilde);
			for (size_t n = 0; n < m_cluster_size - 1; ++n) {
				size_t nearest = *P.begin();
				min_dis = m_pop[X_tilde].variableDistance(m_pop[nearest], env);
				for (size_t i : P) {
					Real dis = m_pop[X_tilde].variableDistance(m_pop[i], env);
					if (dis < min_dis) {
						nearest = i;
						min_dis = dis;
					}
				}
				P.erase(nearest);
				m_clusters.back().push_back(nearest);
			}
		}
		m_clusters.resize(m_clusters.size() + 1);
		for (size_t i : P) {
			m_clusters.back().push_back(i);
		}
	}

	void SelfCSDE::clusteringPartition(Environment *env) {
		std::set<size_t> P;
		for (size_t i = 0; i < m_pop_size; ++i) {
			P.insert(i);
		}
		m_clusters.clear();
		while (P.size() > m_cluster_size) {
			m_clusters.resize(m_clusters.size() + 1);
			size_t best = *P.begin();
			for (size_t i : P) {
				if (dominate(m_pop[i], m_pop[best], env->problem()->optimizeMode())) {
					best = i;
				}
			}
			P.erase(best);
			m_clusters.back().push_back(best);
			for (size_t n = 0; n < m_cluster_size - 1; ++n) {
				size_t nearest = *P.begin();
				Real min_dis = m_pop[best].variableDistance(m_pop[nearest], env);
				for (size_t i : P) {
					Real dis = m_pop[best].variableDistance(m_pop[i], env);
					if (dis < min_dis) {
						nearest = i;
						min_dis = dis;
					}
				}
				P.erase(nearest);
				m_clusters.back().push_back(nearest);
			}
		}
		m_clusters.resize(m_clusters.size() + 1);
		for (size_t i : P) {
			m_clusters.back().push_back(i);
		}
	}
}
