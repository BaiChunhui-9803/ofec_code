/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* FDC: Fitness Density Clustering
* Proposed by Yiya Diao
*********************************************************************************/

#ifndef OFEC_FDC_H
#define OFEC_FDC_H

#include "../../core/definition.h"	
#include <vector>
#include <memory>
#include <algorithm>

namespace ofec {
	template<typename TInd>
	class FDC {
	protected:
		std::vector<Real> m_fitness;
		std::vector<std::vector<Real>> m_distance;
		std::vector<std::vector<size_t>> m_clusters;
		std::vector<std::pair<Real, Real>> m_fitness_distance;
		std::vector<int> m_belong;
		std::vector<int> m_belong2secondBest;
		std::vector<size_t> m_best_order;
		std::vector<size_t> m_center_order;
		std::vector<size_t> m_cluster_num;
		std::vector<Real> m_cut_ratios;
		Real m_dist_threshold;
		bool m_obj_dist_calulated;

	public:
		FDC(const Population<TInd> &pop, Problem *pro) :
			m_fitness(pop.size()),
			m_distance(pop.size()),
			m_fitness_distance(pop.size()),
			m_belong(pop.size(), -1),
			m_best_order(pop.size()),
			m_belong2secondBest(pop.size(), -1),
			m_center_order(pop.size()),
			m_cluster_num(pop.size(), 0),
			m_obj_dist_calulated(false)
		{
			for (int pop_idx(0); pop_idx < pop.size(); ++pop_idx) {
				m_fitness[pop_idx] = pop[pop_idx].fitness();
				m_fitness_distance[pop_idx].first = m_fitness[pop_idx];
				m_distance[pop_idx].resize(pop.size());
				m_distance[pop_idx][pop_idx] = 0;
				for (int other_idx(0); other_idx < pop_idx; ++other_idx) {
					m_distance[other_idx][pop_idx] = pop[pop_idx].variableDistance(pop[other_idx], pro);
					m_distance[pop_idx][other_idx] = pop[pop_idx].variableDistance(pop[other_idx], pro);
				}
			}
			for (int i = 0; i < m_best_order.size(); ++i) {
				m_best_order[i] = i;
			}
			for (int i(0); i < m_center_order.size(); ++i) {
				m_center_order[i] = i;
			}
		}

		FDC(const std::vector<std::unique_ptr<TInd>> &pop, Problem *pro) :
			m_fitness(pop.size()),
			m_distance(pop.size()),
			m_fitness_distance(pop.size()),
			m_belong(pop.size(), -1),
			m_best_order(pop.size()),
			m_belong2secondBest(pop.size(), -1),
			m_center_order(pop.size()),
			m_cluster_num(pop.size(), 0),
			m_obj_dist_calulated(false)
		{
			for (int pop_idx(0); pop_idx < pop.size(); ++pop_idx) {
				m_fitness[pop_idx] = pop[pop_idx]->fitness();
				m_fitness_distance[pop_idx].first = m_fitness[pop_idx];
				m_distance[pop_idx].resize(pop.size());
				m_distance[pop_idx][pop_idx] = 0;
				for (int other_idx(0); other_idx < pop_idx; ++other_idx) {
					m_distance[other_idx][pop_idx] = pop[pop_idx]->variableDistance(*pop[other_idx], pro);
					m_distance[pop_idx][other_idx] = pop[pop_idx]->variableDistance(*pop[other_idx], pro);
				}
			}
			for (int i = 0; i < m_best_order.size(); ++i) {
				m_best_order[i] = i;
			}
			for (int i(0); i < m_center_order.size(); ++i) {
				m_center_order[i] = i;
			}
		}

		FDC(const std::vector<TInd> &pop, Problem *pro) :
			m_fitness(pop.size()),
			m_distance(pop.size()),
			m_fitness_distance(pop.size()),
			m_belong(pop.size(), -1),
			m_best_order(pop.size()),
			m_belong2secondBest(pop.size(), -1),
			m_center_order(pop.size()),
			m_cluster_num(pop.size(), 0),
			m_obj_dist_calulated(false)
		{
			for (int pop_idx(0); pop_idx < pop.size(); ++pop_idx) {
				m_fitness[pop_idx] = pop[pop_idx].fitness();
				m_fitness_distance[pop_idx].first = m_fitness[pop_idx];
				m_distance[pop_idx].resize(pop.size());
				m_distance[pop_idx][pop_idx] = 0;
				for (int other_idx(0); other_idx < pop_idx; ++other_idx) {
					m_distance[other_idx][pop_idx] = pop[pop_idx].variableDistance(pop[other_idx], pro);
					m_distance[pop_idx][other_idx] = pop[pop_idx].variableDistance(pop[other_idx], pro);
				}
			}
			for (int i = 0; i < m_best_order.size(); ++i) {
				m_best_order[i] = i;
			}
			for (int i(0); i < m_center_order.size(); ++i) {
				m_center_order[i] = i;
			}

		}

		const std::vector<std::pair<Real, Real>>& getFitDist() {
			return m_fitness_distance;
		}

		const std::vector<int>& getBelong() {
			return m_belong;
		}

		const std::vector<int>& getBelong2Second() {
			return m_belong2secondBest;
		}

		const std::vector<size_t>& getBestOrder() {
			return m_best_order;
		}

		void calFitDist() {
			if (m_fitness.size() == 1) {
				m_fitness_distance.front().second = 1.0;
				return;
			}
			std::sort(m_best_order.begin(), m_best_order.end(), [&](int a, int b) {
				return m_fitness[a] > m_fitness[b];
				});
			std::vector<int> best_order_idx(m_best_order.size());
			for (int i(0); i < m_best_order.size(); ++i) {
				best_order_idx[m_best_order[i]] = i;
			}
			Real max_dis(0);
			for (size_t cur_idx(1); cur_idx < m_best_order.size(); ++cur_idx) {
				m_fitness_distance[m_best_order[cur_idx]].second = std::numeric_limits<Real>::max();
				m_belong2secondBest[m_best_order[cur_idx]] = m_belong[m_best_order[cur_idx]] = m_best_order[cur_idx];
				for (size_t parent_idx(0); parent_idx < cur_idx; ++parent_idx) {
					if (m_fitness_distance[m_best_order[cur_idx]].second > m_distance[m_best_order[cur_idx]][m_best_order[parent_idx]]) {
						m_fitness_distance[m_best_order[cur_idx]].second = m_distance[m_best_order[cur_idx]][m_best_order[parent_idx]];
						m_belong2secondBest[m_best_order[cur_idx]] = m_belong[m_best_order[cur_idx]];
						m_belong[m_best_order[cur_idx]] = m_best_order[parent_idx];
					}
				}
				max_dis = std::max(max_dis, m_fitness_distance[m_best_order[cur_idx]].second);
			}
 			m_fitness_distance[m_best_order.front()].second = max_dis;
			for (int pop_idx(m_best_order.size() - 1); pop_idx >= 1; --pop_idx) {
				++m_cluster_num[m_best_order[pop_idx]];
				m_cluster_num[m_belong[m_best_order[pop_idx]]] += m_cluster_num[m_best_order[pop_idx]];
			}
			++m_cluster_num[m_best_order.front()];
			std::sort(m_center_order.begin(), m_center_order.end(), [&](int a, int b) {
				if (m_fitness_distance[a].second == m_fitness_distance[b].second)
					return m_cluster_num[a] > m_cluster_num[b];
				else
					return m_fitness_distance[a].second > m_fitness_distance[b].second;

				});
			m_obj_dist_calulated = true;
		}

		void printInfo() const {
			std::cout << "\n\nsolutions_distance" << std::endl;
			for (auto &it : m_center_order) {
				std::cout << m_fitness_distance[it].second << "\t";
			}
			//int solutions_idx(0);
			//for (auto& it : m_best_order) {
			//	std::cout <<solutions_idx++<<":\t"<< "\t" << m_fitness_distance[it].first << "\t" << m_fitness_distance[it].second << "\t" << m_cluster_num[it] << "\t" << "\t";
			//}
			std::cout << std::endl << std::endl << std::endl;
			std::cout << "\n\nsolutions cluster number" << std::endl;
			for (auto &it : m_center_order) {
				std::cout << m_cluster_num[it] << "\t";
			}
			//int solutions_idx(0);
			//for (auto& it : m_best_order) {
			//	std::cout <<solutions_idx++<<":\t"<< "\t" << m_fitness_distance[it].first << "\t" << m_fitness_distance[it].second << "\t" << m_cluster_num[it] << "\t" << "\t";
			//}
			std::cout << std::endl << std::endl << std::endl;
			//cluster(30);
			std::cout << "\n\ncluster size" << std::endl;
			for (auto &it : m_clusters) {
				std::cout << it.size() << "\t";
			}
			//int solutions_idx(0);
			//for (auto& it : m_best_order) {
			//	std::cout <<solutions_idx++<<":\t"<< "\t" << m_fitness_distance[it].first << "\t" << m_fitness_distance[it].second << "\t" << m_cluster_num[it] << "\t" << "\t";
			//}
			std::cout << std::endl << std::endl << std::endl;
		}

		void clusterByDist(Real ratio_distance = 0.1) {
			if (!m_obj_dist_calulated)
				calFitDist();
			if (ratio_distance <= 0 || ratio_distance >= 1)
				ratio_distance = 0.1;
			std::vector<Real> dists(m_fitness_distance.size());
			for (size_t i = 0; i < dists.size(); i++)
				dists[i] = m_fitness_distance[i].second;
			std::sort(dists.begin(), dists.end());
			Real max_dist = dists.back();
			Real min_dist = dists.front();
			m_dist_threshold = min_dist + (max_dist - min_dist) * ratio_distance;
			std::vector<int> belong_cluster(m_fitness_distance.size(), -1);
			int cluster_id(0);
			for (auto &center_id : m_center_order) {
				if (m_fitness_distance[center_id].second > m_dist_threshold)
					belong_cluster[center_id] = cluster_id++;
			}
			m_clusters.resize(cluster_id);
			for (auto &pop_id : m_best_order) {
				if (belong_cluster[pop_id] == -1) {
					belong_cluster[pop_id] = belong_cluster[m_belong[pop_id]];
				}
				m_clusters[belong_cluster[pop_id]].push_back(pop_id);
			}
		}

		void clusterBySize(int min_popsize) {
			if (!m_obj_dist_calulated)
				calFitDist();
			std::vector<int> belong_cluster(m_fitness_distance.size(), -1);
			int cluster_id(0);
			bool before_dis(m_fitness_distance[m_best_order.front()].second);
			for (auto &center_id : m_center_order) {
				if (m_cluster_num[center_id] > min_popsize || before_dis == m_fitness_distance[center_id].second)
					belong_cluster[center_id] = cluster_id++;
				else
					break;
				before_dis = m_fitness_distance[center_id].second;
			}
			m_clusters.resize(cluster_id);
			for (auto &pop_id : m_best_order) {
				if (belong_cluster[pop_id] == -1) {
					belong_cluster[pop_id] = belong_cluster[m_belong[pop_id]];
				}
				m_clusters[belong_cluster[pop_id]].push_back(pop_id);
			}
		}

		const std::vector<std::vector<size_t>>& clusters() const {
			return m_clusters;
		}

		Real distThreshold() const {
			return m_dist_threshold;
		}
	};

}


#endif // !OFEC_FDC_H