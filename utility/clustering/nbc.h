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
* NBC (Nearest-Better Clustering)
* Preuss, Mike. 
* "Niching the CMA-ES via nearest-better clustering." 
* Proceedings of the 12th annual conference companion on Genetic and evolutionary computation. 2010.
*-------------------------------------------------------------------------------
* Implemented by Junchen Wang (wangjunchen.chris@gmail.com) at 2020/9/22
*********************************************************************************/

#ifndef OFEC_NBC_H
#define OFEC_NBC_H

#include <algorithm>
#include <vector>
#include <list>
#include "../../core/problem/solution.h"

namespace ofec {
	class NBC {
	public:
		enum UpdateNBD { kByDistMat, kByKDTree };
		enum CalThrshld { kByMean, kByOutlier };

	private:
		Real m_phi;
		size_t m_N;                 // size of data
		UpdateNBD m_update_nbd;
		CalThrshld m_cal_thrshld;
		bool m_use_rule2;
		std::vector<const SolutionBase*> m_data;
		std::vector<std::vector<Real>> m_distance;
		std::vector<int> m_graph;   // target node, if -1 means no target node
		std::vector<size_t> m_cluster_centers;
		std::vector<std::vector<size_t>> m_clusters;
		std::vector<Real> m_nb_dis_1; // NBDs of solutions excluding the best
		std::vector<Real> m_nb_dis_2; // NBDs of all solutions 
		Real m_mean, m_stddev;

	public:
		NBC(Real phi = 2, UpdateNBD update_nbd = kByDistMat, CalThrshld cal_thrshld = kByMean, bool use_rule2 = false);
		void setPhi(Real phi) { m_phi = phi; }
		void setCalThrshld(CalThrshld ct) { m_cal_thrshld = ct; }
		void setData(const std::vector<const SolutionBase*> &sols);
		template<typename TVariable>
		void setData(const std::vector<const Solution<TVariable>*> &sols);
		template<typename TPopulation>
		void setData(const TPopulation &pop);
		template<typename TSolution>
		void setData(const std::vector<TSolution> &sols);
		template<typename TSolution>
		void setData(const std::vector<std::unique_ptr<TSolution>> &sols);
		void clustering(Environment *env);
		void clustering(size_t num, Environment *env);
		const std::vector<size_t>& clusterCenters() const { return m_cluster_centers; }
		const std::vector<std::vector<size_t>>& clusters() const { return m_clusters; }
		CalThrshld calThrshld() const { return m_cal_thrshld; }
		const std::vector<Real>& nearestBetterDis() const { return m_nb_dis_1; }
		const std::vector<Real>& nearestBetterDis2() const { return m_nb_dis_2; }
		const std::vector<int>& graph() const { return m_graph; }
		Real mean() const { return m_mean; }
		Real stddev() const { return m_stddev; }

		void updateNbDistByDistMat(Environment *env);
		void updateNbDistByKDTree(Environment *env);
		void cutEdgesInGraph();
		void cutEdgesInGraph(size_t num);
		void cutEdgesInGraph(Real thrshld);
		void cutEdgesInGraph(std::vector<size_t> &ids_centers);
		void updateClusters();
	};
	
	template<typename TSolution>
	void NBC::setData(const std::vector<TSolution> &sols) {
		m_N = sols.size();
		m_graph.resize(m_N);
		if (m_data.size() != m_N)
			m_data.resize(m_N);
		for (size_t i = 0; i < m_N; ++i)
			m_data[i] = &sols[i];
		if (!m_clusters.empty())
			m_clusters.clear();
	}

	template<typename TSolution>
	void NBC::setData(const std::vector<std::unique_ptr<TSolution>> &sols) {
		m_N = sols.size();
		m_graph.resize(m_N);
		if (m_data.size() != m_N)
			m_data.resize(m_N);
		for (size_t i = 0; i < m_N; ++i)
			m_data[i] = sols[i].get();
		if (!m_clusters.empty())
			m_clusters.clear();
	}

	template<typename TPopulation>
	void NBC::setData(const TPopulation& pop) {
		m_N = pop.size();
		m_graph.resize(m_N);
		if (m_data.size() != m_N)
			m_data.resize(m_N);
		for (size_t i = 0; i < m_N; ++i)
			m_data[i] = &pop[i];
		if (!m_clusters.empty())
			m_clusters.clear();
	}

	template<typename TVariable>
	void NBC::setData(const std::vector<const Solution<TVariable>*> &sols) {
		m_N = sols.size();
		m_graph.resize(m_N);
		if (m_data.size() != m_N)
			m_data.resize(m_N);
		for (size_t i = 0; i < m_N; ++i)
			m_data[i] = sols[i];
		if (!m_clusters.empty())
			m_clusters.clear();
	}
}

#endif // !OFEC_NBC_H
