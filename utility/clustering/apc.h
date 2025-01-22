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
* APC (affinity propagation clustering)
* Frey, Brendan J., and Delbert Dueck. 
* "Clustering by Passing Messages Between Data Points." 
* Science 315.5814 (2007): 972-976.
*-------------------------------------------------------------------------------
* Constructed by Junchen Wang (wangjunchen.chris@gmail.com) at 2020/5/13
*********************************************************************************/

#ifndef OFEC_APC_H
#define OEFC_APC_H

#include <algorithm>
#include <vector>
#include <set>
#include <numeric>
#include <fstream>
#include "../../core/algorithm/population.h"

namespace ofec {
	class APC {
		const Real m_lambda;//damping factor, [0,1], default 0.5
		const size_t m_Mits, m_Cits;//the last iterations of decisions did not change, 10 
		size_t m_N; // Size of data
		std::vector<std::vector<Real>> m_responsibility;
		std::vector<std::vector<Real>> m_availability;
		std::vector<std::vector<Real>> m_similarity;
		std::vector<const SolutionBase*> m_data;
		std::vector<std::vector<size_t>> m_clusters;
		std::vector<bool> m_is_exemplar;

	public:
		APC(Real lambda, size_t Mits, size_t Cits) : m_lambda(lambda), m_Mits(Mits), m_Cits(Cits) {}
		template<typename TInd>
		void updateData(const Population<TInd>& pop);
		template<typename TInd>
		void updateData(const std::vector<TInd>& inds);
		void clustering(Problem *pro);
		const std::vector<std::vector<size_t>>& clusters() const { return m_clusters; }

	protected:
		void updateSimilarity(Problem *pro);
		void updateResponsibility();
		void updateAvailability();
		void updateClusters();
	};

	template<typename TInd>
	void APC::updateData(const Population<TInd>& pop){
		m_N = pop.size();
		if (m_responsibility.size() != m_N) {
			m_responsibility.resize(m_N);
		}
		for (auto& row : m_responsibility)
			row.assign(m_N, 0);
		if (m_availability.size() != m_N) {
			m_availability.resize(m_N);
		}
		for (auto& row : m_availability)
			row.assign(m_N, 0);
		if (m_similarity.size() != m_N) {
			m_similarity.resize(m_N);
			for (auto& row : m_similarity)
				row.resize(m_N);
		}
		if (m_data.size() != m_N)
			m_data.resize(m_N);
		for (size_t i = 0; i < m_N; ++i)
			m_data[i] = &pop[i];
		m_clusters.clear();
		m_is_exemplar.assign(m_N, false);
	}

	template<typename TInd>
	void APC::updateData(const std::vector<TInd>& inds) {
		m_N = inds.size();
		if (m_responsibility.size() != m_N) {
			m_responsibility.resize(m_N);
			for (auto& row : m_responsibility)
				row.assign(m_N, 0);
		}
		if (m_availability.size() != m_N) {
			m_availability.resize(m_N);
			for (auto& row : m_availability)
				row.assign(m_N, 0);
		}
		if (m_similarity.size() != m_N) {
			m_similarity.resize(m_N);
			for (auto& row : m_similarity)
				row.resize(m_N);
		}
		if (m_data.size() != m_N)
			m_data.resize(m_N);
		for (size_t i = 0; i < m_N; ++i)
			m_data[i] = &inds[i];
		m_clusters.clear();
		m_is_exemplar.assign(m_N, false);
	}
}

#endif // !OFEC_APC_H

