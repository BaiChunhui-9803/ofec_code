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
* APD (Adaptive Peak Detection)
* Ran Cheng, Miqing Li, Ke Li, and Xin Yao
* "Evolutionary Multiobjective Optimization-Based Multimodal Optimization: Fitness Landscape Approximation and Peak Detection."
* IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL. 22, NO. 5, OCTOBER 2018
*-------------------------------------------------------------------------------
* Implemented by Junchen Wang (wangjunchen.chris@gmail.com) at 2020/10/11
*********************************************************************************/

#ifndef OFEC_APD_H
#define OFEC_APD_H

#include <algorithm>
#include <vector>
#include <list>
#include <deque>
#include <unordered_set>

#include "../../core/algorithm/population.h"
#include "../../core/algorithm/Solution.h"
#include "../../core/problem/continuous/continuous.h"

namespace ofec {
	template<typename Solution>
	class APD {
	protected:
		real m_eta;                                     // initial cutting ratio
		std::vector<const solution<>*> m_D;                   // filtered unique solutions 
		std::vector<std::vector<real>> m_manh_dis_mat;
		std::list<const solution<>*> m_peaks;
	public:
		APD(real eta = 0.1) : m_eta(eta) {}
		void updateData(const population<Solution> &pop);
		void updateData(const std::vector<Solution> &inds);
		void clustering();
		const std::list<const solution<>*>& get_peaks() const { return m_peaks; }
	protected:
		void updateDistance(const std::unordered_set<size_t> &D_c);
		void landscapeCutting(std::unordered_set<size_t> &D_c, real eta);
		void peakDetection(std::unordered_set<size_t> D_c);
		void normManhDisMat(const std::unordered_set<size_t> &D_c);
		real minManhDis(const std::unordered_set<size_t>& D_c, size_t i);
	};
	
	template<typename Solution>
	void APD<Solution>::updateData(const population<Solution> &pop) {
		m_D.clear();
		std::set<std::vector<real>> D_x;    // for uniqueness
		for (size_t i = 0; i < pop.size(); i++)	{
			if (D_x.count(pop[i].variable().vect()) == 0) {
				m_D.push_back(&pop[i]);
				D_x.insert(pop[i].variable().vect());
			}
		}
		size_t N = m_D.size();
		if (m_manh_dis_mat.size() != N)
			m_manh_dis_mat.assign(N, std::vector<real>(N));
	}
	
	template<typename Solution>
	void APD<Solution>::updateData(const std::vector<Solution> &inds) {
		m_D.clear();
		std::set<std::vector<real>> D_x;    // for uniqueness
		for (size_t i = 0; i < inds.size(); i++) {
			if (D_x.count(inds[i].variable().vect()) == 0) {
				m_D.push_back(&inds[i]);
				D_x.insert(inds[i].variable().vect());
			}
		}
		size_t N = m_D.size();
		if (m_manh_dis_mat.size() != N)
			m_manh_dis_mat.assign(N, std::vector<real>(N));
	}

	template<typename Solution>
	void APD<Solution>::clustering() {
		std::unordered_set<size_t> D_c;
		for (size_t i = 0; i < m_D.size(); i++)
			D_c.insert(i);
		landscapeCutting(D_c, m_eta);
		updateDistance(D_c);
		while (!D_c.empty()) {
			peakDetection(D_c);
			landscapeCutting(D_c, 0.5);
		}
	}

	template<typename Solution>
	void APD<Solution>::updateDistance(const std::unordered_set<size_t>& D_c) {
		for (auto iter1 = D_c.begin(); iter1 != D_c.end(); iter1++) {
			auto iter2 = iter1;
			iter2++;
			for (; iter2 != D_c.end(); iter2++) {
				m_manh_dis_mat[*iter1][*iter2] = manhattan_distance(m_D[*iter1]->variable().begin(), m_D[*iter1]->variable().end(), m_D[*iter2]->variable().begin());
				m_manh_dis_mat[*iter2][*iter1] = m_manh_dis_mat[*iter1][*iter2];
			}
		}
	}

	template<typename Solution>
	void APD<Solution>::landscapeCutting(std::unordered_set<size_t> &D_c, real eta) {
		real min_obj, max_obj;
		auto iter = D_c.begin();
		min_obj = max_obj = m_D[*iter]->objective(0);
		while (++iter != D_c.end()) {
			if (m_D[*iter]->objective(0) > max_obj)
				max_obj = m_D[*iter]->objective(0);
			if (m_D[*iter]->objective(0) < min_obj)
				min_obj = m_D[*iter]->objective(0);
		}
		real threshold;
		if (CONTINUOUS_CAST->opt_mode(0) == optimization_mode::Maximization) {
			threshold = max_obj - eta * (max_obj - min_obj);
			for (auto iter = D_c.begin(); iter != D_c.end();) {
				if (m_D[*iter]->objective(0) > threshold)
					++iter;
				else
					iter = D_c.erase(iter);
			}
		}
		else {
			threshold = min_obj + eta * (max_obj - min_obj);
			for (auto iter = D_c.begin(); iter != D_c.end();) {
				if (m_D[*iter]->objective(0) < threshold)
					++iter;
				else
					iter = D_c.erase(iter);
			}
		}
		if (D_c.size() > 5000) {
			std::vector<size_t> D_c_cp;
			for (size_t idx_ind : D_c)
				D_c_cp.push_back(idx_ind);
			global::ms_global->m_uniform[caller::Algorithm]->shuffle(D_c_cp.begin(), D_c_cp.end());
			D_c.clear();
			for (size_t i = 0; i < 5000; ++i)
				D_c.insert(D_c_cp[i]);
		}
	}

	template<typename Solution>
	void APD<Solution>::peakDetection(std::unordered_set<size_t> D_c) {
		size_t num_B = D_c.size();
		normManhDisMat(D_c);
		while (num_B > 1) {
			auto iter = D_c.begin();
			auto iter_sigma = iter;
			real min_dis, sigma = minManhDis(D_c, *iter);
			while (++iter != D_c.end()) {
				min_dis = minManhDis(D_c, *iter);
				if (min_dis > sigma) {
					sigma = min_dis;
					iter_sigma = iter;
				}
			}
			std::list<size_t> psi_k;
			std::deque<size_t> neighbors = { *iter_sigma };
			D_c.erase(iter_sigma);
			while (!neighbors.empty()) {
				psi_k.push_back(neighbors.front());
				for (iter = D_c.begin(); iter != D_c.end();) {
					if (m_manh_dis_mat[neighbors.front()][*iter] <= std::max<real>(sigma, 1e-3)) {
						neighbors.push_back(*iter);
						iter = D_c.erase(iter);
					}
					else
						iter++;
				}
				num_B--;
				neighbors.pop_front();
			}
			auto iter_best = psi_k.begin(), iter_psi = iter_best;
			while (++iter_psi != psi_k.end()) {
				if (m_D[*iter_psi]->dominate(*m_D[*iter_best]))
					iter_best = iter_psi;
			}
			bool flag = true;
			for (auto ptr_sol : m_peaks) {
				if (ptr_sol == m_D[*iter_best]) {
					flag = false;
					break;
				}
			}
			if (flag)
				m_peaks.push_back(m_D[*iter_best]);
		}
	}
	
	template<typename Solution>
	void APD<Solution>::normManhDisMat(const std::unordered_set<size_t>& D_c) {
		std::vector<std::vector<real>> normc_x(D_c.size(), std::vector<real>(CONTINUOUS_CAST->variable_size()));
		for (size_t j = 0; j < CONTINUOUS_CAST->variable_size(); j++) {
			real sum = 0;
			for (size_t idx_ind : D_c) {
				sum += pow(m_D[idx_ind]->variable()[j], 2);
			}
			real beta = sqrt(1 / sum);
			size_t i = 0;
			for (size_t idx_ind : D_c) {
				normc_x[i][j] = beta * m_D[idx_ind]->variable()[j];
				i++;
			}
		}
		size_t i = 0;
		for (auto iter1 = D_c.begin(); iter1 != D_c.end(); iter1++) {
			size_t k = i + 1;
			auto iter2 = iter1;
			iter2++;
			for (; iter2 != D_c.end(); iter2++) {
				m_manh_dis_mat[*iter1][*iter2] = manhattan_distance(normc_x[i].begin(), normc_x[i].end(), normc_x[k].begin());
				m_manh_dis_mat[*iter2][*iter1] = m_manh_dis_mat[*iter1][*iter2];
				k++;
			}
			i++;
		}
	}

	template<typename Solution>
	real APD<Solution>::minManhDis(const std::unordered_set<size_t>& D_c, size_t i) {
		auto iter = D_c.begin();
		if (*iter == i) iter++;
		real min_dis = m_manh_dis_mat[i][*iter];
		while (++iter != D_c.end()) {
			if (*iter == i) continue;
			if (m_manh_dis_mat[i][*iter] < min_dis)
				min_dis = m_manh_dis_mat[i][*iter];
		}
		return min_dis;
	}
}

#endif // !OFEC_APD_H
