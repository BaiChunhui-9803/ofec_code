/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Long Xiao
* Email: 917003976@qq.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://gitlab.com/ofec/console for more information
*
*********************************************************************************/
/*
Updated Feb 8, 2019 by Long Xiao
Agglomerative Hierarchical Clustering.
updated Jun 22, 2019 by Long Xiao
*/

#ifndef OFEC_AHC_TSP_H
#define OFEC_AHC_TSP_H

#include "../../../../../problem/realworld/DVRP/dynamic_vrp.h"
#include "../../../../../../core/global.h"
#include "../../lkh/include/lkh.h"

#if defined(WIN32)
#include <windows.h>
#endif

namespace ofec {
	class AHC_TSP {
	public:
		struct Clustering_result {
			std::vector<std::vector<size_t>> member;
		};
		AHC_TSP(Real co_d, const std::vector<size_t> &data_index) :m_co_d(co_d), m_data(data_index) {}
		void clustering(Clustering_result &result, Problem *pro);
	protected:
		LKH::LKHAlg lkh;
		std::vector<size_t> m_data;
		std::vector<std::vector<size_t>> m_cluster;
		std::vector<std::vector<Real>> m_dis_matrix;
		std::vector<std::vector<Real>> m_d_conflict_matrix;
		std::vector<std::vector<Real>> m_tw_conflict_matrix;
		std::vector<Real> m_objective;
		//Real m_weight;
		Real m_ca_cons_co = 1.0;//Capacity constraint coefficient
		Real m_co_d = 0.3;
		bool m_lkh = false;

		//Real cacul_dis(std::vector<size_t>&member);
		void cacul_dis_matrix(std::vector<std::vector<Real>> &dis_temp, 
			const std::vector<std::vector<size_t>> &cluster_temp, const std::vector<size_t> &data, Problem *pro);
		std::pair<Real, Real> tsp_nearest(std::vector<size_t> &member_index, Problem *pro);
		std::pair<Real, Real> tsp_lkh(std::vector<size_t> &member_index, Problem *pro);
		//void tsp_lkh(std::vector<size_t> &member_index, const std::vector<std::vector<Real>> &dis_matrix);
		void cal_tw_conflict(const std::vector<size_t> &member_index, std::vector<std::vector<Real>> &dis_matrix, Problem *pro);
		void get_cost(size_t s, size_t e, Real &present_time, Problem *pro);
		const bool is_valid(const std::vector<size_t> &cluster, const std::vector<size_t> &data, Problem *pro) const;
	};
}
#endif // !OFEC_AHC_TSP_H