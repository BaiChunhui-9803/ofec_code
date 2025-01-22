/*
Agglomerative Hierarchical Clustering.
updated Jun 22, 2019 by Long Xiao
*/

#ifndef AHC_H
#define AHC_H

#include <iostream>
#include <vector>
#include "../../../../../problem/realworld/DVRP/dynamic_vrp.h"
#include "../../../../../../core/global.h"

namespace ofec {
	class AHC {
	public:
		struct Clustering_result {
			std::vector<std::vector<size_t>> member;
		};
		AHC(int num_cluster, const std::vector<size_t> &data_index) :m_num_cluster(num_cluster), m_data(data_index) {}
		void clustering(Clustering_result &result);
	protected:
		std::vector<size_t> m_data;
		std::vector<std::vector<size_t>> m_cluster;
		std::vector<std::vector<Real>> m_dis_matrix;
		std::vector<std::vector<Real>> m_d_conflict_matrix;
		std::vector<std::vector<Real>> m_tw_conflict_matrix;
		std::vector<Real> m_objective;
		Real m_weight = 0.6;
		int m_num_cluster;
		Real m_ca_cons_co = 0.75;//Capacity constraint coefficient

		std::vector<std::vector<Real>> cacul_dis_matrix(const std::vector<size_t> &data);
		Real cacul_d_conflict(size_t i, size_t j);
		Real cacul_tw_conflict(size_t i, size_t j);
		//void cacul_objective(std::vector<std::vector<size_t>> &cluster);
		//std::vector<size_t> tsp(std::vector<size_t> &member_index);
		//Real get_cost(size_t s, size_t e);
		bool is_valid(const std::vector<size_t> &cluster, const std::vector<size_t> &data);
	};
}
#endif // !AHC_H