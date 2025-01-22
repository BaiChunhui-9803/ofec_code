/*
Agglomerative Hierarchical Clustering.
updated Jun 22, 2019 by Long Xiao
*/

#ifndef AHC_v2_H
#define AHC_v2_H


#include "../../../../../problem/realworld/DVRP/dynamic_vrp.h"
#include "../../../../../../core/global.h"

namespace ofec {
	class AHC_v2 {
	public:
		struct Clustering_result {
			std::vector<std::vector<size_t>> member;
		};
		AHC_v2(Real weight, const std::vector<size_t> &data_index) :m_weight(weight), m_data(data_index) {}
		void clustering(Clustering_result &result, Problem *pro);
	protected:
		std::vector<size_t> m_data;
		std::vector<std::vector<size_t>> m_cluster;
		std::vector<std::vector<Real>> m_dis_matrix;
		std::vector<std::vector<Real>> m_d_conflict_matrix;
		std::vector<std::vector<Real>> m_tw_conflict_matrix;
		std::vector<Real> m_objective;
		Real m_weight;
		Real m_ca_cons_co = 1.0;//Capacity constraint coefficient

		const std::vector<std::vector<Real>>& cacul_dis_matrix(const std::vector<size_t> &data, Problem *pro);
		Real cacul_d_conflict(size_t i, size_t j, Problem *pro);
		Real cacul_tw_conflict(size_t i, size_t j, Problem *pro);
		void cacul_objective(const std::vector<std::vector<size_t>> &cluster);
		//std::vector<size_t> tsp(std::vector<size_t> &member_index);
		//Real get_cost(size_t s, size_t e);
		const bool is_valid(const std::vector<size_t> &cluster, const std::vector<size_t> &data, Problem *pro) const;
	};
}
#endif // !AHC_v2_H