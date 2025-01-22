#ifndef FCM_TWL_V2_H
#define FCM_TWL_V2_H


#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <locale>
#include <tuple>
#include "../../../../../problem/realworld/DVRP/dynamic_vrp.h"
#include "../../../../../../core/global.h"
//#include <windows.h>

namespace ofec {
	class fcm_twl_v2
	{
	public:
		struct Clustering_result {
			std::vector<std::vector<Real>> membership;
			std::vector<std::vector<Real>> clusters;
			std::vector<size_t> clusters_index;
			int numClst;
			std::vector<std::vector<size_t>> member;
		};
		fcm_twl_v2(Real fuzziness, int num_cluster, Real error, const std::vector<size_t> & Data_points);
		void clustering(fcm_twl_v2::Clustering_result &);



	private:
		size_t m_dim;
		//std::vector<size_t> m_data_points;
		std::vector<size_t> m_data_points;
		size_t m_num_data_points = 0;
		size_t m_num_clusters;
		Real m_fuzziness;
		Real m_error = 0.0;
		std::pair<Real, Real> m_d_min_max;
		std::pair<Real, Real> m_tw_min_max;
		std::vector<std::vector<Real>> m_degree_of_memb;
		std::vector<std::vector<Real>> m_membership_d;
		std::vector<std::vector<Real>> m_membership_tw;
		std::vector<std::vector<Real>> m_cluster_centre;
		std::vector<size_t> m_cluster_centre_index;
		Clustering_result m_result;

		void initialize_membership();
		void get_member(Clustering_result &result);
		void update_membership();
		Real cacul_tw_conflict(size_t i, size_t j);
		Real cacul_d_conflict(size_t i, size_t j);
		Real get_dis(size_t i, size_t j);
		std::pair<Real, Real> get_tw_conflict_min_max();
		std::pair<Real, Real> get_d_min_max();
		void cacul_cluster_centre();
		Real get_new_value(int i, int j);
		Real update_degree_of_membership();
		void normalize(std::vector<std::vector<Real>> &);

	};
}



#endif // !FCM_TWL_V2_H

