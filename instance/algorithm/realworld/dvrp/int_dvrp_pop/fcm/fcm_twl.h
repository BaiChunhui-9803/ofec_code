#ifndef FCM_TWL_H
#define FCM_TWL_H

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
	class fcm_twl
	{
	public:
		struct Clustering_result {
			std::vector<std::vector<Real>> membership;
			std::vector<std::vector<Real>> clusters;
			std::vector<size_t> clusters_index;
			int numClst;
			std::vector<std::vector<size_t>> member;
		};
		fcm_twl(int num_cluster, const std::vector<size_t> & Data_points);
		void get_result(fcm_twl::Clustering_result &);



	private:
		size_t m_dim;
		//std::vector<size_t> m_data_points;
		std::vector<size_t> m_data_points;
		size_t m_num_data_points = 0;
		size_t m_num_clusters;
		std::vector<std::vector<Real>> m_degree_of_memb;
		std::vector<std::vector<Real>> m_cluster_centre;
		std::vector<size_t> m_cluster_centre_index;
		//Clustering_result m_result;

		void initialize_membership(fcm_twl::Clustering_result &);
		void get_member(Clustering_result &result);
		std::vector<size_t> tsp(std::vector<size_t>&);
		void get_cost(size_t, size_t, Real &);

	};
}
#endif // !FCM_TWL_H

