#ifndef FCM_M_H
#define FCM_M_H
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <locale>
#include "../../../../../problem/realworld/DVRP/dynamic_vrp.h"
#include "../../../../../../core/global.h"


namespace ofec {
		class fcm_m
		{
		public:
			struct Clustering_result {
				std::vector<std::vector<Real>> membership;
				std::vector<std::vector<Real>> clusters;
				std::vector<size_t> clusters_index;
				int numClst;
				std::vector<std::vector<size_t>> member;
			};
			fcm_m() = default;
			fcm_m(const size_t Dim, const Real Fuzziness, const Real Epsilon, const size_t Num_clusters) :m_dim(Dim), m_fuzziness(Fuzziness), m_epsilon(Epsilon), m_num_clusters(Num_clusters) {
				if (m_fuzziness <= 1.0)
				{
					std::cout << "The fuzziness must be > 1.0 !" << std::endl;
					exit(1);
				}
				if (m_epsilon > 1.0 || m_epsilon <= 0.0)
				{
					std::cout << "The epsilon must be <= 1.0 !" << std::endl;
					exit(1);
				}
				/*if (ofec::global::ms_global->m_problem.get() != nullptr && ofec::global::ms_global->m_problem->has_tag(ofec::problem_tag::DVRP)) {
					m_dvrp = dynamic_cast<ofec::dynamic_vrp*>(ofec::global::ms_global->m_problem.get());
				}*/
			}
			void set_parameter(const size_t Dim, const Real Fuzziness, const Real Epsilon, const size_t Num_clusters) {
				m_dim = Dim;
				m_fuzziness = Fuzziness;
				m_epsilon = Epsilon; 
				m_num_clusters = Num_clusters;
			}

			void clustering(const std::vector<size_t> &, fcm_m::Clustering_result &);
			int gnuplot_membership_matrix();
			void print_membership_matrix(const char *);


		private:
			//ofec::dynamic_vrp* m_dvrp = nullptr;
			size_t m_dim;
			Real m_fuzziness;
			Real m_epsilon;
			std::vector<std::pair<Real, Real>> m_min_max;
			std::vector<size_t> m_data_points;
			std::vector<std::vector<Real>> m_normal_pos;
			std::vector<std::vector<Real>> m_matrix_dis;
			size_t m_num_data_points = 0;
			size_t m_num_clusters;
			std::vector<std::vector<Real>> m_degree_of_memb;
			std::vector<std::vector<Real>> m_cluster_centre;
			std::vector<size_t> m_cluster_centre_index;

			void init(const std::vector<size_t> &Data_points);
			void initialize_membership();
			void calculate_centre_vectors();
			Real get_norm(int i, int j);
			Real get_new_value(int i, int j);
			Real update_degree_of_membership();


		};
}
#endif
