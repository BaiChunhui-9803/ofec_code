#define _CRT_SECURE_NO_DEPRECATE
#include "fcm_m.h"

namespace ofec {
		void fcm_m::init(const std::vector<size_t> &Data_points)
		{
			m_data_points = Data_points;
			m_num_data_points = m_data_points.size();


			//m_min_max.resize(m_dim);
			//for (int i = 0; i < m_dim; ++i) {
			//	m_min_max[i].first = DBL_MAX;
			//}
			//for (int i = 0; i < DVRP_CAST->get_net().get_road_net().size(); i++) {
			//	
			//	if (DVRP_CAST->get_net().get_road_net()[i]->datum_t.x < m_min_max[0].first)
			//		m_min_max[0].first = DVRP_CAST->get_net().get_road_net()[i]->datum_t.x;
			//	if (DVRP_CAST->get_net().get_road_net()[i]->datum_t.x > m_min_max[0].second)
			//		m_min_max[0].second = DVRP_CAST->get_net().get_road_net()[i]->datum_t.x;
			//	if (DVRP_CAST->get_net().get_road_net()[i]->datum_t.y < m_min_max[1].first)
			//		m_min_max[1].first = DVRP_CAST->get_net().get_road_net()[i]->datum_t.y;
			//	if (DVRP_CAST->get_net().get_road_net()[i]->datum_t.y > m_min_max[1].second)
			//		m_min_max[1].second = DVRP_CAST->get_net().get_road_net()[i]->datum_t.y;
			//	
			//}
			////Normalized
			//m_normal_pos.resize(DVRP_CAST->get_net().get_road_net().size());
			//for (int i = 0; i < DVRP_CAST->get_net().get_road_net().size(); ++i) {
			//	m_normal_pos[i].push_back((DVRP_CAST->get_net().get_road_net()[i]->datum_t.x - m_min_max[0].first) / (m_min_max[0].second-m_min_max[0].first));
			//	m_normal_pos[i].push_back((DVRP_CAST->get_net().get_road_net()[i]->datum_t.y - m_min_max[1].first) / (m_min_max[1].second-m_min_max[1].first));
			//}

			initialize_membership();
			//initial membership
			//m_degree_of_memb.resize(m_num_data_points);
			//Real s;
			//int r, rval;
			//for (int i = 0; i < m_num_data_points; ++i) {
			//	s = 0.0;
			//	r = 100;
			//	m_degree_of_memb[i].resize(m_num_clusters);
			//	for (int j = 1; j < m_num_clusters; ++j) {
			//		rval = ofec::global::ms_global->m_uniform[ofec::caller::Problem]->next_non_standard<size_t>(0, r + 1);
			//		//rval = rand() % (r + 1);
			//		r -= rval;
			//		m_degree_of_memb[i][j] = rval / 100.0;
			//		s += m_degree_of_memb[i][j];
			//	}
			//	m_degree_of_memb[i][0] = 1.0 - s;
			//}
		}

		void fcm_m::initialize_membership()
		{
			const Real pi = 3.14159265;
			Real s = 0.0;
			Real r = 0.02;
			Real sum = 0.0;
			Real angle = 360 / m_num_clusters;
			std::pair<Real, Real> depot;
			depot.first = DVRP_CAST->get_net().get_road_net()[DVRP_CAST->get_depot()]->datum_t.x;
			depot.second = DVRP_CAST->get_net().get_road_net()[DVRP_CAST->get_depot()]->datum_t.y;
			/*depot.first = m_normal_pos[DVRP_CAST->get_depot()][0];
			depot.second = m_normal_pos[DVRP_CAST->get_depot()][1];*/
			m_cluster_centre.resize(m_num_clusters);
			for (int i = 0; i < m_num_clusters; ++i) {
				m_cluster_centre[i].push_back(depot.first + r * cos((angle*i)*pi/180));
				m_cluster_centre[i].push_back(depot.second + r * sin((angle*i)*pi / 180));
			}

			std::vector<std::vector<Real>> dis_centreTOpoints(m_num_clusters);
			for (int i = 0; i < m_cluster_centre.size(); ++i) {
				sum = 0.0;
				for (int j = 0; j < DVRP_CAST->get_net().get_road_net().size(); ++j) {
					sum = pow(m_cluster_centre[i][0] - DVRP_CAST->get_net().get_road_net()[j]->datum_t.x, 2) + pow(m_cluster_centre[i][1] - DVRP_CAST->get_net().get_road_net()[j]->datum_t.y, 2);
					dis_centreTOpoints[i].push_back(sqrt(sum));
				}
			}
			m_cluster_centre_index.clear();
			for (auto &i : dis_centreTOpoints) {
				auto it = std::min_element(i.begin(), i.end());
				auto index = std::distance(i.begin(), it);
				m_cluster_centre_index.push_back(index);
			}

			m_degree_of_memb.resize(m_num_data_points);
			for (int i = 0; i < m_num_data_points; ++i) {
				m_degree_of_memb[i].resize(m_num_clusters);
				sum = 0.0;
				s = 0.0;
				for (int j = 0; j < m_num_clusters; ++j) {
					if (m_data_points[i] == m_cluster_centre_index[j])
						s = std::numeric_limits<double>::min();
					else
						s= pow((1/DVRP_CAST->get_shortest_dis(DVRP_CAST->get_depart_time(), m_data_points[i], m_cluster_centre_index[j])),1/(m_fuzziness-1));
					sum += s;
				}
				for (int k = 0; k < m_num_clusters; ++k) {
					if (m_data_points[i] == m_cluster_centre_index[k])
						m_degree_of_memb[i][k] = std::numeric_limits<double>::min();
					else
						m_degree_of_memb[i][k] = pow((1 / DVRP_CAST->get_shortest_dis(DVRP_CAST->get_depart_time(), m_data_points[i], m_cluster_centre_index[k])), 1 / (m_fuzziness - 1)) / sum;
					//s += m_degree_of_memb[i][k];
				}
				//m_degree_of_memb[i][0] = 1.0 - s;
			}

			for (int i = 0; i < m_num_data_points; ++i) {
				auto max_it = std::max_element(m_degree_of_memb[i].begin(), m_degree_of_memb[i].end());
				auto max_pos = std::distance(m_degree_of_memb[i].begin(), max_it);
				s = 0.0;
				for (int j = 0; j < m_num_clusters; ++j) {
					if(j!=max_pos)
						s += m_degree_of_memb[i][j];
				}
				m_degree_of_memb[i][max_pos] = 1 - s;
			}

		}

		void fcm_m::calculate_centre_vectors()
		{
			Real numerator, denominator;
			std::vector<std::vector<Real>> t;
			t.resize(m_num_data_points);
			for (int i = 0; i < m_num_data_points; ++i) {
				t[i].resize(m_num_clusters);
				for (int j = 0; j < m_num_clusters; ++j) {
					t[i][j] = pow(m_degree_of_memb[i][j], m_fuzziness);
				}
			}
			m_cluster_centre.resize(m_num_clusters);
			for (int j = 0; j < m_num_clusters; j++) {
				m_cluster_centre[j].resize(m_dim);
				numerator = 0.0;
				denominator = 0.0;
				for (int i = 0; i < m_num_data_points; i++) {
					numerator += t[i][j] * DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.x;//��һ������,���ԣ���Ч��
					denominator += t[i][j];
				}
				m_cluster_centre[j][0] = numerator / denominator;
				numerator = 0.0;
				denominator = 0.0;
				for (int i = 0; i < m_num_data_points; i++) {
					numerator += t[i][j] * DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.y;
					denominator += t[i][j];
				}
				m_cluster_centre[j][1] = numerator / denominator;
			}

			std::vector<std::vector<Real>> dis_centreTOpoints(m_num_clusters);
			for (int i = 0; i < m_cluster_centre.size(); ++i) {
				for (int j = 0; j < DVRP_CAST->get_net().get_road_net().size(); ++j) {
					Real sum = pow(m_cluster_centre[i][0] - DVRP_CAST->get_net().get_road_net()[j]->datum_t.x, 2) + pow(m_cluster_centre[i][1] - DVRP_CAST->get_net().get_road_net()[j]->datum_t.y, 2);
					dis_centreTOpoints[i].push_back(sqrt(sum));
				}
			}
			m_cluster_centre_index.clear();
			for (auto &i:dis_centreTOpoints) {
				auto it = std::min_element(i.begin(), i.end());
				auto index = std::distance(i.begin(), it);
				m_cluster_centre_index.push_back(index);//����ͬ����
			}
			/*for (int i = 0; i < m_cluster_centre_index.size(); ++i) {
				int pos = 0;
				auto it = std::find(m_cluster_centre_index.begin(), m_cluster_centre_index.end(), m_cluster_centre_index[i]);
				pos = std::distance(m_cluster_centre_index.begin(), it);
				if(pos==i){
					continue;
				}
				else {
					auto temp = m_cluster_centre_index[i] + ofec::global::ms_global->m_uniform[ofec::caller::Problem]->next_non_standard<size_t>(1, 3);
					if (temp > DVRP_CAST->get_net().get_road_net().size() - 1)
						temp = DVRP_CAST->get_net().get_road_net().size() - 1;
					m_cluster_centre_index[pos] = temp;
				}

			}*/
		}

		void fcm_m::clustering(const std::vector<size_t> &Data_points, fcm_m::Clustering_result &result)
		{
			//bool flag = true;
			Real max_diff;
			init(Data_points);
			do {
				calculate_centre_vectors();
				max_diff = update_degree_of_membership();
			} while (max_diff > m_epsilon);//max_diff > m_epsilon

			result.membership = m_degree_of_memb;
			result.numClst = m_num_clusters;
			result.clusters = m_cluster_centre;
			result.clusters_index = m_cluster_centre_index;
			int cluster;
			Real highest;
			result.member.resize(m_num_clusters);
			for (int i = 0; i < m_num_data_points; i++) {
				cluster = 0;
				highest = 0.0;
				for (int j = 0; j < m_num_clusters; j++) {
					if (m_degree_of_memb[i][j] > highest) {
						highest = m_degree_of_memb[i][j];
						cluster = j;
					}
				}
				result.member[cluster].push_back(i);
			}

			
		}

		int fcm_m::gnuplot_membership_matrix()
		{
			int i, j, cluster;
			char fname[100];
			Real highest;
			FILE * f[100];
			if (m_dim != 2) {
				printf("Plotting the cluster only works when the\n");
				printf("number of dimensions is two. This will create\n");
				printf("a two-dimensional plot of the cluster points.\n");
				exit(1);
			}
			for (j = 0; j < m_num_clusters; j++) {
				sprintf(fname, "cluster.%d", j);
				if ((f[j] = fopen(fname, "w")) == NULL) {
					printf("Could not create %s\n", fname);
					for (i = 0; i < j; i++) {
						fclose(f[i]);
						sprintf(fname, "cluster.%d", i);
						remove(fname);
					}
					return -1;
				}
				fprintf(f[j], "#Data points for cluster: %d\n", j);
			}
			for (i = 0; i < m_num_data_points; i++) {
				cluster = 0;
				highest = 0.0;
				for (j = 0; j < m_num_clusters; j++) {
					if (m_degree_of_memb[i][j] > highest) {
						highest = m_degree_of_memb[i][j];
						cluster = j;
					}
				}
				fprintf(f[cluster], "%lf %lf\n", DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.x, DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.y);
			}
			for (j = 0; j < m_num_clusters; j++) {
				fclose(f[j]);
			}
			if ((f[0] = fopen("gnuplot.script", "w")) == NULL) {
				printf("Could not create gnuplot.script.\n");
				for (i = 0; i < j; i++) {
					fclose(f[i]);
					sprintf(fname, "cluster.%d", i);
					remove(fname);
				}
				return -1;
			}
			fprintf(f[0], "set terminal png medium\n");
			fprintf(f[0], "set output \"cluster_plot.png\"\n");
			fprintf(f[0], "set title \"FCM clustering\"\n");
			fprintf(f[0], "set xlabel \"x-coordinate\"\n");
			fprintf(f[0], "set ylabel \"y-coordinate\"\n");
			fprintf(f[0], "set xrange [%lf : %lf]\n", m_min_max[0].first, m_min_max[0].second);
			fprintf(f[0], "set yrange [%lf : %lf]\n", m_min_max[1].first, m_min_max[1].second);
			fprintf(f[0],
				"plot 'cluster.0' using 1:2 with points pt 7 ps 1 lc 1 notitle");
			for (j = 1; j < m_num_clusters; j++) {
				sprintf(fname, "cluster.%d", j);
				fprintf(f[0],
					",\\\n'%s' using 1:2 with points  pt 7 ps 1 lc %d notitle",
					fname, j + 1);
			}
			fprintf(f[0], "\n");
			fclose(f[0]);
			return 0;
		}

		void fcm_m::print_membership_matrix(const char *fname)
		{
			int i, j;
			FILE *f;
			if (fname == NULL)
				f = stdout;
			else if ((f = fopen(fname, "w")) == NULL) {
				printf("Cannot create output file.\n");
				exit(1);
			}
			fprintf(f, "Membership matrix:\n");
			for (i = 0; i < m_num_data_points; i++) {
				fprintf(f, "Data[%d]: ", i);
				for (j = 0; j < m_num_clusters; j++) {
					fprintf(f, "%lf ", m_degree_of_memb[i][j]);
				}
				fprintf(f, "\n");
			}
			if (fname == NULL)
				fclose(f);
		}

		Real fcm_m::get_norm(int i, int j)
		{
			Real sum = 0.0;
			/*sum = pow(DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.x - m_cluster_centre[j][0], 2)+ pow(DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.y - m_cluster_centre[j][1], 2);
			return sqrt(sum);*/
			//std::vector<size_t> path2centre;
			if (m_data_points[i] == m_cluster_centre_index[j]) { 
				sum = std::numeric_limits<double>::min();
			}
			else {
				/*DVRP_CAST->shortest_path(m_data_points[i], m_cluster_centre_index[j], path2centre, DVRP_CAST->get_depart_time());
				for (int k = 0; k < path2centre.size() - 1; ++k) {
					sum += DVRP_CAST->get_net().get_node_dis()[path2centre[k]][path2centre[k + 1]];
				}*/
				sum = DVRP_CAST->get_shortest_dis(DVRP_CAST->get_depart_time(), m_data_points[i], m_cluster_centre_index[j]);
			}
			return sum;
			
		}

		Real fcm_m::get_new_value(int i, int j)
		{
			int k;
			Real t, p, sum;
			sum = 0.0;
			p = 2 / (m_fuzziness - 1);
			for (k = 0; k < m_num_clusters; k++) {
				t = get_norm(i, j) / get_norm(i, k);
				t = pow(t, p);
				sum += t;
			}
			return 1.0 / sum;
		}

		Real fcm_m::update_degree_of_membership()
		{
			int i, j;
			Real new_uij;
			Real max_diff = 0.0, diff;
			for (j = 0; j < m_num_clusters; j++) {
				for (i = 0; i < m_num_data_points; i++) {
					new_uij = get_new_value(i, j);
					diff = new_uij - m_degree_of_memb[i][j];
					if (diff > max_diff)
						max_diff = diff;
					m_degree_of_memb[i][j] = new_uij;
				}
			}
			return max_diff;
		}
}