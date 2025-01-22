#include "fcm_twl_v2.h"

namespace ofec {
	fcm_twl_v2::fcm_twl_v2(Real fuzziness, int num_cluster, Real error, const std::vector<size_t>& Data_points) :m_fuzziness(fuzziness), m_num_clusters(num_cluster), m_error(error), m_data_points(Data_points)
	{
		m_num_data_points = m_data_points.size();
		m_dim = 2;
		m_cluster_centre_index.resize(num_cluster);
		m_membership_tw.resize(m_num_data_points);
		for (auto &i : m_membership_tw)
			i.resize(m_num_clusters);
		m_membership_d.resize(m_num_data_points);
		for (auto &i : m_membership_d)
			i.resize(m_num_clusters);
		m_degree_of_memb.resize(m_num_data_points);
		for (auto &i : m_degree_of_memb)
			i.resize(m_num_clusters);

	}

	void fcm_twl_v2::clustering(fcm_twl_v2::Clustering_result &result)
	{
		Real max_diff = 0.0, min_diff=std::numeric_limits<double>::max();
		initialize_membership();
		
		do {
			m_d_min_max.first = std::numeric_limits<double>::max();
			m_d_min_max.second = 0;
			m_tw_min_max.first = std::numeric_limits<double>::max();
			m_tw_min_max.second = 0;
			get_tw_conflict_min_max();
			get_d_min_max();
			cacul_cluster_centre();
			/*max_diff = update_degree_of_membership();*/
			update_degree_of_membership();
			get_member(m_result);
			max_diff = 0.0;
			min_diff = std::numeric_limits<double>::max();
			for (int i = 0; i < m_result.member.size(); ++i) {
				if (max_diff < m_result.member[i].size())
					max_diff = m_result.member[i].size();
				if (min_diff > m_result.member[i].size())
					min_diff = m_result.member[i].size();
			}
			max_diff = max_diff - min_diff;
		} while (max_diff > 10);

		get_member(result);
	}

	void fcm_twl_v2::initialize_membership()
	{
		for (int i = 0; i < m_num_clusters; ++i) {
			size_t index = 0;
			do {
				index = ofec::global::ms_global->m_uniform[ofec::caller::Problem]->next_non_standard<size_t>(0, m_data_points.size());
			} while (std::find(m_cluster_centre_index.begin(), m_cluster_centre_index.end(), m_data_points[index]) != m_cluster_centre_index.end());
			m_cluster_centre_index[i] = m_data_points[index];
		}

		Real sum = 0.0, s = 0.0;
		for (int i = 0; i < m_num_data_points; ++i) {
			sum = 0.0;
			s = 0.0;
			for (int j = 0; j < m_num_clusters; ++j) {
				if (m_data_points[i] == m_cluster_centre_index[j])
					s = 0;
				else
					s = pow((1 / DVRP_CAST->get_shortest_dis(DVRP_CAST->get_depart_time(), m_data_points[i], m_cluster_centre_index[j])), 1 / (m_fuzziness - 1));
				sum += s;
			}
			for (int k = 0; k < m_num_clusters; ++k) {
				if (m_data_points[i] == m_cluster_centre_index[k])
					m_degree_of_memb[i][k] = 1;
				else
					m_degree_of_memb[i][k] = pow((1 / DVRP_CAST->get_shortest_dis(DVRP_CAST->get_depart_time(), m_data_points[i], m_cluster_centre_index[k])), 1 / (m_fuzziness - 1)) / sum;
			}
		}
		cacul_cluster_centre();
		get_member(m_result);
	}

	void fcm_twl_v2::get_member(Clustering_result & result)
	{
		result.membership = m_degree_of_memb;
		result.numClst = m_num_clusters;
		result.clusters = m_cluster_centre;
		result.clusters_index = m_cluster_centre_index;
		int cluster;
		Real highest;
		result.member.clear();
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

	void fcm_twl_v2::update_membership()
	{
		for (int i = 0; i < m_num_data_points; ++i) {
			for (int j = 0; j < m_num_clusters; ++j) {
				m_degree_of_memb[i][j] = 0.5*m_membership_d[i][j] + 0.5*m_membership_tw[i][j];
			}
		}
	}

	Real fcm_twl_v2::cacul_tw_conflict(size_t i, size_t j)
	{
		Real conflict = 0.0;
		Real diff_tw = 0.0;
		size_t cus1 = m_data_points[i];
		size_t cus2 = 0;
		for (int k = 0; k < m_result.member[j].size(); ++k) {
			cus2 = m_data_points[m_result.member[j][k]];
			if (cus1 == cus2) continue;
			Real cus1_beginTime = DVRP_CAST->get_net().get_road_net()[cus1]->datum_t.begin_time;
			Real cus1_endTime = DVRP_CAST->get_net().get_road_net()[cus1]->datum_t.end_time;
			Real cus2_beginTime = DVRP_CAST->get_net().get_road_net()[cus2]->datum_t.begin_time;
			Real cus2_endTime = DVRP_CAST->get_net().get_road_net()[cus2]->datum_t.end_time;
			Real beginTime = cus1_beginTime > cus2_beginTime ? cus1_beginTime : cus2_beginTime;
			Real endTime = cus1_endTime < cus2_endTime ? cus1_endTime : cus2_endTime;
			if (endTime >= beginTime) {
				diff_tw = endTime - beginTime;
				if (diff_tw <= 0)
					diff_tw = 0.1;
				diff_tw = diff_tw / (cus1_endTime - cus1_beginTime);
			}
			else
			{
				diff_tw = beginTime - endTime;
				diff_tw = 1 / diff_tw / (cus1_endTime - cus1_beginTime);
			}
			conflict += diff_tw;
		}
		return conflict;
	}

	Real fcm_twl_v2::cacul_d_conflict(size_t i, size_t j)
	{
		Real s = 0.0;
		if (m_data_points[i] == m_cluster_centre_index[j])
			s = 0;
		else
			s = DVRP_CAST->get_shortest_dis(DVRP_CAST->get_depart_time(), m_data_points[i], m_cluster_centre_index[j]);

		/*Real diff_x = DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.x - DVRP_CAST->get_net().get_road_net()[m_cluster_centre_index[j]]->datum_t.x;
		Real diff_y = DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.y - DVRP_CAST->get_net().get_road_net()[m_cluster_centre_index[j]]->datum_t.y;*/

		/*Real diff_x = DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.x -m_cluster_centre[j][0];
		Real diff_y = DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.y - m_cluster_centre[j][1];
		s = sqrt(pow(diff_x, 2) + pow(diff_y, 2));*/
		return s;
	}

	Real fcm_twl_v2::get_dis(size_t i, size_t j)
	{
		Real s = 0.0, s_d = 0.0, s_tw = 0.0;
		s_d = cacul_d_conflict(i, j) / (m_d_min_max.second - m_d_min_max.first);
		s_tw= cacul_tw_conflict(i, j) / (m_tw_min_max.second - m_tw_min_max.first);
		s = 0.8*s_d + 0.2*s_tw;
		//s = cacul_d_conflict(i, j);
		return s;
	}

	std::pair<Real, Real> fcm_twl_v2::get_tw_conflict_min_max()
	{
		/*std::pair<Real, Real> tw_min_max;
		tw_min_max.first = std::numeric_limits<double>::max();
		tw_min_max.second = 0;*/
		for (int i = 0; i < m_num_data_points; ++i) {
			for (int j = 0; j < m_num_clusters; ++j) {
				if (m_tw_min_max.first > cacul_tw_conflict(i, j))
					m_tw_min_max.first = cacul_tw_conflict(i, j);
				if (m_tw_min_max.second < cacul_tw_conflict(i, j))
					m_tw_min_max.second = cacul_tw_conflict(i, j);
			}
		}
		return m_tw_min_max;
	}

	std::pair<Real, Real> fcm_twl_v2::get_d_min_max()
	{
		/*std::pair<Real, Real> d_min_max;
		d_min_max.first = std::numeric_limits<double>::max();
		d_min_max.second = 0;*/
		for (int i = 0; i < m_num_data_points; ++i) {
			for (int j = 0; j < m_num_clusters; ++j) {
				if (m_d_min_max.first > cacul_d_conflict(i, j))
					m_d_min_max.first = cacul_d_conflict(i, j);
				if (m_d_min_max.second < cacul_d_conflict(i, j))
					m_d_min_max.second = cacul_d_conflict(i, j);
			}
		}
		return m_d_min_max;
	}

	void fcm_twl_v2::cacul_cluster_centre()
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
		for (auto &i : dis_centreTOpoints) {
			auto it = std::min_element(i.begin(), i.end());
			auto index = std::distance(i.begin(), it);
			m_cluster_centre_index.push_back(index);//����ͬ����
		}


		//std::vector<std::vector<Real>> dis_temp(m_result.member.size());
		//for (int i = 0; i < m_result.member.size(); ++i) {
		//	if (m_result.member[i].size() == 0) continue;
		//	for (int j = 0; j < m_result.member[i].size(); ++j) {
		//		Real sum = 0.0;
		//		size_t ii = m_result.member[i][j], jj = i;
		//		Real s = 0.0, s_d = 0.0, s_tw = 0.0;
		//		s_d = cacul_d_conflict(ii, jj) / (m_d_min_max.second - m_d_min_max.first);
		//		s_tw = cacul_tw_conflict(ii, jj) / (m_tw_min_max.second - m_tw_min_max.first);
		//		s = 0.6*s_d + 0.4*s_tw;
		//		sum += s;
		//		/*for (int k = 0; k < m_result.member[i].size(); ++k) {
		//			if (i == k) continue;
		//			size_t c1 = m_data_points[m_result.member[i][j]];
		//			size_t c2 = m_data_points[m_result.member[i][k]];
		//			sum += DVRP_CAST->get_shortest_dis()[c1][c2];
		//		}*/
		//		dis_temp[i].push_back(sum);
		//	}
		//}
		//for (int i = 0; i < dis_temp.size();++i) {
		//	auto it = std::min_element(dis_temp[i].begin(), dis_temp[i].end());
		//	auto index = std::distance(dis_temp[i].begin(), it);
		//	m_cluster_centre_index[i] = m_data_points[m_result.member[i][index]];
		//}

	}

	Real fcm_twl_v2::get_new_value(int i, int j)
	{
		int k;
		Real t, p, sum;
		sum = 0.0;
		p = 2 / (m_fuzziness - 1);
		auto temp = get_dis(i, j);
		for (k = 0; k < m_num_clusters; k++) {
			t = temp / get_dis(i, k);
			t = pow(t, p);
			sum += t;
		}
		return 1.0 / sum;
	}

	Real fcm_twl_v2::update_degree_of_membership()
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

	void fcm_twl_v2::normalize(std::vector<std::vector<Real>>& membership)
	{
		Real min = std::numeric_limits<double>::max();
		Real max = -1;
		for (int i = 0; i < m_num_data_points; ++i) {
			if (max <= *std::max_element(membership[i].begin(), membership[i].end()))
				max = *std::max_element(membership[i].begin(), membership[i].end());
			if (min >= *std::min_element(membership[i].begin(), membership[i].end()))
				min = *std::min_element(membership[i].begin(), membership[i].end());
		}
		for (int i = 0; i < m_num_data_points; ++i) {
			for (int j = 0; j < m_num_clusters; ++j) {
				membership[i][j] = (membership[i][j] - min) / (max - min);
			}
		}
	}


}