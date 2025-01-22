#include "AHC_v2.h"

//#include "../../../../../problem/realworld/DVRP/dynamic_vrp.h"

namespace ofec {
	void AHC_v2::clustering(Clustering_result &result, Problem *pro) {
		cacul_dis_matrix(m_data, pro);
		auto dis_matrix = m_dis_matrix;
		auto dis_temp = dis_matrix;
		//std::vector<size_t> member_index;
		std::vector<std::vector<size_t>> cluster_temp(m_data.size());
		for (int i = 0; i < m_data.size(); ++i) {
			cluster_temp[i].push_back(i);
		}
		auto data = m_data;
		auto data_temp = m_data;

		while (cluster_temp.size() != 1) {
			double minDt = std::numeric_limits<double>::max();
			size_t mi = 0, mj = 0;
			bool flag_erase = false;
			for (int i = 0; i < dis_temp.size(); ++i) {
				for (int j = 0; j < dis_temp[i].size(); ++j) {
					if (i == j) continue;
					if (minDt > dis_temp[i][j]) {
						minDt = dis_temp[i][j];
						mi = i;
						mj = j;
					}
				}
			}

			for (int i = 0; i < cluster_temp[mj].size(); ++i) {
				cluster_temp[mi].push_back(cluster_temp[mj][i]);
				if (is_valid(cluster_temp[mi], data, pro)) {
					flag_erase = true;
					cluster_temp[mi].pop_back();
					//break;
				}
				else
				{
					flag_erase = false;
					cluster_temp[mj].erase(cluster_temp[mj].begin() + i);
					i = i - 1;
				}
			}

			if (flag_erase == false) {
				if (cluster_temp[mj].size() == 0)
					cluster_temp.erase(cluster_temp.begin() + mj);
				dis_temp.clear();
				dis_temp.resize(cluster_temp.size());
				for (auto &i : dis_temp)
					i.resize(cluster_temp.size());

				for (int i = 0; i < cluster_temp.size(); ++i) {
					for (int j = 0; j < cluster_temp.size(); ++j) {
						if (i == j) {
							dis_temp[i][j] = std::numeric_limits<double>::max();
						}
						else
						{
							double sum_dis = 0.0;
							for (auto &k : cluster_temp[i]) {
								for (auto &l : cluster_temp[j]) {
									sum_dis += (dis_matrix[k][l] < dis_matrix[l][k] ? dis_matrix[k][l] : dis_matrix[l][k]);
								}
							}
							if (cluster_temp[i].size() == 0 || cluster_temp[j].size() == 0) {
								std::cout << "cluster_temp[i].size() == 0" << std::endl;
							}
							dis_temp[i][j] = sum_dis / (cluster_temp[i].size() * cluster_temp[j].size());

							/*double sum_dis = std::numeric_limits<double>::max();
							for (auto &k : cluster_temp[i]) {
								for (auto &l : cluster_temp[j]) {
									if (sum_dis > dis_matrix[k][l])
										sum_dis = dis_matrix[k][l];
								}
							}
							dis_temp[i][j] = sum_dis;*/
						}
					}
				}
			}
			else
			{
 				m_cluster.push_back({});
				data_temp = data;
				for (int i = 0; i < cluster_temp[mi].size(); ++i) {
					m_cluster.back().push_back(data_temp[cluster_temp[mi][i]]);
					auto it = std::find(data.begin(), data.end(), data_temp[cluster_temp[mi][i]]);
					data.erase(it);
				}
				cluster_temp.clear();
				cluster_temp.resize(data.size());
				for (int i = 0; i < data.size(); ++i) {
					cluster_temp[i].push_back(i);
				}
				dis_matrix.clear();
				dis_matrix.resize(data.size());
				for (auto &i : dis_matrix)
					i.resize(data.size());
				for (int i = 0; i < data.size(); ++i) {
					for (int j = 0; j < data.size(); ++j) {
						if (i == j)
						{
							dis_matrix[i][j] = std::numeric_limits<double>::max();
						}
						else
						{
							auto ii = std::distance(m_data.begin(), std::find(m_data.begin(), m_data.end(), data[i]));
							auto jj = std::distance(m_data.begin(), std::find(m_data.begin(), m_data.end(), data[j]));
							dis_matrix[i][j] = m_dis_matrix[ii][jj];
						}
					}
				}
				dis_temp = dis_matrix;
			}
		}

		m_cluster.push_back({});
		for (int i = 0; i < cluster_temp[0].size(); ++i) {
			m_cluster.back().push_back(data[cluster_temp[0][i]]);
		}
		result.member = m_cluster;

		for (int i = 0; i < m_cluster.size(); ++i) {
			for (int j = 0; j < m_cluster[i].size(); ++j) {
				for (int k = 0; k < m_cluster.size(); ++k) {
					if (i == k) continue;
					if (std::find(m_cluster[k].begin(), m_cluster[k].end(), m_cluster[i][j]) != m_cluster[k].end()) {
						std::cout << "The AHC clustering was wrong. There is same customer in different vehicle!" << std::endl;
					}
				}
			}
		}


		for (size_t i = 0; i < m_cluster.size(); ++i) {
			bool is_over_load = 0;
			double deliveryWeight = 0;
			for (int j = 0; j < m_cluster[i].size(); ++j) {
				deliveryWeight = deliveryWeight + DVRP_CAST(pro)->get_net().get_road_net()[m_cluster[i][j]]->datum_t.delivery_weight;
			}
			if (deliveryWeight > DVRP_CAST(pro)->get_vehicle_property()[2].capacity * m_ca_cons_co) {
				is_over_load = true;
			}
			else is_over_load = false;
		}

		cacul_objective(m_cluster);
	}
	
	const std::vector<std::vector<Real>>& AHC_v2::cacul_dis_matrix(const std::vector<size_t> &data, Problem *pro) {
		std::pair<Real, Real> d_min_max, tw_min_max, tw_min_max_1;
		Real d_conflict = 0.0, tw_conflict = 0.0;
		d_min_max.first = std::numeric_limits<double>::max();
		d_min_max.second = 0;
		tw_min_max.first = std::numeric_limits<double>::max();
		tw_min_max.second = 0;
		/*std::vector<std::vector<Real>> d_conflict_matrix;
		std::vector<std::vector<Real>> tw_conflict_matrix;
		std::vector<std::vector<Real>> dis_matrix;*/

		m_d_conflict_matrix.resize(data.size());
		m_tw_conflict_matrix.resize(data.size());
		//tw_conflict_matrix_1.resize(data.size());
		m_dis_matrix.resize(data.size());

		for (int i = 0; i < data.size(); ++i) {
			m_d_conflict_matrix[i].resize(data.size());
			m_tw_conflict_matrix[i].resize(data.size());
			//tw_conflict_matrix_1[i].resize(data.size());
			m_dis_matrix[i].resize(data.size());
		}
		for (int i = 0; i < data.size(); ++i) {
			for (int j = 0; j < data.size(); ++j) {
				if (i == j) {
					m_d_conflict_matrix[i][j] = std::numeric_limits<double>::max();
					m_tw_conflict_matrix[i][j] = std::numeric_limits<double>::max();
					//tw_conflict_matrix_1[i][j] = std::numeric_limits<double>::max();
				}
				else
				{
					d_conflict = cacul_d_conflict(data[i], data[j], pro);
					tw_conflict = cacul_tw_conflict(data[i], data[j], pro);
					//tw_conflict_1 = 1 / cacul_tw_conflict(data[i], data[j]);

					if (d_min_max.first > d_conflict)
						d_min_max.first = d_conflict;
					if (d_min_max.second < d_conflict)
						d_min_max.second = d_conflict;
					if (tw_min_max.first > tw_conflict)
						tw_min_max.first = tw_conflict;
					if (tw_min_max.second < tw_conflict)
						tw_min_max.second = tw_conflict;
					/*if (tw_min_max_1.first > tw_conflict_1)
						tw_min_max_1.first = tw_conflict_1;
					if (tw_min_max_1.second < tw_conflict_1)
						tw_min_max_1.second = tw_conflict_1;*/

					m_d_conflict_matrix[i][j] = d_conflict;
					m_tw_conflict_matrix[i][j] = tw_conflict;
					//tw_conflict_matrix_1[i][j] = tw_conflict_1;
				}
				/*if (d_min_max.first > d_conflict_matrix[i][j])
					d_min_max.first = d_conflict_matrix[i][j];
				if (d_min_max.second < d_conflict_matrix[i][j])
					d_min_max.second = d_conflict_matrix[i][j];
				if (tw_min_max.first > tw_conflict_matrix[i][j])
					tw_min_max.first = tw_conflict_matrix[i][j];
				if (tw_min_max.second < tw_conflict_matrix[i][j])
					tw_min_max.second = tw_conflict_matrix[i][j];*/
			}
		}

		for (int i = 0; i < data.size(); ++i) {
			for (int j = 0; j < data.size(); ++j) {
				if (i == j) m_dis_matrix[i][j] = std::numeric_limits<double>::max();
				else {
					d_conflict = (m_d_conflict_matrix[i][j] - d_min_max.first) / (d_min_max.second - d_min_max.first);
					
					tw_conflict = (m_tw_conflict_matrix[i][j] - tw_min_max.first) / (tw_min_max.second - tw_min_max.first);
					//tw_conflict_1 = (tw_conflict_matrix_1[i][j] - tw_min_max_1.first) / (tw_min_max_1.second - tw_min_max_1.first);
					if (tw_conflict > std::numeric_limits<double>::max()) {
						std::cout << "tw_conflict == inf" << std::endl;
					}
					if (d_conflict > std::numeric_limits<double>::max()) {
						std::cout << "d_conflict == inf" << std::endl;
					}
					m_d_conflict_matrix[i][j] = d_conflict;
					m_tw_conflict_matrix[i][j] = tw_conflict;
					//tw_conflict_matrix_1[i][j] = tw_conflict_1;
					m_dis_matrix[i][j] = m_weight * d_conflict + (1 - m_weight) * tw_conflict;//�˴�����Ҫ��Щ����
					//dis_matrix[i][j] = m_weight * d_conflict + (1 - m_weight) * tw_conflict_1;//�˷�ʽ�����ʣ�
					//if(p>=0.5)
					//	dis_matrix[i][j] = m_weight * d_conflict + (1 - m_weight) * tw_conflict;//�˴�����Ҫ��Щ����
					//else
					//{
					//	dis_matrix[i][j] = m_weight * d_conflict + (1 - m_weight) * tw_conflict_1;
					//}
				}
			}
		}
		return m_dis_matrix;
	}
	
	Real AHC_v2::cacul_d_conflict(size_t i, size_t j, Problem *pro) {
		return DVRP_CAST(pro)->get_shortest_dis(DVRP_CAST(pro)->get_depart_time(),i,j);
	}

	Real AHC_v2::cacul_tw_conflict(size_t i, size_t j, Problem *pro) {
		Real diff_tw = 0.0;
		size_t cus1 = i;
		size_t cus2 = j;

		Real cus1_beginTime = DVRP_CAST(pro)->get_net().get_road_net()[cus1]->datum_t.begin_time;
		Real cus1_endTime = DVRP_CAST(pro)->get_net().get_road_net()[cus1]->datum_t.end_time;
		Real cus2_beginTime = DVRP_CAST(pro)->get_net().get_road_net()[cus2]->datum_t.begin_time;
		Real cus2_endTime = DVRP_CAST(pro)->get_net().get_road_net()[cus2]->datum_t.end_time;
		Real beginTime = cus1_beginTime > cus2_beginTime ? cus1_beginTime : cus2_beginTime;
		Real endTime = cus1_endTime < cus2_endTime ? cus1_endTime : cus2_endTime;
		if (endTime >= beginTime) {
			diff_tw = endTime - beginTime;
			if (diff_tw <= 0)
				diff_tw = 0.1;
			diff_tw = diff_tw / (cus2_endTime - cus2_beginTime);
		}
		else
		{
			diff_tw = beginTime - endTime;
			diff_tw = 1 / diff_tw / (cus2_endTime - cus2_beginTime);
		}


		return diff_tw;
	}
	void AHC_v2::cacul_objective(const std::vector<std::vector<size_t>> &cluster)
	{
		Real sum = 0.0;
		Real sum_d = 0.0, sum_tw = 0.0;
		for (int i = 0; i < m_cluster.size(); ++i) {
			for (int j = 0; j < m_cluster[i].size(); ++j) {
				for (int k = 0; k < m_cluster[i].size(); ++k) {
					if (j == k) continue;
					auto ii = std::distance(m_data.begin(), std::find(m_data.begin(), m_data.end(), m_cluster[i][j]));
					auto jj = std::distance(m_data.begin(), std::find(m_data.begin(), m_data.end(), m_cluster[i][k]));
					Real dis_temp = m_dis_matrix[ii][jj] < m_dis_matrix[jj][ii] ? m_dis_matrix[ii][jj] : m_dis_matrix[jj][ii];
					sum += dis_temp;
					if (dis_temp == m_dis_matrix[ii][jj]) {
						sum_d += m_d_conflict_matrix[ii][jj];
						sum_tw += m_tw_conflict_matrix[ii][jj];
					}
					else
					{
						sum_d += m_d_conflict_matrix[jj][ii];
						sum_tw += m_tw_conflict_matrix[jj][ii];
					}
				}
			}
		}
		m_objective.push_back(sum_d);
		m_objective.push_back(sum_tw);
		m_objective.push_back(sum);
	}

	const bool AHC_v2::is_valid(const std::vector<size_t> &cluster, const std::vector<size_t> &data, Problem *pro) const{
		bool is_over_load = 0;
		double deliveryWeight = 0;
		for (size_t i = 0; i < cluster.size(); ++i) {

			deliveryWeight = deliveryWeight + DVRP_CAST(pro)->get_net().get_road_net()[data[cluster[i]]]->datum_t.delivery_weight;

		}
		if (deliveryWeight > DVRP_CAST(pro)->get_vehicle_property()[2].capacity * m_ca_cons_co) {
			is_over_load = true;
		}
		else is_over_load = false;
		return is_over_load;
	}
}