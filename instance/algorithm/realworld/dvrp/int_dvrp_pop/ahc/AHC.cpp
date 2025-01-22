#include "AHC.h"
namespace ofec {
	void AHC::clustering(Clustering_result &result)
	{
		auto dis_matrix = cacul_dis_matrix(m_data);
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
			size_t mi=0, mj=0;
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
				if (is_valid(cluster_temp[mi], data)) {
					flag_erase = true;
					break;
				}
				cluster_temp[mi].push_back(cluster_temp[mj][i]);
			}
			//if (cluster_temp.size() == 1) break;
			if (flag_erase ==false) {
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
				//is_valid(cluster_temp[mi], data);
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
				dis_matrix = cacul_dis_matrix(data);
				dis_temp = dis_matrix;
			}
		}

		m_cluster.push_back({});
		for (int i = 0; i < cluster_temp[0].size(); ++i) {
			m_cluster.back().push_back(data_temp[cluster_temp[0][i]]);
		}
		result.member = m_cluster;
		
		for (size_t i = 0; i < m_cluster.size(); ++i) {
			bool is_over_load = 0;
			double deliveryWeight = 0;
			for (int j = 0; j < m_cluster[i].size(); ++j) {
				deliveryWeight = deliveryWeight + DVRP_CAST->get_net().get_road_net()[m_cluster[i][j]]->datum_t.delivery_weight;
			}
			if (deliveryWeight > DVRP_CAST->get_vehicle_property()[2].capacity*0.75) {
				is_over_load = true;
			}
			else is_over_load = false;
		}
		
		//cacul_objective(m_cluster);
	}
	std::vector<std::vector<Real>> AHC::cacul_dis_matrix(const std::vector<size_t> &data)
	{
		std::pair<Real, Real> d_min_max, tw_min_max, tw_min_max_1;
		Real d_conflict = 0.0, tw_conflict = 0.0, tw_conflict_1 = 0.0;
		d_min_max.first = std::numeric_limits<double>::max();
		d_min_max.second = 0;
		tw_min_max.first = std::numeric_limits<double>::max();
		tw_min_max.second = 0;
		tw_min_max_1.first = std::numeric_limits<double>::max();
		tw_min_max_1.second = 0;
		std::vector<std::vector<Real>> d_conflict_matrix;
		std::vector<std::vector<Real>> tw_conflict_matrix;
		//std::vector<std::vector<Real>> tw_conflict_matrix_1;
		std::vector<std::vector<Real>> dis_matrix;

		d_conflict_matrix.resize(data.size());
		tw_conflict_matrix.resize(data.size());
		//tw_conflict_matrix_1.resize(data.size());
		dis_matrix.resize(data.size());

		for (int i = 0; i < data.size(); ++i) {
			d_conflict_matrix[i].resize(data.size());
			tw_conflict_matrix[i].resize(data.size());
			//tw_conflict_matrix_1[i].resize(data.size());
			dis_matrix[i].resize(data.size());
		}
		for (int i = 0; i < data.size(); ++i) {
			for (int j = 0; j < data.size(); ++j) {
				if (i == j) {
					d_conflict_matrix[i][j] = std::numeric_limits<double>::max();
					tw_conflict_matrix[i][j] = std::numeric_limits<double>::max();
					//tw_conflict_matrix_1[i][j] = std::numeric_limits<double>::max();
				}
				else
				{
					d_conflict = cacul_d_conflict(data[i], data[j]);
					tw_conflict = cacul_tw_conflict(data[i], data[j]);
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

					d_conflict_matrix[i][j] = d_conflict;
					tw_conflict_matrix[i][j] = tw_conflict;
					//tw_conflict_matrix_1[i][j] = tw_conflict_1;
				}
			}
		}

		for (int i = 0; i < data.size(); ++i) {
			for (int j = 0; j < data.size(); ++j) {
				if (i == j) dis_matrix[i][j] = std::numeric_limits<double>::max();
				else {
					d_conflict = (d_conflict_matrix[i][j] - d_min_max.first) / (d_min_max.second - d_min_max.first);
					if(tw_min_max.first==0&&tw_min_max.first==tw_min_max.second){
						tw_conflict = 0;
					}
					else
					{
						tw_conflict = (tw_conflict_matrix[i][j] - tw_min_max.first) / (tw_min_max.second - tw_min_max.first);
					}
					//tw_conflict_1 = (tw_conflict_matrix_1[i][j] - tw_min_max_1.first) / (tw_min_max_1.second - tw_min_max_1.first);
					d_conflict_matrix[i][j] = d_conflict;
					tw_conflict_matrix[i][j] = tw_conflict;
					dis_matrix[i][j] = m_weight * d_conflict + (1 - m_weight) * tw_conflict;//�˴�����Ҫ��Щ����
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
		return dis_matrix;
	}
	Real AHC::cacul_d_conflict(size_t i, size_t j)
	{
		return DVRP_CAST->get_shortest_dis(DVRP_CAST->get_depart_time(),i,j);
	}
	Real AHC::cacul_tw_conflict(size_t i, size_t j)
	{
		Real diff_tw = 0.0;
		size_t cus1 = i;
		size_t cus2 = j;

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
			diff_tw = diff_tw / (cus2_endTime - cus2_beginTime);
		}
		else
		{
			diff_tw = beginTime - endTime;
			diff_tw = 1 / diff_tw / (cus2_endTime - cus2_beginTime);
		}

		//Real cus1_beginTime = DVRP_CAST->get_net().get_road_net()[cus1]->datum_t.begin_time;
		//Real cus1_endTime = DVRP_CAST->get_net().get_road_net()[cus1]->datum_t.end_time;
		//Real cus2_beginTime = DVRP_CAST->get_net().get_road_net()[cus2]->datum_t.begin_time;
		//Real cus2_endTime = DVRP_CAST->get_net().get_road_net()[cus2]->datum_t.end_time;

		//auto t = DVRP_CAST->get_cost_time()[cus1][cus2];
		//Real arrive_time = cus1_beginTime + DVRP_CAST->get_cost_time()[cus1][cus2] + DVRP_CAST->get_net().get_road_net()[cus1]->server_time;
		////Real arrive_time = cus1_beginTime + DVRP_CAST->get_cost_time()[cus1][cus2];
		//
		//if (arrive_time < cus2_beginTime) {
		//	diff_tw = (cus2_beginTime - arrive_time)/(cus2_endTime- cus2_beginTime);
		//}
		//if (arrive_time > cus2_endTime) {
		//	diff_tw = (arrive_time - cus2_endTime) / (cus2_endTime - cus2_beginTime);
		//}
		//if (arrive_time >= cus2_beginTime && arrive_time <= cus2_endTime) {
		//	diff_tw = 0;
		//}

		return diff_tw;
	}
	//void AHC::cacul_objective(std::vector<std::vector<size_t>> &cluster)
	//{
	//	m_objective.resize(2);
	//	std::vector<Real> objective;
	//	Real sum = 0.0;
	//	Real sum_d = 0.0, sum_tw = 0.0;
	//	//std::vector<std::vector<size_t>> cluster;
	//	for (int i = 0; i < cluster.size(); ++i) {
	//		for (int j = 0; j < cluster[i].size(); ++j) {
	//			for (int k = 0; k < cluster[i].size(); ++k) {
	//				if (j == k) continue;
	//				Real dis_temp = cacul_d_conflict(cluster[i][j],cluster[i][k]) < cacul_d_conflict(cluster[i][k], cluster[i][j]) ? cacul_d_conflict(cluster[i][j], cluster[i][k]) : cacul_d_conflict(cluster[i][k], cluster[i][j]);
	//				sum_d += dis_temp;
	//				if (dis_temp == cacul_d_conflict(cluster[i][j], cluster[i][k])) {
	//					sum_tw += cacul_tw_conflict(cluster[i][j], cluster[i][k]);
	//				}
	//				else
	//				{
	//					sum_tw += cacul_tw_conflict(cluster[i][k], cluster[i][j]);
	//				}
	//			}
	//		}
	//	}
	//	sum = m_weight*sum_d + (1-m_weight)*sum_tw;
	//	m_objective.push_back(sum_d);
	//	m_objective.push_back(sum_tw);
	//	objective.push_back(sum);
	//}

	bool AHC::is_valid(const std::vector<size_t> &cluster, const std::vector<size_t> &data) {
		bool is_over_load = 0;
		double deliveryWeight = 0;
		for (size_t i = 0; i < cluster.size(); ++i) {

			deliveryWeight = deliveryWeight + DVRP_CAST->get_net().get_road_net()[data[cluster[i]]]->datum_t.delivery_weight;

		}
		if (deliveryWeight > DVRP_CAST->get_vehicle_property()[2].capacity*m_ca_cons_co) {
			is_over_load = true;
		}
		else is_over_load = false;
		return is_over_load;
	}
}