#include "fcm_twl.h"
namespace ofec {

	fcm_twl::fcm_twl(int num_cluster, const std::vector<size_t> & Data_points) :m_num_clusters(num_cluster) {
		m_data_points = Data_points;
		m_num_data_points = m_data_points.size();
		m_dim = 2;
	}
	void fcm_twl::get_result(fcm_twl::Clustering_result &result)
	{
		initialize_membership(result);
	}

	void fcm_twl::initialize_membership(fcm_twl::Clustering_result & res_temp)
	{
		const Real pi = 3.14159265;
		Real s = 0.0;
		Real r = 0.02;
		Real sum = 0.0;
		Real angle = 360 / m_num_clusters;
		std::pair<Real, Real> depot;
		depot.first = DVRP_CAST->get_net().get_road_net()[DVRP_CAST->get_depot()]->datum_t.x;
		depot.second = DVRP_CAST->get_net().get_road_net()[DVRP_CAST->get_depot()]->datum_t.y;
		m_cluster_centre.resize(m_num_clusters);
		for (int i = 0; i < m_num_clusters; ++i) {
			m_cluster_centre[i].push_back(depot.first + r * cos((angle*i)*pi / 180));
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
		std::vector<std::vector<Real>> membership_d(m_num_data_points);
		for (auto &i : membership_d)
			i.resize(m_num_clusters);
		m_degree_of_memb.resize(m_num_data_points);
		for (int i = 0; i < m_num_data_points; ++i) {
			m_degree_of_memb[i].resize(m_num_clusters);
			sum = 0.0;
			s = 0.0;
			for (int j = 0; j < m_num_clusters; ++j) {
				if (m_data_points[i] == m_cluster_centre_index[j])
					s = std::numeric_limits<double>::min();
				else
					s = DVRP_CAST->get_shortest_dis(DVRP_CAST->get_depart_time(), m_data_points[i], m_cluster_centre_index[j]);
				sum += s;
			}
			for (int k = 0; k < m_num_clusters; ++k) {
				if (m_data_points[i] == m_cluster_centre_index[k])
					membership_d[i][k] = std::numeric_limits<double>::min();
				else
					membership_d[i][k] = (1 / DVRP_CAST->get_shortest_dis(DVRP_CAST->get_depart_time(), m_data_points[i], m_cluster_centre_index[k])) / sum;
			}
		}

		m_degree_of_memb = membership_d;
		get_member(res_temp);

		std::vector<std::vector<Real>> membership_tw(m_num_data_points);
		std::vector<std::vector<Real>> membership_l(m_num_data_points);
		for (auto &i : membership_tw)
			i.resize(m_num_clusters);
		for (auto &i : membership_l)
			i.resize(m_num_clusters);
		int max_diff = 0;
		std::vector<std::pair<size_t, size_t>> index_archive;
		std::vector<size_t> tsp_member;
		std::vector<size_t> tsp_temp;
		Real total_load = 0.0;
		for (int i = 0; i < m_num_data_points; ++i) {
			total_load += DVRP_CAST->get_net().get_road_net()[m_data_points[i]]->datum_t.delivery_weight;
		}



		/*std::vector<size_t> m_cus_order;
		m_cus_order.push_back(DVRP_CAST->get_depot());
		for (int j = 0; j < res_temp.member[0].size(); ++j) {
			m_cus_order.push_back(m_data_points[res_temp.member[0][j]]);
		}
		std::string tsp_filename = "D:/MyCloud/OneDrive/ʵ���ҹ���/�㱨ppt/20190529/tsp_lkh.atsp";
		std::ofstream tsp_lkh(tsp_filename);
		tsp_lkh << "NAME : " + tsp_filename + "\n" + "COMMENT : ATSP" + "\n" + "TYPE : ATSP" + "\n" + "DIMENSION : " + std::to_string(m_cus_order.size()) + "\n" + "EDGE_WEIGHT_TYPE : EXPLICIT" + "\n"\
			+ "EDGE_WEIGHT_FORMAT : FULL_MATRIX" + "\n" + "EDGE_WEIGHT_SECTION" + "\n";
		for (int j = 0; j < m_cus_order.size(); ++j) {
			for (int k = 0; k < m_cus_order.size(); ++k) {
				if (j == k) {
					tsp_lkh << "9999999 ";
				}
				else {

					tsp_lkh << DVRP_CAST->get_shortest_dis()[m_cus_order[j]][m_cus_order[k]] << " ";
				}
			}
			tsp_lkh << "\n";
		}
		tsp_lkh << "EOF";
		tsp_lkh.close();
		std::string par_filename;
		par_filename = "D:/MyCloud/OneDrive/ʵ���ҹ���/�㱨ppt/20190529/tsp_lkh.par";
		std::ofstream par(par_filename);
		par << "PROBLEM_FILE = ./DVRP/tsp_lkh.atsp" << std::endl;
		par << "MOVE_TYPE = 5" << std::endl;
		par << "PATCHING_C = 3" << std::endl;
		par << "PATCHING_A = 2" << std::endl;
		par << "RUNS = 10" << std::endl;
		par << "TOUR_FILE = ./DVRP/result_lkh.txt" << std::endl;
		par.close();*/


		do {
			for (int i = 0; i < res_temp.member.size(); ++i) {
				for (int j = 0; j < res_temp.member[i].size(); ++j) {
					Real sum_tw = 0.0;
					index_archive.clear();
					
					for (int z = 0; z < res_temp.member.size(); ++z) {
						//DWORD star_time_item = GetTickCount();
						if (res_temp.member[z].size() == 0)
							continue;

						tsp_temp.clear();
						tsp_member.clear();
						if (i == z) {
							tsp_temp = res_temp.member[z];
						}
						else
						{
							tsp_temp = res_temp.member[z];
							tsp_temp.push_back(res_temp.member[i][j]);
						}

						//DWORD star_time = GetTickCount();
						tsp_member = tsp(tsp_temp);
						//DWORD end_time = GetTickCount();
						//std::cout << std::endl << "TSP: " << (end_time - star_time) << "ms" << std::endl;

						auto it = std::find(tsp_member.begin(), tsp_member.end(), res_temp.member[i][j]);
						if (it == tsp_member.end()) std::cout << "can not find element" << std::endl;
						auto it2 = tsp_member.begin();
						if (it == tsp_member.begin()) {
							it2 = std::find(res_temp.member[z].begin(), res_temp.member[z].end(), *(it + 1));
						}
						else
						{
							it2 = std::find(res_temp.member[z].begin(), res_temp.member[z].end(), *(it - 1));
						}
						auto index = std::distance(res_temp.member[z].begin(), it2);


						/*int index = 0;
						do {
							index = ofec::global::ms_global->m_uniform[ofec::caller::Problem]->next_non_standard<size_t>(0, res_temp.member[z].size());
						} while (res_temp.member[z][index] == res_temp.member[i][j]);*/

						index_archive.push_back(std::make_pair(z,index));

						//size_t cus_1 = m_data_points[res_temp.member[i][j]];
						//size_t cus_2 = m_data_points[res_temp.member[z][index]];

						//Real cus1_beginTime = DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.begin_time;
						//Real cus1_endTime = DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.end_time; 
						//Real cus2_beginTime = DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.begin_time;
						//Real cus2_endTime = DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.end_time;
						//Real beginTime = cus1_beginTime > cus2_beginTime ? cus1_beginTime : cus2_beginTime;
						//Real endTime = cus1_endTime < cus2_endTime ? cus1_endTime : cus2_endTime;
						//Real diff_tw = endTime - beginTime;
						////Real diff_tw = (cus1_endTime - cus1_beginTime) / (endTime - beginTime);
						//sum_tw += abs(diff_tw);

						/*Real diff_tw = DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.begin_time + (DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.end_time - DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.begin_time) \
							- (DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.begin_time + (DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.end_time - DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.begin_time));
						sum_tw += pow(abs(diff_tw) - DVRP_CAST->get_average_cost(), 2);*/
					
						//DWORD end_time_item = GetTickCount();
						//std::cout << std::endl << "ITEM: " << (end_time_item - star_time_item) << "ms" << std::endl;
					}
					

					for (int k = 0; k < index_archive.size(); ++k) {
						if (res_temp.member[index_archive[k].first].size() == 0) {
							membership_tw[res_temp.member[i][j]][index_archive[k].first] = 1;
						}
						else
						{
							size_t cus_1 = m_data_points[res_temp.member[i][j]];
							size_t cus_2 = m_data_points[res_temp.member[index_archive[k].first][index_archive[k].second]];
							
							Real cus1_beginTime = DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.begin_time;
							Real cus1_endTime = DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.end_time;
							Real cus2_beginTime = DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.begin_time;
							Real cus2_endTime = DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.end_time;
							Real beginTime = cus1_beginTime > cus2_beginTime ? cus1_beginTime : cus2_beginTime;
							Real endTime = cus1_endTime < cus2_endTime ? cus1_endTime : cus2_endTime;
							if (endTime >= beginTime) {
								Real diff_tw = endTime - beginTime;
								if (diff_tw <= 0)
									diff_tw = 0.1;
								membership_tw[res_temp.member[i][j]][index_archive[k].first] = 1 / diff_tw / (cus1_endTime - cus1_beginTime);
								//membership_tw[res_temp.member[i][j]][index_archive[k].first] = 1 / diff_tw / sum_tw;
							}
							else
							{
								Real diff_tw = beginTime - endTime;
								membership_tw[res_temp.member[i][j]][index_archive[k].first] = diff_tw / (cus1_endTime - cus1_beginTime);
								//membership_tw[res_temp.member[i][j]][index_archive[k].first] = diff_tw / sum_tw;
							}
							/*Real diff_tw = (cus1_endTime - cus1_beginTime) / (endTime - beginTime);
							membership_tw[res_temp.member[i][j]][index_archive[k].first] = 1 / diff_tw / sum_tw;*/
							
							/*Real diff_tw = DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.begin_time + (DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.end_time - DVRP_CAST->get_net().get_road_net()[cus_1]->datum_t.begin_time) \
								- (DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.begin_time + (DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.end_time - DVRP_CAST->get_net().get_road_net()[cus_2]->datum_t.begin_time));
							if (diff_tw == DVRP_CAST->get_average_cost())
								membership_tw[res_temp.member[i][j]][index_archive[k].first] = 1;
							else
							{
								membership_tw[res_temp.member[i][j]][index_archive[k].first] = 1 / pow(abs(diff_tw) - DVRP_CAST->get_average_cost(), 2) / sum_tw;
							}*/
						}
					}
					for (int x = 0; x < m_num_clusters; ++x) {
						Real load = 0.0;
						for (auto &k : res_temp.member[x]) {
							load += DVRP_CAST->get_net().get_road_net()[m_data_points[k]]->datum_t.delivery_weight;
						}
						if (load == 0)
							membership_l[res_temp.member[i][j]][x] = 1;
						else
							membership_l[res_temp.member[i][j]][x] = 1 / load / total_load;
						//membership_l[res_temp.member[i][j]][x] = 1 / Real(res_temp.member[x].size()) / Real(m_num_data_points);
					}

				}
			}
			//��һ��
			for (auto &i : membership_tw) {
				for (auto &j : i) {
					int distur = ofec::global::ms_global->m_uniform[ofec::caller::Problem]->next_non_standard<size_t>(0, 10);
					j = j * (1 + distur * 0.01);
				}
			}
			Real ms_d_min = std::numeric_limits<double>::max();
			Real ms_d_max = -1;
			Real ms_tw_min = std::numeric_limits<double>::max();
			Real ms_tw_max = -1;
			Real ms_l_min = std::numeric_limits<double>::max();
			Real ms_l_max = -1;
			for (int i = 0; i < m_num_data_points; ++i) {
				if (ms_d_max <= *std::max_element(membership_d[i].begin(), membership_d[i].end()))
					ms_d_max = *std::max_element(membership_d[i].begin(), membership_d[i].end());
				if (ms_d_min >= *std::min_element(membership_d[i].begin(), membership_d[i].end()))
					ms_d_min = *std::min_element(membership_d[i].begin(), membership_d[i].end());

				 
				if (ms_tw_max <= *std::max_element(membership_tw[i].begin(), membership_tw[i].end()))
					ms_tw_max = *std::max_element(membership_tw[i].begin(), membership_tw[i].end());
				if (ms_tw_min >= *std::min_element(membership_tw[i].begin(), membership_tw[i].end()))
					ms_tw_min = *std::min_element(membership_tw[i].begin(), membership_tw[i].end());

				if (ms_l_max <= *std::max_element(membership_l[i].begin(), membership_l[i].end()))
					ms_l_max = *std::max_element(membership_l[i].begin(), membership_l[i].end());
				if (ms_l_min >= *std::min_element(membership_l[i].begin(), membership_l[i].end()))
					ms_l_min = *std::min_element(membership_l[i].begin(), membership_l[i].end());
			}
			for (int i = 0; i < m_num_data_points; ++i) {
				for (int j = 0; j < m_num_clusters; ++j) {
					membership_d[i][j]= (membership_d[i][j] - ms_d_min) / (ms_d_max - ms_d_min);
					membership_tw[i][j] = (membership_tw[i][j] - ms_tw_min) / (ms_tw_max - ms_tw_min);
					membership_l[i][j] = (membership_l[i][j] - ms_l_min) / (ms_l_max - ms_l_min);
				}
			}

			for (int i = 0; i < m_num_data_points; ++i) {
				for (int j = 0; j < m_num_clusters; ++j) {
					//m_degree_of_memb[i][j] = membership_tw[i][j];
					m_degree_of_memb[i][j] = membership_d[i][j] + membership_tw[i][j];
					//m_degree_of_memb[i][j] = membership_d[i][j] + membership_tw[i][j] + 0.01*membership_l[i][j];
				}
			}
			res_temp.member.clear();
			get_member(res_temp);
			int member_max = 0;
			int member_min = std::numeric_limits<int>::max();
			for (auto &i : res_temp.member) {
				if (member_max < i.size())
					member_max = i.size();
				if (member_min > i.size())
					member_min = i.size();
			}
			//max_diff = member_max - member_min;
			max_diff ++;
		}while (max_diff < 1);

		
	}

	void fcm_twl::get_member(Clustering_result &result)
	{
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

	std::vector<size_t> fcm_twl::tsp(std::vector<size_t>& member_index)
	{
		std::vector<bool> member_visited(member_index.size(),false);
		std::vector<size_t> tsp;
		std::vector<std::tuple<size_t, Real, Real>> member_cost;
		size_t start = DVRP_CAST->get_depot();
		Real Time = DVRP_CAST->get_depart_time();
		Real present_time = Time;
		size_t next = 0;
		while(!member_index.empty()){
		//while (tsp.size() != member_index.size()) {
			member_cost.clear();
			for (int i = 0; i < member_index.size(); ++i) {
				if (m_data_points[member_index[i]] == start) continue;
				//if (member_visited[i] == true) continue;
				present_time = Time;

				//DWORD start_time = GetTickCount();
				get_cost(start, m_data_points[member_index[i]], present_time);//�˴���ʱ����̫��
				//DWORD end_time = GetTickCount();
				//std::cout << std::endl << "get_cost(): " << (end_time - start_time) << "ms" << std::endl;

				member_cost.push_back(std::make_tuple(i, present_time - Time, present_time));
				//member_cost.push_back(std::make_tuple(i, DVRP_CAST->get_shortest_dis()[start][m_data_points[member_index[i]]], present_time));
			}
			Real sum_cost = 0.0;
			Real sum_tw = 0.0;
			for (int i = 0; i < member_cost.size();++i) {
				sum_cost += std::get<1>(member_cost[i]);//
				Real diff_tw = 0.0;
				Real beginTime = DVRP_CAST->get_net().get_road_net()[m_data_points[member_index[std::get<0>(member_cost[i])]]]->datum_t.begin_time;
				Real endTime = DVRP_CAST->get_net().get_road_net()[m_data_points[member_index[std::get<0>(member_cost[i])]]]->datum_t.end_time;
				if (beginTime >= std::get<2>(member_cost[i])) {
					diff_tw = beginTime - std::get<2>(member_cost[i]);
				}
				else if(beginTime < std::get<2>(member_cost[i]) && endTime >= std::get<2>(member_cost[i]))
				{
					diff_tw = 0;
				}
				else if(std::get<2>(member_cost[i]) > endTime)
				{
					diff_tw = std::get<2>(member_cost[i]) - endTime;
				}
				sum_tw += diff_tw;
			}
			std::vector<Real> member_ship_tw(member_cost.size());
			std::vector<Real> member_ship_cost(member_cost.size());
			std::vector<Real> member_ship(member_cost.size());
			for (int i = 0; i < member_cost.size(); ++i) {
				Real beginTime = DVRP_CAST->get_net().get_road_net()[m_data_points[member_index[std::get<0>(member_cost[i])]]]->datum_t.begin_time;
				Real endTime = DVRP_CAST->get_net().get_road_net()[m_data_points[member_index[std::get<0>(member_cost[i])]]]->datum_t.end_time;
				if (beginTime >= std::get<2>(member_cost[i])) 
				{
					member_ship_tw[i] = 1 / (beginTime - std::get<2>(member_cost[i])) / sum_tw;
				}
				else if (beginTime < std::get<2>(member_cost[i]) && endTime >= std::get<2>(member_cost[i]))
				{
					member_ship_tw[i] = 1;
				}
				else if (std::get<2>(member_cost[i]) > endTime)
				{
					member_ship_tw[i] = 1 / (std::get<2>(member_cost[i]) - endTime) / sum_tw; 
				}

				member_ship_cost[i] = 1 / std::get<1>(member_cost[i]) / sum_cost;
			}
			//��һ��
			Real ms_cost_min = *std::min_element(member_ship_cost.begin(),member_ship_cost.end());
			Real ms_cost_max = *std::max_element(member_ship_cost.begin(), member_ship_cost.end());
			Real ms_tw_min = *std::min_element(member_ship_tw.begin(), member_ship_tw.end());
			Real ms_tw_max = *std::max_element(member_ship_tw.begin(), member_ship_tw.end());
			
			for (int i = 0; i < member_cost.size(); ++i) {
				if (member_cost.size() <= 1) {
					ms_tw_min = 0;
					ms_cost_min = 0;
				}
				member_ship[i] = ((member_ship_cost[i] - ms_cost_min) / (ms_cost_max - ms_cost_min));
				//member_ship[i] = ((member_ship_cost[i] - ms_cost_min) / (ms_cost_max - ms_cost_min)) + ((member_ship_tw[i] - ms_tw_min) / (ms_tw_max - ms_tw_min));
			}
			auto index = std::distance(member_ship.begin(), std::max_element(member_ship.begin(), member_ship.end()));
			next = member_index[std::get<0>(member_cost[index])];
			member_index.erase(member_index.begin() + std::get<0>(member_cost[index]));
			//member_index.erase(member_index.begin() + index);
			//member_visited[std::get<0>(member_cost[index])] = true;
			if (start == DVRP_CAST->get_depot())
				Time = std::get<2>(member_cost[index]);
			else
				Time = std::get<2>(member_cost[index]) + DVRP_CAST->get_net().get_road_net()[DVRP_CAST->get_current_cus_order()[0]]->server_time;
			start = m_data_points[next];
			tsp.push_back(next);
			
		}
		return tsp;
	}

	void fcm_twl::get_cost(size_t s, size_t e,Real &present_time)
	{
		std::vector<size_t> path;
		if (present_time >= 1440) present_time = 1440;
		//DWORD start_time = GetTickCount();
		DVRP_CAST->shortest_path(s, e, path, present_time);
		//DWORD end_time = GetTickCount();
		//std::cout << std::endl << "Shortest_path: " << (end_time - start_time) << "ms" << std::endl;
		Real sum = 0.0;
		for (int i = 0; i < path.size() - 1; ++i) {
			sum += DVRP_CAST->get_net().get_node_dis()[path[i]][path[i + 1]];
			present_time = present_time + (DVRP_CAST->get_net().get_node_dis()[path[i]][path[i + 1]] / DVRP_CAST->get_net().get_road_net()[i + 1]->datum_t.v[static_cast<int>(present_time)] * 60);
		}
		//return sum;
	}

}