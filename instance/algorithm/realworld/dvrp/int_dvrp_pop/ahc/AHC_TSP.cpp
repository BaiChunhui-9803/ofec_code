#include "AHC_TSP.h"


namespace ofec {
	void AHC_TSP::clustering(Clustering_result &result, Problem *pro) {
		std::vector<std::vector<size_t>> cluster_temp(m_data.size());
		for (int i = 0; i < m_data.size(); ++i) {
			cluster_temp[i].push_back(i);
		}
		std::vector<std::vector<Real>> dis_temp;

		cacul_dis_matrix(dis_temp, cluster_temp, m_data, pro);
		auto data = m_data;
		auto data_temp = m_data;

		while (cluster_temp.size() != 1) {
			double minDt = DBL_MAX;
			size_t mi = 0, mj = 0;
			bool flag_erase = false;
			for (int i = 0; i < dis_temp.size(); ++i) {
				auto it = std::min_element(dis_temp[i].begin(), dis_temp[i].end());
				if (minDt > (*it)) {
					minDt = (*it);
					mi = i;
					mj = std::distance(dis_temp[i].begin(),it);
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
				cacul_dis_matrix(dis_temp, cluster_temp, data, pro);
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
				dis_temp.clear();
				cacul_dis_matrix(dis_temp, cluster_temp, data, pro);
			}
		}

		m_cluster.push_back({});
		for (int i = 0; i < cluster_temp[0].size(); ++i) {
			m_cluster.back().push_back(data[cluster_temp[0][i]]);
		}
		for (int i = 0; i < m_cluster.size(); ++i) {
			if (m_cluster[i].size() <= 1)continue;
			if (m_lkh == true)
				tsp_lkh(m_cluster[i], pro);
			else
				tsp_nearest(m_cluster[i], pro);
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
	}

	void AHC_TSP::cacul_dis_matrix(std::vector<std::vector<Real>> &dis_temp, 
		const std::vector<std::vector<size_t>> &cluster_temp, const std::vector<size_t> &data, Problem *pro) 
	{
		dis_temp.resize(cluster_temp.size());
		for (auto &i : dis_temp)
			i.resize(cluster_temp.size());
		std::vector<size_t> member;
		std::vector<std::vector<std::pair<Real, Real>>> d;
		d.resize(cluster_temp.size());
		for (auto &i : d)
			i.resize(cluster_temp.size());
		std::pair<Real, Real> d_min_max, tw_min_max;
		d_min_max.first = DBL_MAX;
		d_min_max.second = 0;
		tw_min_max.first = DBL_MAX;
		tw_min_max.second = 0;

		for (int i = 0; i < cluster_temp.size(); ++i) {
			for (int j = 0; j < cluster_temp.size(); ++j) {
				if (i == j) {
					dis_temp[i][j] = DBL_MAX;
					d[i][j].first = DBL_MAX;
					d[i][j].second = DBL_MAX;
				}
				else if(j > i)
				{
					member.clear();
					for (auto &k : cluster_temp[i])
						member.push_back(data[k]);
					
					for (auto &l : cluster_temp[j])
						member.push_back(data[l]);
					if(m_lkh == true)
						d[i][j]=tsp_lkh(member, pro);
					else
						d[i][j] = tsp_nearest(member, pro);//计算cluster i和cluster j之间的距离

					if (d_min_max.first > d[i][j].first)
						d_min_max.first = d[i][j].first;
					if (d_min_max.second < d[i][j].first)
						d_min_max.second = d[i][j].first;

					if (tw_min_max.first > d[i][j].second)
						tw_min_max.first = d[i][j].second;
					if (tw_min_max.second < d[i][j].second)
						tw_min_max.second = d[i][j].second;
				}
				else
				{
					d[i][j] = d[j][i];
				}
			}
		}

		for (int i = 0; i < dis_temp.size(); ++i) {
			for (int j = 0; j < dis_temp[i].size(); ++j) {
				if (i == j) continue;
				if(j>i)
				{
					if (d_min_max.first == d_min_max.second) {
						dis_temp[i][j] = 0;
					}
					else
					{
						dis_temp[i][j] = m_co_d * (d[i][j].first - d_min_max.first) / (d_min_max.second - d_min_max.first) + \
							(1 - m_co_d) * (d[i][j].second - tw_min_max.first) / (tw_min_max.second - tw_min_max.first);

					}
				}
				else if(j<i)
				{
					dis_temp[i][j] = dis_temp[j][i];
				}
				if (dis_temp[i][j] == NAN|| dis_temp[i][j] > DBL_MAX) {
					std::cout << "dis_temp[" << i << "][" << j << "] means nothing" << std::endl;
				}
			}
		}
	}
	
	std::pair<Real, Real> AHC_TSP::tsp_nearest(std::vector<size_t>& member_index, Problem *pro)
	{
		std::vector<size_t> tsp;
		std::vector<std::tuple<size_t, Real, Real>> member_cost;
		size_t start = DVRP_CAST(pro)->get_depot();
		Real Time = DVRP_CAST(pro)->get_depart_time();
		Real present_time = Time;
		size_t next = 0;
		std::pair<Real, Real>conflict_d_tw;
		Real sum_conflict_d = 0.0, sum_conflict_tw = 0.0;
		int cnt = 0;
		while (!member_index.empty()) {
			member_cost.clear();
			std::vector<Real> conflict_tw(member_index.size());
			std::vector<Real> conflict_d(member_index.size());
			std::vector<Real> conflict(member_index.size());
			for (int i = 0; i < member_index.size(); ++i) {
				if (member_index[i] == start) continue;
				present_time = Time;

				get_cost(start, member_index[i], present_time, pro);//此处耗时可能太长

				member_cost.push_back(std::make_tuple(i, present_time - Time, present_time));

				conflict_d[i] = present_time - Time;

				Real diff_tw = 0.0;
				Real beginTime = DVRP_CAST(pro)->get_net().get_road_net()[member_index[i]]->datum_t.begin_time;
				Real endTime = DVRP_CAST(pro)->get_net().get_road_net()[member_index[i]]->datum_t.end_time;
				if (beginTime >= present_time) {
					diff_tw = beginTime - present_time;
				}
				else if (beginTime < present_time && endTime >= present_time)
				{
					diff_tw = 0;
				}
				else if (present_time > endTime)
				{
					diff_tw = present_time - endTime;
				}
				conflict_tw[i] = diff_tw;

			}
			Real conflict_d_min = *std::min_element(conflict_d.begin(), conflict_d.end());
			Real conflict_d_max = *std::max_element(conflict_d.begin(), conflict_d.end());
			Real conflict_tw_min = *std::min_element(conflict_tw.begin(), conflict_tw.end());
			Real conflict_tw_max = *std::max_element(conflict_tw.begin(), conflict_tw.end());

			for (int i = 0; i < member_index.size(); ++i) {
				if (member_index.size() <= 1) {
					conflict_tw_min = 0;
					conflict_d_min = 0;
				}
				conflict[i] = m_co_d * ((conflict_d[i] - conflict_d_min) / (conflict_d_max - conflict_d_min)) + \
					(1 - m_co_d) * ((conflict_tw[i] - conflict_tw_min) / (conflict_tw_max - conflict_tw_min));
			}

			auto index = std::distance(conflict.begin(), std::min_element(conflict.begin(), conflict.end()));
			next = member_index[index];
			if (start == DVRP_CAST(pro)->get_depot()) {
				Time = std::get<2>(member_cost[index]);
			}
				//Time = present_time;
			else {
				Time = std::get<2>(member_cost[index]) + DVRP_CAST(pro)->get_net().get_road_net()[member_index[index]]->server_time;
			}
				//Time = present_time + DVRP_CAST(pro)->get_net().get_road_net()[DVRP_CAST(pro)->get_current_cus_order()[member_index[index]]]->server_time;
			start = next;
			cnt++;
			sum_conflict_d += conflict_d[index];
			sum_conflict_tw += conflict_tw[index];
			tsp.push_back(next);
			auto it_check=std::find(member_index.begin(), member_index.end(), *(member_index.begin() + index));
			if (it_check == member_index.end()) {
				system("pause");
			}
			member_index.erase(member_index.begin() + index);
		}
		member_index = tsp;
		conflict_d_tw.first = sum_conflict_d;
		conflict_d_tw.second = sum_conflict_tw;
		return conflict_d_tw;
	}

	std::pair<Real, Real> AHC_TSP::tsp_lkh(std::vector<size_t>& member_index, Problem *pro)
	{
		std::vector<size_t> tsp;
		//std::vector<std::tuple<size_t, Real, Real>> member_cost;
		size_t start = DVRP_CAST(pro)->get_depot();
		Real Time = DVRP_CAST(pro)->get_depart_time();
		Real present_time = Time;
		std::pair<Real, Real>conflict_d_tw;
		Real sum_conflict_d = 0.0, sum_conflict_tw = 0.0;

		auto member_t = member_index;
 		member_t.insert(member_t.begin(), start);
		std::string tsp_filename;
		tsp_filename = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/" + "LKH/temp file/cluster_temp.atsp";
		std::ostringstream Totsp;
		std::ofstream clusterTotsp(tsp_filename);
		Totsp << "NAME : " + tsp_filename + "\n" + "COMMENT : ATSP" + "\n" + "TYPE : ATSP" + "\n" + "DIMENSION : " + std::to_string(member_t.size()) + "\n" + "EDGE_WEIGHT_TYPE : EXPLICIT" + "\n"\
			+ "EDGE_WEIGHT_FORMAT : FULL_MATRIX" + "\n" + "EDGE_WEIGHT_SECTION" + "\n";
		
		std::vector<std::vector<Real>> dis_matrix;

		cal_tw_conflict(member_t,dis_matrix, pro);

		for (int j = 0; j < member_t.size(); ++j) {
			for (int k = 0; k < member_t.size(); ++k) {
				Totsp << dis_matrix[j][k] * 100000 << " ";
				/*if (j == k) {
					Totsp << "9999999 ";
				}
				else {

					Totsp << DVRP_CAST(pro)->get_cost_time(present_time, member_t[j], member_t[k]) << " ";
				}*/
			}
			Totsp << "\n";
		}
		Totsp << "EOF";
		clusterTotsp << Totsp.str();
		clusterTotsp.close();

		std::string par_filename = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/" + "LKH/temp file/cluster_temp.par";
		std::ostringstream Topar;
		std::ofstream par(par_filename);
		Topar << "PROBLEM_FILE = "<<"cluster_temp.atsp" <<"\n" <<"MOVE_TYPE = 5" << "\n" << "PATCHING_C = 3" << "\n" << "PATCHING_A = 2" << "\n" << "RUNS = 10" << "\n" << "TOUR_FILE = "<< static_cast<std::string>(g_working_dir) << "/instance/algorithm/realworld/DVRP/LKH/temp file/cluster_result.txt"<< "\n";
		par << Topar.str();
		par.close();
		//cus_order[i] = tsp(cus_order[i]);

		ostringstream FileName;
		FileName << static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/" + "LKH/temp file/cluster_temp.par";
		
		lkh.set_test(9999);
		lkh.run(FileName.str(), tsp);

		member_index.clear();
		for (int i = 0; i < tsp.size(); ++i) {
			if (tsp[i] - 1 == 0) {
				member_index.push_back(start);
				continue;
			}
			member_index.push_back(member_t[tsp[i] - 1]);
		}

		//get conflict
		Real diff_tw = 0.0;
		Real beginTime = 0;
		Real endTime = 0;
		for (int i = 0; i < member_index.size()-1; ++i) {
			//if (member_index[i] == start) continue;
			//present_time = Time;

			get_cost(member_index[i], member_index[i + 1], present_time, pro);

			//member_cost.push_back(std::make_tuple(i, present_time - Time, present_time));

			sum_conflict_d += (present_time - Time); 

			diff_tw = 0.0; 
			beginTime = DVRP_CAST(pro)->get_net().get_road_net()[member_index[i + 1]]->datum_t.begin_time;
			endTime = DVRP_CAST(pro)->get_net().get_road_net()[member_index[i + 1]]->datum_t.end_time;
			if (beginTime >= present_time) {
				diff_tw = beginTime - present_time;
			}
			else if (beginTime < present_time && endTime >= present_time)
			{
				diff_tw = 0;
			}
			else if (present_time > endTime)
			{
				diff_tw = present_time - endTime;
			}
			sum_conflict_tw += diff_tw;
			Time = present_time;
		}
		member_index.erase(member_index.begin());
		conflict_d_tw.first = sum_conflict_d;
		conflict_d_tw.second = sum_conflict_tw;
		
		return conflict_d_tw;
	}

	
	void AHC_TSP::cal_tw_conflict(const std::vector<size_t> &member_index, std::vector<std::vector<Real>> &dis_matrix, Problem *pro)
	{
		//get conflict
		//member_index.insert(member_index.begin(),DVRP_CAST(pro)->get_depot());
		dis_matrix.clear();
		dis_matrix.resize(member_index.size());
		for (auto &i : dis_matrix)
			i.resize(member_index.size());
		Real Time = DVRP_CAST(pro)->get_depart_time();
		Real present_time = Time;
		Real diff_tw = 0.0;
		Real beginTime = 0;
		Real endTime = 0;
		//Real conflict_d = 0, conflict_tw = 0;
		std::vector<std::vector<Real>> conflict_tw(member_index.size());
		std::vector<std::vector<Real>> conflict_d(member_index.size());
		for (auto &i : conflict_tw) {
			i.resize(member_index.size());
		}
		for (auto &i : conflict_d) {
			i.resize(member_index.size());
		}
		for (int i = 0; i < member_index.size(); ++i) {
			//if (member_index[i] == start) continue;
			//present_time = Time;
			for (int j = 0; j < member_index.size(); ++j) {
				if (i == j) {
					//dis_matrix[i][j] = 99999;
					conflict_d[i][j] = 9999;
					conflict_tw[i][j] = 9999;
					continue;
				}
				get_cost(member_index[i], member_index[j], present_time, pro);
				conflict_d[i][j] = present_time - Time;
				diff_tw = 0.0;
				
				beginTime = DVRP_CAST(pro)->get_net().get_road_net()[member_index[j]]->datum_t.begin_time;
				endTime = DVRP_CAST(pro)->get_net().get_road_net()[member_index[j]]->datum_t.end_time;
				
				if (beginTime >= present_time) {
					diff_tw = beginTime - present_time;
				}
				else if (beginTime < present_time && endTime >= present_time)
				{
					diff_tw = 0;
				}
				else if (present_time > endTime)
				{
					diff_tw = present_time - endTime;
				}
				conflict_tw[i][j] = diff_tw;
				present_time = Time;
			}
			//dis_matrix[i][j] = conflict_d[i][j] + conflict_tw[i][j];
		}
		std::pair<Real, Real>d_min_max;
		std::pair<Real, Real>tw_min_max;
		d_min_max.first = 9999;
		d_min_max.second = 0;
		tw_min_max.first = 9999;
		tw_min_max.second = 0;
		for (int i = 0; i < conflict_tw.size(); ++i) {
			for (int j = 0; j < conflict_tw[i].size(); ++j) {
				if (tw_min_max.first > conflict_tw[i][j])
					tw_min_max.first = conflict_tw[i][j];
				if (tw_min_max.second < conflict_tw[i][j])
					tw_min_max.second = conflict_tw[i][j];

				if (d_min_max.first > conflict_d[i][j])
					d_min_max.first = conflict_d[i][j];
				if (d_min_max.second < conflict_d[i][j])
					d_min_max.second = conflict_d[i][j];
			}
		}
		for (int i = 0; i < dis_matrix.size(); ++i) {
			for (int j = 0; j < dis_matrix[i].size(); ++j) {
				if (i == j) {
					dis_matrix[i][j] = 9999;
					continue;
				}
				if (d_min_max.first == d_min_max.second) {
					dis_matrix[i][j] = 0;
				}
				else
				{
					dis_matrix[i][j] = m_co_d * (conflict_d[i][j] - d_min_max.first) / (d_min_max.second - d_min_max.first) + \
						(1 - m_co_d) * (conflict_tw[i][j] - tw_min_max.first) / (tw_min_max.second - tw_min_max.first);

				}
				if (dis_matrix[i][j] == NAN || dis_matrix[i][j] > DBL_MAX) {
					std::cout << "dis_temp[" << i << "][" << j << "] means nothing" << std::endl;
				}
			}
		}
			
	}

	void AHC_TSP::get_cost(size_t s, size_t e, Real &present_time, Problem *pro)
	{
		//present_time = present_time + DVRP_CAST(pro)->get_cost_time(present_time, s, e);
		//present_time = present_time + DVRP_CAST(pro)->get_cost_time()[s][e];
		present_time = present_time + DVRP_CAST(pro)->get_cost_time(present_time, s, e);
	}

	const bool AHC_TSP::is_valid(const std::vector<size_t> &cluster, const std::vector<size_t> &data, Problem *pro) const {
		/*bool is_over_load = false;
		double deliveryWeight = 0;
		for (size_t i = 0; i < cluster.size(); ++i) {

			deliveryWeight = deliveryWeight + DVRP_CAST(pro)->get_net().get_road_net()[data[cluster[i]]]->datum_t.delivery_weight;

		}
		if (deliveryWeight > DVRP_CAST(pro)->get_vehicle_property()[2].capacity * m_ca_cons_co) {
			is_over_load = true;
		}
		else is_over_load = false;
		return is_over_load;*/

		bool is_over_load = false;
		double deliveryWeight = 0;
		double over_load = 0;
		for (size_t i = 0; i < cluster.size(); ++i) {
			if (DVRP_CAST(pro)->get_net().get_road_net()[data[cluster[i]]]->datum_t.is_served == false)
				deliveryWeight += DVRP_CAST(pro)->get_net().get_road_net()[data[cluster[i]]]->datum_t.delivery_weight;
		}
		if (deliveryWeight > DVRP_CAST(pro)->get_vehicle_property()[2].capacity * m_ca_cons_co) {
			over_load += deliveryWeight - DVRP_CAST(pro)->get_vehicle_property()[2].capacity * m_ca_cons_co;
		}
		for (size_t i = 0; i < cluster.size(); ++i) {
			if (DVRP_CAST(pro)->get_net().get_road_net()[data[cluster[i]]]->datum_t.is_served == false)
				deliveryWeight = deliveryWeight - DVRP_CAST(pro)->get_net().get_road_net()[data[cluster[i]]]->datum_t.delivery_weight + DVRP_CAST(pro)->get_net().get_road_net()[data[cluster[i]]]->datum_t.pick_weight;
			else
				deliveryWeight = deliveryWeight - 0 + DVRP_CAST(pro)->get_net().get_road_net()[data[cluster[i]]]->datum_t.pick_weight;

			if (deliveryWeight > DVRP_CAST(pro)->get_vehicle_property()[2].capacity * m_ca_cons_co) {
				over_load += deliveryWeight - DVRP_CAST(pro)->get_vehicle_property()[2].capacity * m_ca_cons_co;
			}
		}
		if (over_load > 0) {
			is_over_load = true;
		}
		else is_over_load = false;
		return is_over_load;
	}
}