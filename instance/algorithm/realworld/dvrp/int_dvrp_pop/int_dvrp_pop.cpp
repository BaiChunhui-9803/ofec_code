#include "int_dvrp_pop.h"
#include <thread>

namespace ofec {
	void INIT_DVRP_POP::ahc_tsp_init_pop(std::vector<DVRP_ind> &pop, Problem *pro, Random *rnd) {
		//auto num_core = thread::hardware_concurrency() - 6;//获取cpu核心个数  
		//
		//for (int i = 0; i < pop.size(); ) {
		//	system("pause");
		//	std::vector<thread> init_pop;
		//	for (int j = 0; j < num_core; ++j) {
		//		if (i >= pop.size())
		//			break;
		//		init_pop.push_back(std::thread(&INIT_DVRP_POP::ahc_tsp_init_pop_, this, std::ref(pop[i]), i));
		//		i++;
		//	}
		//	for (auto &k : init_pop)
		//		k.join();
		//	init_pop.clear();
		//}

		for (int i = 0; i < pop.size();++i ) {
			ahc_tsp_init_pop_(pop[i], i, pro, rnd);
		}
		//for (int i = 0; i < pop.size(); ++i) {//加并行？
		//	std::cout << "initialize pop " << std::to_string(i) << std::endl;
		//	auto &x = pop[i]->variable();
		//	auto &obj = pop[i]->objective();
		//	//随机生成权重
		//	Real weight = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<Real>(0, 1);
		//	//DWORD start_time = GetTickCount();
		//	
		//	AHC_TSP::Clustering_result result;
		//	AHC_TSP ahc(weight, DVRP_CAST(pro)->get_current_cus_order());
		//	ahc.clustering(result);

		//	//DWORD end_time = GetTickCount();
		//	//std::cout << std::endl << "The run time of cluster is:" << (end_time - start_time) << "ms!" << std::endl;

		//	std::vector<std::vector<size_t>> m_cus_order;
		//	m_cus_order = result.member;
		//
		//	for (int j = 0; j < m_cus_order.size(); ++j) {
		//		x.m_cus_order.push_back({});
		//		x.m_cus_order.back().push_back(DVRP_CAST(pro)->get_depot());
		//		for (int k = 0; k < m_cus_order[j].size(); ++k) {
		//			x.m_cus_order.back().push_back(m_cus_order[j][k]);
		//		}
		//		x.m_cus_order.back().push_back(DVRP_CAST(pro)->get_depot());
		//	}

		//	x.m_vehicle_type.resize(x.m_cus_order.size());
		//	for (size_t z = 0; z < x.m_vehicle_type.size(); z++) {
		//		x.m_vehicle_type[z] = 2;
		//	}
		//	DVRP_CAST(pro)->construct_path(*pop[i]);
			
		//}
		
	}
	void INIT_DVRP_POP::ahc_tsp_init_pop_(DVRP_ind& ind, int i, Problem *pro, Random *rnd) {
		std::cout << "initialize pop " << std::to_string(i) << std::endl;
		auto &x = ind->variable();
		//随机生成权重
		Real weight = rnd->uniform.next();
		

		AHC_TSP::Clustering_result result;
		AHC_TSP ahc(weight, DVRP_CAST(pro)->get_current_cus_order());
		ahc.clustering(result, pro);


		std::vector<std::vector<size_t>> m_cus_order;
		m_cus_order = result.member;

		for (int j = 0; j < m_cus_order.size(); ++j) {
			x.m_cus_order.push_back({});
			x.m_cus_order.back().push_back(DVRP_CAST(pro)->get_depot());
			for (int k = 0; k < m_cus_order[j].size(); ++k) {
				x.m_cus_order.back().push_back(m_cus_order[j][k]);
			}
			x.m_cus_order.back().push_back(DVRP_CAST(pro)->get_depot());
		}

		x.m_vehicle_type.resize(x.m_cus_order.size());
		for (size_t z = 0; z < x.m_vehicle_type.size(); z++) {
			x.m_vehicle_type[z] = 2;
		}
		DVRP_CAST(pro)->construct_path(*ind);
	}

	void INIT_DVRP_POP::randomly_init_pop(std::vector<DVRP_ind>& pop, Problem *pro, Random *rnd) {
		for (int i = 0; i < pop.size(); ++i) {
			std::cout << "initialize pop " << std::to_string(i) << std::endl;
			DVRP_CAST(pro)->initializeSolution(*pop[i], rnd);
		}
	}

	void INIT_DVRP_POP::ahc_initialize_pop(std::vector<DVRP_ind>& pop, Problem *pro, Random *rnd) {
		for (int i = 0; i < pop.size(); ++i) {
			auto &x = pop[i]->variable();
			Real weight = rnd->uniform.next();

			AHC_v2::Clustering_result result;
			AHC_v2 ahc(weight, DVRP_CAST(pro)->get_current_cus_order());
			ahc.clustering(result, pro);

			std::vector<std::vector<size_t>> m_cus_order;
			m_cus_order = result.member;

			calcu_cus_order(m_cus_order, pro);

			for (int j = 0; j < m_cus_order.size(); ++j) {
				x.m_cus_order.push_back({});
				x.m_cus_order.back().push_back(DVRP_CAST(pro)->get_depot());
				for (int k = 0; k < m_cus_order[j].size(); ++k) {
					x.m_cus_order.back().push_back(m_cus_order[j][k]);
				}
				x.m_cus_order.back().push_back(DVRP_CAST(pro)->get_depot());
			}
			x.m_vehicle_type.resize(x.m_cus_order.size());
			for (size_t z = 0; z < x.m_vehicle_type.size(); z++) {
				x.m_vehicle_type[z] = 2;
			}
			DVRP_CAST(pro)->construct_path(*pop[i]);
		}
	}

	void INIT_DVRP_POP::calcu_cus_order(std::vector<std::vector<size_t>>& m_cus_order, Problem *pro) {
		for (int i = 0; i < m_cus_order.size(); ++i) {
			m_cus_order[i] = tsp(m_cus_order[i], pro);
		}
	}
	
	std::vector<size_t> INIT_DVRP_POP::tsp(std::vector<size_t>& member_index, Problem *pro) {
		std::vector<size_t> tsp;
		std::vector<std::tuple<size_t, Real, Real>> member_cost;
		size_t start = DVRP_CAST(pro)->get_depot();
		Real Time = DVRP_CAST(pro)->get_depart_time();
		Real present_time = Time;
		size_t next = 0;
		std::pair<Real, Real>conflict_d_tw;
		//Real sum_conflict_d = 0.0, sum_conflict_tw = 0.0;
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
					(1-m_co_d) * ((conflict_tw[i] - conflict_tw_min) / (conflict_tw_max - conflict_tw_min));
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
			tsp.push_back(next);
			member_index.erase(member_index.begin() + index);
		}
		member_index = tsp;
		return tsp;
	}

	void INIT_DVRP_POP::get_cost(size_t s, size_t e, Real &present_time, Problem *pro)
	{
		//present_time = present_time + DVRP_CAST(pro)->get_cost_time(present_time,s,e);
		//present_time = present_time + DVRP_CAST(pro)->get_cost_time()[s][e];
		present_time = present_time + DVRP_CAST(pro)->get_cost_time(present_time, s, e);
	}








}