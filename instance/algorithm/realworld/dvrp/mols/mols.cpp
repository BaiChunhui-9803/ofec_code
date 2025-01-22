#include "MOLS.h"
//#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include <numeric>

using namespace std::placeholders;

#ifdef OFEC_DEMO
#include "../../../../../../ui/buffer/realworld/DVRP/buffer_DVRP.h"
extern unique_ptr<ofec_demo::scene> ofec_demo::msp_buffer;
#endif

namespace ofec {
	void MOLS_pop::initialize() {
		//test_lkh();
		//m_ahc_tsp.ahc_tsp_init_pop(m_individuals);
		m_ahc_tsp.randomly_init_pop(m_individuals);
		m_pop_size = m_individuals.size();
		evaluate();
		for (auto &i : m_individuals) {
			i->set_vio_obj(0, 0);
		}
		sort();
		int rank_tmp = 0;
		int ind_id = 0;
		while (1) {
			for (; ind_id < m_individuals.size(); ) {
				if (!check_order(*m_individuals[ind_id]))
					std::cout << "order is not enough" << "\n";
				if (m_individuals[ind_id]->fitness() == rank_tmp) {
					m_archive.emplace_back(new Solution_type(*m_individuals[ind_id]));
				}
				ind_id++;
				if (m_archive.size() == DVRP_CAST->objective_size())
					break;
			}
			if (m_archive.size() == DVRP_CAST->objective_size())
				break;
			else {
				m_archive.clear();
				ind_id = 0;
				rank_tmp++;
			}
			if (rank_tmp >= m_individuals.size()) {
				m_archive.clear();
			}
		}
		if (m_archive.size() != DVRP_CAST->objective_size()) {
			std::cout << "elements of Archive is not enough!" << "\n";
			system("pause");
		}
		m_individuals.clear();
		for (int i = 0; i < m_archive.size(); ++i) {
			m_individuals.emplace_back(new Solution_type(*m_archive[i]));
		}
		for (auto &i : m_individuals) {
			i->evaluate();
		}
		/****************conference paper experiment********************************/
		//记录订单，对比订单是否相同
		std::string current_order_Fname;
		current_order_Fname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/current_order" + to_string(DVRP_CAST->get_current_cus_order().size()) + "MOLS_MOVRPRTC_AHC_initial.txt";
		std::ofstream Corder_(current_order_Fname);
		std::stringstream Corder_str;
		for (auto &i : DVRP_CAST->get_current_cus_order()) {
			Corder_str << i << " ";
		}
		Corder_ << Corder_str.str();
		Corder_.close();
		/**************************************************************************/
	}

	MOLS_pop::Solution_type & MOLS_pop::select_sol(size_t obj_id)
	{
		/*Real max_obj = 0.0;
		int sol_id = 0;
		for (int i = 0; i < m_archive.size(); ++i) {
			if (max_obj <= m_archive[i]->objective()[obj_id]) {
				max_obj = m_archive[i]->objective()[obj_id];
				sol_id = i;
			}
		}*/
		int sol_id = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, m_archive.size());
		if (!check_order(*m_archive[sol_id]))
			std::cout << "order is not enough" << "\n";
		Solution_type &sol_tmp = *m_archive[sol_id];
		return sol_tmp;
	}


	void MOLS_pop::sort() {
		std::vector<int> ranks_;
		std::vector<Solution_type*> pop_;
		std::function<ofec::dominationship(Solution_type* const&, Solution_type* const&)> comp = std::bind(&MOLS_pop::Pareto_compare, this, _1, _2);
		for (auto& i : m_individuals)
			pop_.emplace_back(i.get());
		nd_sort::fast_sort<Solution_type*>(pop_, ranks_, comp);
		for (size_t i = 0; i < m_individuals.size(); i++) {
			m_individuals[i]->set_rank(ranks_[i]);
		}
	}

	void MOLS_pop::Online()
	{
		if (DVRP_CAST->objective_size() == 3)
			std::cout << "Offline WT: " << m_best_ind_Offline.objective()[2] << "\n";
		if (DVRP_CAST->objective_size() == 1)
			std::cout << "Offline WT: " << m_best_ind_Offline.constraint_value()[3] << "\n";
		DVRP_CAST->set_Online(true);

		m_best_ind_Online = m_best_ind_Offline;

		//以一个最优方案为基础，处理新的订单。最优订单的定义为个指标归一化之后

		for (int i = 0; i < m_best_ind_Online.variable().m_cus_order.size(); ++i) {
			double presentTime = DVRP_CAST->get_depart_time() + DVRP_CAST->get_cost_time(DVRP_CAST->get_depart_time(), m_best_ind_Online.variable().m_cus_order[i][0], m_best_ind_Online.variable().m_cus_order[i][1]);
			for (int j = 1; j < m_best_ind_Online.variable().m_cus_order[i].size() - 1; ++j) {
				if (DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][j]]->datum_t.is_served == true)
					continue;
				DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][j]]->datum_t.is_served = true;
				//updateSeqOnline(m_best_ind_Online, i);
				auto customer = m_best_ind_Online.variable().m_cus_order[i][j];
				auto presentTime_curOrder = presentTime;
				auto presentTime_newOrder = presentTime;
				std::pair<int, int> new_order;
				m_insert_success = false;
				for (int k = 0; k < DVRP_CAST->get_future_cus_order().size(); ++k) {
					if (presentTime >= DVRP_CAST->get_future_cus_order()[k].second) {
						new_order = DVRP_CAST->get_future_cus_order()[k];
						if (DVRP_CAST->get_net().get_road_net()[new_order.first]->type != DVRP::node_type::future_customer || new_order.first == 175) {
							std::cout << "error" << "\n";
						}
						size_t next_cur_order = m_best_ind_Online.variable().m_cus_order[i][j + 1];
						presentTime = presentTime + DVRP_CAST->get_net().get_road_net()[customer]->server_time + DVRP_CAST->get_cost_time(presentTime, customer, new_order.first);
						if (presentTime >= DVRP_CAST->get_net().get_road_net()[new_order.first]->datum_t.begin_time && presentTime <= DVRP_CAST->get_net().get_road_net()[new_order.first]->datum_t.end_time) {

							m_best_ind_Online.variable().m_cus_order[i].insert(m_best_ind_Online.variable().m_cus_order[i].begin() + j + 1, new_order.first);
							DVRP_CAST->get_future_cus_order().erase(DVRP_CAST->get_future_cus_order().begin() + k);
							m_insert_success = true;
							if (m_best_ind_Online.variable().m_cus_order[i][0] != DVRP_CAST->get_depot()) {
								std::cout << "error" << "\n";
							}
							updateSeqOnline_inserted(m_best_ind_Online, i, new_order);

							presentTime_newOrder = presentTime;
							presentTime = presentTime + DVRP_CAST->get_net().get_road_net()[new_order.first]->server_time;
							presentTime = presentTime + DVRP_CAST->get_cost_time(presentTime, new_order.first, next_cur_order);

							/*double delay_time = 0;
							if (presentTime > DVRP_CAST->get_net().get_road_net()[next_cur_order]->datum_t.end_time) {
								delay_time = presentTime - DVRP_CAST->get_net().get_road_net()[next_cur_order]->datum_t.end_time;

								auto it= std::find(m_best_ind_Offline.variable().m_cus_order[i].begin(), m_best_ind_Offline.variable().m_cus_order[i].end(), next_cur_order);
								auto nextOrder_pos = std::distance(std::begin(m_best_ind_Offline.variable().m_cus_order[i]), it);
								if (delay_time > m_best_ind_Offline.variable().m_hi.delay_time[i][nextOrder_pos - 1]) {
									m_best_ind_Online.variable().m_cus_order[i].erase(m_best_ind_Online.variable().m_cus_order[i].begin() + j + 1);
									DVRP_CAST->get_future_cus_order().push_back(new_order);
									m_insert_success = false;
								}
							}*/
							/*if (m_insert_success == true) {
								auto it = std::find(m_best_ind_Offline.variable().m_cus_order[i].begin() + 1, m_best_ind_Offline.variable().m_cus_order[i].end(), next_cur_order);
								auto nextOrder_pos = std::distance(std::begin(m_best_ind_Offline.variable().m_cus_order[i]), it);

								if (presentTime > m_best_ind_Offline.variable().m_hi.arrive_time[i][nextOrder_pos-1]) {
									auto time_t = presentTime - m_best_ind_Offline.variable().m_hi.arrive_time[i][nextOrder_pos - 1];
									if (time_t > 60) {
										m_best_ind_Online.variable().m_cus_order[i].erase(m_best_ind_Online.variable().m_cus_order[i].begin() + j + 1);
										DVRP_CAST->get_future_cus_order().push_back(new_order);
										m_insert_success = false;
									}
								}
							}*/
							if (m_insert_success == true) {
								bool is_over_load = false;
								double deliveryWeight = 0;
								double over_load = 0;
								for (size_t m = 0; m < m_best_ind_Online.variable().m_cus_order[i].size(); ++m) {
									if (DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.is_served == false)
										deliveryWeight += DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.delivery_weight;
								}
								if (deliveryWeight > DVRP_CAST->get_vehicle_property()[2].capacity) {
									over_load += deliveryWeight - DVRP_CAST->get_vehicle_property()[2].capacity;
								}
								for (size_t m = 0; m < m_best_ind_Online.variable().m_cus_order[i].size(); ++m) {
									if (DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.is_served == false)
										deliveryWeight = deliveryWeight - DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.delivery_weight + DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.pick_weight;
									else
										deliveryWeight = deliveryWeight - 0 + DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.pick_weight;

									if (deliveryWeight > DVRP_CAST->get_vehicle_property()[2].capacity) {
										over_load += deliveryWeight - DVRP_CAST->get_vehicle_property()[2].capacity;
									}
								}
								if (over_load > 0) {
									m_best_ind_Online.variable().m_cus_order[i].erase(m_best_ind_Online.variable().m_cus_order[i].begin() + j + 1);
									DVRP_CAST->get_future_cus_order().push_back(new_order);
									m_insert_success = false;
								}
							}
						}
						if (m_insert_success == false)
							presentTime = presentTime_curOrder;
						if (m_insert_success == true) {
							//updateSeqOnline(m_best_ind_Online, i);
							presentTime = presentTime_newOrder;
							break;
						}
						//std::cout << m_best_ind_Online.objective()[0] << "\t" << m_best_ind_Online.constraint_value()[3] << "\t" << m_best_ind_Online.constraint_value()[2] << "\t" << m_best_ind_Online.constraint_value()[0] << "\t" << m_best_ind_Online.constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << m_best_ind_Online.variable().m_cus_order.size() << std::endl;

					}

				}
				if (m_insert_success == false) {
					presentTime = presentTime_curOrder;
					presentTime = presentTime + DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][j]]->server_time + DVRP_CAST->get_cost_time(presentTime, m_best_ind_Online.variable().m_cus_order[i][j], m_best_ind_Online.variable().m_cus_order[i][j + 1]);
				}
			}
			//标记已被服务
			for (int m = 1; m < m_best_ind_Online.variable().m_cus_order[i].size() - 1; ++m) {
				DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.is_served = false;
			}
			for (int n = 1; n < m_best_ind_Online.variable().m_cus_order[i].size() - 1; ++n) {
				if (DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][n]]->type == DVRP::node_type::future_customer)
					break;
				DVRP_CAST->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][n]]->datum_t.is_served = true;
			}
			DVRP_CAST->set_Online(true);
			construct_path_Online(m_best_ind_Online, i);
		}
		m_best_ind_Online.evaluate();
		if (DVRP_CAST->objective_size() == 1)
			std::cout << "Online WT: " << m_best_ind_Online.constraint_value()[3] << "\n";
		if (DVRP_CAST->objective_size() == 3)
			std::cout << "Online WT: " << m_best_ind_Online.objective()[2] << "\n";
	}

	void MOLS_pop::construct_path_Online(Solution_type & ind, size_t i)
	{
		ind.variable().m_path[i].clear();
		ind.variable().m_path_node_type[i].clear();

		int cnt = 0;
		std::vector<std::vector<size_t>> path_temp(ind.variable().m_cus_order[i].size() - 1);
		std::vector<std::vector<DVRP::node_type>> path_node_type_temp(ind.variable().m_cus_order[i].size() - 1);
		//path_temp[0].push_back(m_depot);
		Real present_time = DVRP_CAST->get_depart_time();
		for (size_t j = 0; j < ind.variable().m_cus_order[i].size() - 1; j++) {
			std::vector<size_t> Shortest_path;

			DVRP_CAST->shortest_path(ind.variable().m_cus_order[i][j], ind.variable().m_cus_order[i][j + 1], Shortest_path, present_time);//时间要变


			if (DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::current_customer || DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::future_customer) {
				present_time = present_time + DVRP_CAST->get_cost_time(present_time, ind.variable().m_cus_order[i][j], ind.variable().m_cus_order[i][j + 1]) + DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->server_time;
			}
			if (DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::depot) {
				present_time = present_time + DVRP_CAST->get_cost_time(present_time, ind.variable().m_cus_order[i][j], ind.variable().m_cus_order[i][j + 1]);
			}
			auto tt = DVRP_CAST->get_cost_time(present_time, ind.variable().m_cus_order[i][j], ind.variable().m_cus_order[i][j + 1]);
			Shortest_path.pop_back();
			path_temp[j] = Shortest_path;
			for (auto &m : path_temp[j]) {
				ind.variable().m_path[i].push_back(m);
				path_node_type_temp[j].push_back(DVRP::node_type::general_node);
			}
			path_node_type_temp[j][0] = DVRP_CAST->get_net().get_road_net()[path_temp[j][0]]->type;
			for (auto &m : path_node_type_temp[j]) {
				ind.variable().m_path_node_type[i].push_back(m);
			}
		}
		ind.variable().m_path[i].push_back(DVRP_CAST->get_depot());
		ind.variable().m_path_node_type[i].push_back(DVRP_CAST->get_net().get_road_net()[DVRP_CAST->get_depot()]->type);
	}

	void MOLS_pop::updateSeqOnline(Solution_type & ind, size_t i)
	{
		DVRP_CAST->set_Online(false);

		DVRP_CAST->construct_path(ind);
		ind.evaluate();

		m_score_Online.clear();
		m_score_Online.resize(m_num_operators_Online, 1);
		m_p_select_Online.clear();
		m_p_select_Online.resize(m_num_operators_Online);
		m_cp_Online.clear();
		m_cp_Online.resize(m_num_operators_Online);

		m_best_inds.clear();
		m_best_inds.emplace_back(new Solution_type(ind));
		calculate_initial_max_violation(m_best_inds);
		m_e = m_max_G;
		m_k = 0;
		std::cout << "updateSeqOnline() before: " << "\n";
		if (DVRP_CAST->objective_size() == 1)
			std::cout << ind.objective()[0] << "\t" << ind.constraint_value()[3] << "\t" << ind.constraint_value()[2] << "\t" << ind.constraint_value()[0] << "\t" << ind.constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		if (DVRP_CAST->objective_size() == 3)
			std::cout << ind.objective()[0] << "\t" << ind.objective()[2] << "\t" << ind.objective()[1] << "\t" << ind.constraint_value()[0] << "\t" << ind.constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;

		int iter_cnt = 0;
		while (iter_cnt < m_OnlineIters) {
			Online_ALS(ind, i);
			iter_cnt++;
		}
		std::cout << "updateSeqOnline() end: " << "\n";
		if (DVRP_CAST->objective_size() == 1)
			std::cout << ind.objective()[0] << "\t" << ind.constraint_value()[3] << "\t" << ind.constraint_value()[2] << "\t" << ind.constraint_value()[0] << "\t" << ind.constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		if (DVRP_CAST->objective_size() == 3)
			std::cout << ind.objective()[0] << "\t" << ind.objective()[2] << "\t" << ind.objective()[1] << "\t" << ind.constraint_value()[0] << "\t" << ind.constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;

		DVRP_CAST->set_Online(true);
	}

	void MOLS_pop::updateSeqOnline_inserted(Solution_type & ind, size_t i, std::pair<int, int> new_order)
	{
		DVRP_CAST->set_Online(false);
		DVRP_CAST->construct_path(ind);
		ind.evaluate();

		m_score_Online.clear();
		m_score_Online.resize(m_num_operators_Online, 1);
		m_p_select_Online.clear();
		m_p_select_Online.resize(m_num_operators_Online);
		m_cp_Online.clear();
		m_cp_Online.resize(m_num_operators_Online);

		m_best_inds.clear();
		m_best_inds.emplace_back(new Solution_type(ind));
		calculate_initial_max_violation(m_best_inds);
		m_e = m_max_G;
		m_k = 0;

		int iter_cnt = 0;
		//auto is_dominated = Pareto_compare(&m_best_ind_Online, &best_ind_Online_tmp);//m_best_ind_Offline
		double max_delayTime = 0;
		auto averge_delayTime = (std::numeric_limits<Real>::max)();
		std::cout << "Before Online ALS" << "\n";
		if (DVRP_CAST->objective_size() == 3) {
			std::cout << ind.objective()[0] << "\t" << ind.objective()[2] << "\t" << ind.objective()[1] << "\t" << ind.constraint_value()[0] << "\t" << ind.constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		}
		if (DVRP_CAST->objective_size() == 1) {
			std::cout << ind.objective()[0] << "\t" << ind.constraint_value()[3] << "\t" << ind.constraint_value()[2] << "\t" << ind.constraint_value()[0] << "\t" << ind.constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		}
		while (iter_cnt < m_OnlineIters) {
			Online_ALS(ind, i);
			iter_cnt++;
			/*max_delayTime = 0;
			for (auto &i : m_best_ind_Online.variable().m_hi.delay_time) {
				if (max_delayTime <= *std::max_element(i.begin(), i.end()))
					max_delayTime = *std::max_element(i.begin(), i.end());
			}*/
			if (DVRP_CAST->objective_size() == 3) {
				averge_delayTime = ind.objective()[1] / ind.variable().m_cus_order.size();//(DVRP_CAST->get_current_cus_order().size() + (DVRP_CAST->get_num_futureOrder() - DVRP_CAST->get_future_cus_order().size()));;//m_best_ind_Online.constraint_value()[2] / (DVRP_CAST->get_current_cus_order().size() + (DVRP_CAST->get_num_futureOrder() - DVRP_CAST->get_future_cus_order().size()));
			}
			if (DVRP_CAST->objective_size() == 1) {
				averge_delayTime = ind.constraint_value()[2] / ind.variable().m_cus_order.size();//(DVRP_CAST->get_current_cus_order().size() + (DVRP_CAST->get_num_futureOrder() - DVRP_CAST->get_future_cus_order().size()));;//m_best_ind_Online.constraint_value()[2] / (DVRP_CAST->get_current_cus_order().size() + (DVRP_CAST->get_num_futureOrder() - DVRP_CAST->get_future_cus_order().size()));
			}
			if (averge_delayTime > 45 && iter_cnt > 0.9 * m_OnlineIters) {
				if (DVRP_CAST->get_num_futureOrder() - DVRP_CAST->get_future_cus_order().size() == 0)
					break;
				auto it = std::find(ind.variable().m_cus_order[i].begin(), ind.variable().m_cus_order[i].end(), new_order.first);
				auto newOrder_pos = std::distance(std::begin(ind.variable().m_cus_order[i]), it);
				ind.variable().m_cus_order[i].erase(ind.variable().m_cus_order[i].begin() + newOrder_pos);
				DVRP_CAST->get_future_cus_order().push_back(new_order);
				auto fSize = DVRP_CAST->get_future_cus_order().size();
				if (ind.variable().m_cus_order[i][0] == ind.variable().m_cus_order[i][1]) {
					std::cout << "error" << "\n";
				}
				DVRP_CAST->construct_path(ind);
				ind.evaluate();
				m_insert_success = false;
				break;
			}
			if (averge_delayTime <= 45) {//&& iter_cnt > 0.5 * m_OnlineIters
				m_insert_success = true;
				std::cout << "Insert Success" << "\n";
				break;
			}
		}
		if (averge_delayTime > 45) {
			m_insert_success = false;
			std::cout << "Insert Fail" << "\n";
		}
		if (DVRP_CAST->objective_size() == 3) {
			std::cout << ind.objective()[0] << "\t" << ind.objective()[2] << "\t" << ind.objective()[1] << "\t" << ind.constraint_value()[0] << "\t" << ind.constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		}
		if (DVRP_CAST->objective_size() == 1) {
			std::cout << ind.objective()[0] << "\t" << ind.constraint_value()[3] << "\t" << ind.constraint_value()[2] << "\t" << ind.constraint_value()[0] << "\t" << ind.constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		}
		DVRP_CAST->set_Online(true);
	}

	void MOLS_pop::Online_ALS(Solution_type & ind, size_t i)
	{
		//adaptive local search
		Real sum_score = 0;
		for (auto &j : m_score_Online)
			sum_score += j;
		for (int j = 0; j < m_num_operators_Online; ++j) {
			m_p_select_Online[j] = Real(Real(m_score_Online[j]) / sum_score);
		}
		for (auto &i : m_cp_Online)
			i = 0;
		m_cp_Online[0] = m_p_select_Online[0];
		for (int j = 1; j < m_num_operators_Online; ++j) {
			m_cp_Online[j] += (m_p_select_Online[j] + m_cp_Online[j - 1]);
		}
		DCMOEA_ind<Solution<DVRP::routes, Real>> best_ind_tmp = ind;

		dominationship is_dominated;
		size_t operator_id = 0;
		m_evo_success_Online = false;


		Real p = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<Real>(0, 1);
		if (p <= m_cp_Online[0]) {
			ls_Rswap2_Online(ind, i);
			operator_id = 0;
		}
		if (p > m_cp_Online[0] && p <= m_cp_Online[1]) {
			ls_maxDT_Online(ind, i);
			operator_id = 1;
		}
		if (p > m_cp_Online[1] && p <= m_cp_Online[2]) {
			ls_2opt_Online(ind, i);
			operator_id = 2;
		}
		DVRP_CAST->construct_path(ind);
		ind.evaluate();

		m_best_inds.clear();
		m_best_inds.emplace_back(new Solution_type(ind));
		calculate_violation_objective(m_best_inds);
		mark_Solution_efeasible(m_best_inds);
		ind = *m_best_inds[0];

		is_dominated = Pareto_compare(&ind, &best_ind_tmp);
		//if (m_convergence_Online) {
		//	if (is_dominated == dominationship::Dominated || is_dominated == dominationship::Equal || is_dominated == dominationship::Non_comparable) {
		//		ind = best_ind_tmp;
		//		//m_evo_success = false;
		//	}
		//	else {
		//		m_evo_success_Online = true;
		//		m_cnt_fail_Online = 0;
		//	}
		//	if (is_dominated == dominationship::Dominating)
		//		m_score_Online[operator_id] = m_score_Online[operator_id] + 1;
		//}
		//else {
		if (is_dominated != dominationship::Dominating) {
			ind = best_ind_tmp;
			//m_evo_success = false;
		}
		else {
			m_evo_success_Online = true;
			m_cnt_fail_Online = 0;
			m_score_Online[operator_id] = m_score_Online[operator_id] + 1;
		}

		//}
		is_dominated = Pareto_compare(&ind, &best_ind_tmp);
		if (is_dominated == dominationship::Equal)
			m_cnt_fail_Online++;
		//auto b = exp(-pow(int(m_cnt_fail),0.5));
		auto decrease_level = 50;//50 * exp(-pow(int(m_cnt_fail_Online), 0.5));//0.3*int(m_individuals.size()) + 1 * 
		if (double(m_cnt_fail_Online) >= decrease_level)
			m_convergence_Online = true;

		if (m_convergence_Online) {
			m_score_Online.clear();
			m_score_Online.resize(m_num_operators_Online, 1);
			if (m_evo_success_Online) {
				//std::cout << "jump out local optimua successfully" << std::endl;
				m_evo_success_Online = false;
				m_convergence_Online = false;
				m_cnt_fail_Online = 0;
			}
		}

	}

	void MOLS_pop::ls_Rswap2_Online(Solution_type & ind, size_t i)
	{
		auto i_t = ind;
		int selected_car = i;
		int begin_pos = -1;
		for (int j = 1; j < ind.variable().m_cus_order[i].size(); ++j) {
			if (DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::depot)
				continue;
			if (DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->datum_t.is_served == false)
			{
				begin_pos = j + 1;
				break;
			}
		}
		if (begin_pos == -1 || begin_pos == 0 || begin_pos == 1 || begin_pos == ind.variable().m_cus_order[selected_car].size() - 1 || begin_pos == ind.variable().m_cus_order[selected_car].size() - 2)
			return;
		//随机选1个顾客
		int selected_pos1 = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);

		int selected_pos2 = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);
		while (selected_pos1 == selected_pos2)
			selected_pos2 = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);

		if (selected_pos1 > selected_pos2) {
			selected_pos1 = selected_pos1 + selected_pos2;
			selected_pos2 = selected_pos1 - selected_pos2;
			selected_pos1 = selected_pos1 - selected_pos2;
		}

		if (selected_pos2 < ind.variable().m_cus_order[selected_car].size() - 1) {
			int customer = ind.variable().m_cus_order[selected_car][selected_pos2];
			ind.variable().m_cus_order[selected_car][selected_pos2] = ind.variable().m_cus_order[selected_car][selected_pos1];
			ind.variable().m_cus_order[selected_car][selected_pos1] = customer;
		}
		if (ind.variable().m_cus_order[i][0] != 175) {
			std::cout << "ls_Rswap2_Online() error" << "\n";
		}
		if (ind.variable().m_cus_order[i][0] == ind.variable().m_cus_order[i][1]) {
			std::cout << "ls_Rswap2_Online() error" << "\n";
		}
	}

	void MOLS_pop::ls_maxDT_Online(Solution_type & ind, size_t i)
	{
		auto i_t = ind;
		bool flag = true;
		if (ind.objective().size() <= 2) {
			if (ind.constraint_value()[2] == 0)
				flag = false;
		}
		if (ind.objective().size() > 2) {
			if (ind.objective()[1] == 0)
				flag = false;
		}
		if (flag == true) {
			int max_delay_car = i;
			Real max_delay_time = 0.0;
			int max_delay_pos = -1;
			int begin_pos = -1;
			for (int j = 1; j < ind.variable().m_cus_order[i].size(); ++j) {
				if (DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::depot)
					continue;
				if (DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->datum_t.is_served == false)
				{
					begin_pos = j + 1;
					break;
				}
			}
			if (begin_pos == -1 || begin_pos == 0 || begin_pos == 1 || begin_pos == ind.variable().m_cus_order[i].size() - 1 || begin_pos == ind.variable().m_cus_order[i].size() - 2)
				return;
			int begin_pos1 = begin_pos;
			int begin_pos2 = begin_pos;
			//for (int j = 0; j < ind.variable().m_hi.delay_time.size(); ++j) {
			for (; begin_pos1 < ind.variable().m_hi.delay_time[i].size() - 1; ++begin_pos1) {
				if (max_delay_time < ind.variable().m_hi.delay_time[i][begin_pos1]) {
					max_delay_time = ind.variable().m_hi.delay_time[i][begin_pos1];
					max_delay_pos = begin_pos1 + 1;
				}
			}
			//}
			int max_wait_car = i;
			Real max_wait_time = 0.0;
			int max_wait_pos = -1;
			//for (int j = 0; j < offspring.variable().m_hi.wait_time.size(); ++j) {
			for (; begin_pos2 < ind.variable().m_hi.wait_time[i].size() - 1; ++begin_pos2) {
				if (max_wait_time < ind.variable().m_hi.wait_time[i][begin_pos2]) {
					max_wait_time = ind.variable().m_hi.wait_time[i][begin_pos2];
					max_wait_pos = begin_pos2 + 1;
				}
			}
			//}
			if (max_wait_time == 0)
				return;
			if (max_delay_time != 0) {
				int max_delay_cus = ind.variable().m_cus_order[max_delay_car][max_delay_pos];
				ind.variable().m_cus_order[max_delay_car].erase(ind.variable().m_cus_order[max_delay_car].begin() + max_delay_pos);
				ind.variable().m_cus_order[max_wait_car].insert(ind.variable().m_cus_order[max_wait_car].begin() + max_wait_pos, max_delay_cus);
			}
			if (ind.variable().m_cus_order[i][0] != DVRP_CAST->get_depot()) {
				std::cout << "ls_maxDT_Online() error" << "\n";
			}
			if (ind.variable().m_cus_order[i][0] == ind.variable().m_cus_order[i][1]) {
				std::cout << "ls_maxDT_Online() error" << "\n";
			}
		}


	}

	void MOLS_pop::ls_2opt_Online(Solution_type & ind, size_t i)
	{
		auto i_t = ind;
		if (ind.variable().m_cus_order[i][0] != DVRP_CAST->get_depot()) {
			std::cout << "ls_2opt_Online() error" << "\n";
		}
		int selected_car = i;
		int begin_pos = -1;
		for (int j = 1; j < ind.variable().m_cus_order[i].size(); ++j) {
			if (DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::depot)
				continue;
			if (DVRP_CAST->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->datum_t.is_served == false)
			{
				begin_pos = j + 1;
				break;
			}
		}
		if (begin_pos == -1 || begin_pos == 0 || begin_pos == 1 || begin_pos == ind.variable().m_cus_order[selected_car].size() - 1 || begin_pos == ind.variable().m_cus_order[selected_car].size() - 2)
			return;
		//随机选1个顾客
		int selected_pos1 = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);

		int selected_pos2 = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);
		while (selected_pos1 == selected_pos2)
			selected_pos2 = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);

		if (selected_pos1 > selected_pos2) {
			selected_pos1 = selected_pos1 + selected_pos2;
			selected_pos2 = selected_pos1 - selected_pos2;
			selected_pos1 = selected_pos1 - selected_pos2;
		}

		std::vector<size_t> seq_1, seq_2, seq_3;

		for (int j = 1; j < selected_pos1; ++j) {
			seq_1.push_back(ind.variable().m_cus_order[selected_car][j]);
		}
		for (int j = selected_pos2; j >= selected_pos1; --j) {
			seq_2.push_back(ind.variable().m_cus_order[selected_car][j]);
		}
		for (int j = selected_pos2 + 1; j < ind.variable().m_cus_order[selected_car].size() - 1; ++j) {
			seq_3.push_back(ind.variable().m_cus_order[selected_car][j]);
		}
		ind.variable().m_cus_order[selected_car].clear();
		ind.variable().m_cus_order[selected_car].push_back(DVRP_CAST->get_depot());

		if (seq_1.size() != 0) {
			for (auto &i : seq_1)
				ind.variable().m_cus_order[selected_car].push_back(i);
		}
		if (seq_2.size() != 0) {
			for (auto &i : seq_2)
				ind.variable().m_cus_order[selected_car].push_back(i);
		}
		if (seq_3.size() != 0) {
			for (auto &i : seq_3)
				ind.variable().m_cus_order[selected_car].push_back(i);
		}
		ind.variable().m_cus_order[selected_car].push_back(DVRP_CAST->get_depot());
		if (ind.variable().m_cus_order[i][0] != DVRP_CAST->get_depot()) {
			std::cout << "ls_2opt_Online() error" << "\n";
		}
		if (ind.variable().m_cus_order[selected_car][0] == ind.variable().m_cus_order[selected_car][1]) {
			std::cout << "ls_2opt_Online() error" << "\n";
		}
	}

	MOLS_pop::MOLS_pop(size_t size_pop) : population(size_pop) {
		for (auto& i : m_individuals) {
			i->resize_vio_obj(m_num_vio_obj);
		}
	}

	int MOLS_pop::evolve() {
		int tag = kNormalEval;

		//MOLS
		int max_eva = 1000;
		int eva = 0;
		/*while (eva <= max_eva) {*/
		Solution_type &x_t = select_sol(0);
		for (int i = 0; i < DVRP_CAST->objective_size(); ++i) {
			//Solution_type &x_t = select_sol(i);
			LS(x_t, i);
		}
		m_individuals.clear();
		for (int i = 0; i < m_archive.size(); ++i) {
			m_individuals.emplace_back(new Solution_type(*m_archive[i]));
		}
		//std::cout << "eva: " << eva << "\n";
		for (auto &i : m_individuals) {
			i->evaluate();
			//std::cout << (*i).objective()[0] << "\t"<< (*i).objective()[1] << "\t" << (*i).objective()[2] << "\t" << (*i).constraint_value()[0] << "\t" << (*i).constraint_value()[1] << "\n";
		}
		/*	eva++;
		}*/
		if (tag == kNormalEval)
			m_iteration++;
		return tag;
	}

	void MOLS_pop::LS(Solution_type & offspring, size_t obj_id)
	{
		//randomly select (N1,N2,N3)
		int nei_operator = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, 3);

		if (nei_operator == 0) {
			check_order(offspring);
			int selected_car = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, offspring.variable().m_cus_order.size());
			std::vector<Solution_type> sol_tmp;
			int selected_customer_id = 1;
			int selected_customer;
			while (offspring.variable().m_cus_order[selected_car].size() <= 3) {
				selected_car = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, offspring.variable().m_cus_order.size());
			}
			selected_customer_id = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);
			selected_customer = offspring.variable().m_cus_order[selected_car][selected_customer_id];
			//insert it into best position
			auto offspring_tmp = offspring;
			offspring.variable().m_cus_order[selected_car].erase(offspring.variable().m_cus_order[selected_car].begin() + selected_customer_id);
			auto offspring_tmp_erased = offspring;
			for (int i = 0; i < offspring.variable().m_cus_order.size(); ++i) {
				for (int j = 1; j < offspring.variable().m_cus_order[i].size() - 1; ++j) {
					/*if (i == selected_car)
						continue;*/
					offspring.variable().m_cus_order[i].insert(offspring.variable().m_cus_order[i].begin() + j, selected_customer);
					check_order(offspring);
					for (int k = 0; k < offspring.variable().m_cus_order.size(); ++k) {
						auto it = *(offspring.variable().m_cus_order[k].end() - 1);
						if (offspring.variable().m_cus_order[k].size() == 2 && offspring.variable().m_cus_order[k][0] == DVRP_CAST->get_depot() && *(offspring.variable().m_cus_order[k].end() - 1) == DVRP_CAST->get_depot()) {
							offspring.variable().m_cus_order.erase(offspring.variable().m_cus_order.begin() + k);
						}
					}
					DVRP_CAST->construct_path(offspring);
					offspring.evaluate();
					sol_tmp.push_back(offspring);
					offspring = offspring_tmp_erased;
				}
			}
			std::multimap<double, Solution_type>sol_tmp_map;

			for (int i = 0; i < sol_tmp.size(); ++i) {
				if (feasible_check(sol_tmp[i])) {//feasible check
					sol_tmp_map.insert(std::make_pair(sol_tmp[i].objective()[obj_id], sol_tmp[i]));
				}
			}
			if (!sol_tmp_map.empty()) {
				if (sol_tmp_map.begin()->first <= offspring_tmp.objective()[obj_id]) {
					offspring = sol_tmp_map.begin()->second;
				}
				else {
					updateA(sol_tmp_map.begin()->second);
					offspring = offspring_tmp;
				}
			}
			else
				offspring = offspring_tmp;
			check_order(offspring);
		}


		if (nei_operator == 1) {
			std::vector<Solution_type> sol_tmp_t;
			std::vector<Solution_type> sol_tmp;
			int selected_car = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, offspring.variable().m_cus_order.size());
			while (offspring.variable().m_cus_order[selected_car].size() <= 3) {
				selected_car = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, offspring.variable().m_cus_order.size());
			}
			int selected_customer_num = 0;
			std::vector<size_t> selected_customers;
			std::vector<size_t> selected_customers_id;
			selected_customer_num = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);
			int num_cus = 1;
			auto offspring_tmp = offspring;
			while (num_cus <= selected_customer_num) {
				int customer_id = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);
				int customer = offspring.variable().m_cus_order[selected_car][customer_id];
				selected_customers_id.push_back(customer_id);
				selected_customers.push_back(customer);
				offspring.variable().m_cus_order[selected_car].erase(offspring.variable().m_cus_order[selected_car].begin() + customer_id);
				num_cus++;
			}
			auto offspring_tmp_erased = offspring;
			for (int m = 0; m < selected_customers.size(); ++m) {
				auto offspring_tmp_erased_t = offspring;
				sol_tmp.clear();
				for (int i = 0; i < offspring.variable().m_cus_order.size(); ++i) {
					for (int j = 1; j < offspring.variable().m_cus_order[i].size() - 1; ++j) {
						offspring.variable().m_cus_order[i].insert(offspring.variable().m_cus_order[i].begin() + j, selected_customers[m]);
						//check_order(offspring);
						for (int k = 0; k < offspring.variable().m_cus_order.size(); ++k) {
							auto it = *(offspring.variable().m_cus_order[k].end() - 1);
							if (offspring.variable().m_cus_order[k].size() == 2 && offspring.variable().m_cus_order[k][0] == DVRP_CAST->get_depot() && *(offspring.variable().m_cus_order[k].end() - 1) == DVRP_CAST->get_depot()) {
								offspring.variable().m_cus_order.erase(offspring.variable().m_cus_order.begin() + k);
							}
						}
						DVRP_CAST->construct_path(offspring);
						offspring.evaluate();
						sol_tmp.push_back(offspring);
						offspring = offspring_tmp_erased_t;
					}
				}
				std::multimap<double, Solution_type>sol_tmp_map;

				for (int i = 0; i < sol_tmp.size(); ++i) {
					sol_tmp_map.insert(std::make_pair(sol_tmp[i].objective()[obj_id], sol_tmp[i]));
				}
				offspring = sol_tmp_map.begin()->second;

			}
			check_order(offspring);
			if (!feasible_check(offspring)) {//feasible check
				offspring = offspring_tmp;
			}
			if (feasible_check(offspring) && offspring.objective()[obj_id] > offspring_tmp.objective()[obj_id]) {
				updateA(offspring);
				offspring = offspring_tmp;
			}
			check_order(offspring);
		}


		if (nei_operator == 2) {
			check_order(offspring);
			int selected_car = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, offspring.variable().m_cus_order.size());
			while (offspring.variable().m_cus_order[selected_car].size() <= 3) {
				selected_car = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(0, offspring.variable().m_cus_order.size());
			}
			std::vector<Solution_type> sol_tmp;
			int selected_customer_id = 1;
			std::vector<size_t> selected_customers1;
			std::vector<size_t> selected_customers1_id;
			selected_customer_id = global::ms_global->m_uniform[caller::Algorithm]->next_non_standard<int>(1, offspring.variable().m_cus_order[selected_car].size() - 2);
			//insert it into best position
			auto offspring_tmp = offspring;
			for (int i = selected_customer_id + 1; i < offspring_tmp.variable().m_cus_order[selected_car].size() - 1; ++i) {
				if (offspring_tmp.variable().m_cus_order[selected_car][i] == DVRP_CAST->get_depot())
					break;
				selected_customers1_id.push_back(i);
				selected_customers1.push_back(offspring_tmp.variable().m_cus_order[selected_car][i]);
				auto it = std::find(offspring.variable().m_cus_order[selected_car].begin(), offspring.variable().m_cus_order[selected_car].end(), offspring_tmp.variable().m_cus_order[selected_car][i]);
				if (it != offspring.variable().m_cus_order[selected_car].end())
					offspring.variable().m_cus_order[selected_car].erase(it);

			}
			auto offspring_tmp_erased = offspring;
			for (int i = 0; i < offspring.variable().m_cus_order.size(); ++i) {
				for (int j = 1; j < offspring.variable().m_cus_order[i].size() - 2; ++j) {
					if (i == selected_car)
						continue;
					auto offspring_tmp2 = offspring;
					std::vector<size_t> selected_customers2;
					std::vector<size_t> selected_customers2_id;
					for (int k = j + 1; k < offspring_tmp2.variable().m_cus_order[i].size() - 1; ++k) {
						if (offspring_tmp2.variable().m_cus_order[i][k] == DVRP_CAST->get_depot())
							break;
						selected_customers2_id.push_back(k);
						selected_customers2.push_back(offspring_tmp2.variable().m_cus_order[i][k]);
						auto it = std::find(offspring.variable().m_cus_order[i].begin(), offspring.variable().m_cus_order[i].end(), offspring_tmp2.variable().m_cus_order[i][k]);
						if (it != offspring.variable().m_cus_order[i].end())
							offspring.variable().m_cus_order[i].erase(it);
					}
					auto offspring_tmp2_erased = offspring;
					offspring.variable().m_cus_order[i].insert(offspring.variable().m_cus_order[i].begin() + j + 1, selected_customers1.begin(), selected_customers1.end());
					offspring.variable().m_cus_order[selected_car].insert(offspring.variable().m_cus_order[selected_car].begin() + selected_customer_id + 1, selected_customers2.begin(), selected_customers2.end());
					check_order(offspring);
					for (int k = 0; k < offspring.variable().m_cus_order.size(); ++k) {
						if (offspring.variable().m_cus_order[k].size() == 2 && offspring.variable().m_cus_order[k][0] == DVRP_CAST->get_depot() && *(offspring.variable().m_cus_order[k].end() - 1) == DVRP_CAST->get_depot()) {
							offspring.variable().m_cus_order.erase(offspring.variable().m_cus_order.begin() + k);
						}
					}
					check_order(offspring);
					DVRP_CAST->construct_path(offspring);
					offspring.evaluate();
					sol_tmp.push_back(offspring);
					offspring = offspring_tmp_erased;
				}
			}
			std::multimap<double, Solution_type>sol_tmp_map;
			for (int i = 0; i < sol_tmp.size(); ++i) {
				if (feasible_check(sol_tmp[i])) {//feasible check
					sol_tmp_map.insert(std::make_pair(sol_tmp[i].objective()[obj_id], sol_tmp[i]));
				}
			}
			if (!sol_tmp_map.empty()) {
				if (sol_tmp_map.begin()->first <= offspring_tmp.objective()[obj_id]) {
					offspring = sol_tmp_map.begin()->second;
				}
				else {
					updateA(sol_tmp_map.begin()->second);
					offspring = offspring_tmp;
				}
			}
			else
				offspring = offspring_tmp;
			check_order(offspring);
		}
	}

	void MOLS_pop::updateA(Solution_type offspring)
	{
		if (m_archive.size() >= m_pop_size)
			return;
		m_archive.emplace_back(new Solution_type(offspring));
	}

	bool MOLS_pop::feasible_check(Solution_type & offspring)
	{
		for (int i = 0; i < offspring.constraint_value().size(); ++i) {
			if (offspring.constraint_value()[i] > 0)
				return false;
		}
		return true;
	}

	bool MOLS_pop::check_order(Solution_type & offspring)
	{
		int num_order = 0;
		for (auto &i : offspring.variable().m_cus_order) {
			num_order += (int(i.size()) - 2);
		}
		if (num_order == DVRP_CAST->get_current_cus_order().size())
			return true;
		else if (num_order < DVRP_CAST->get_current_cus_order().size()) {
			std::cout << "order is not enough" << "\n";
			return false;
		}
		else {
			std::cout << "order is repeat" << "\n";
			return false;
		}
	}




	dominationship MOLS_pop::e_Pareto_compare(Solution_type* const&s1, Solution_type* const&s2) {
		/* One efeasible one in-efeasible */
		if (s1->get_efeasible() != s2->get_efeasible()) {
			if (s1->get_efeasible())
				return dominationship::Dominating;
			else
				return dominationship::Dominated;
		}

		/* Both efeasible */
		else if (s1->get_efeasible() && s2->get_efeasible()) {
			auto nor_obj_result = objective_compare<Real>(s1->objective(), s2->objective(), DVRP_CAST->opt_mode());
			auto vio_obj_result = objective_compare<Real>(s1->get_vio_obj(), s2->get_vio_obj(), optimization_mode::Minimization);

			if (nor_obj_result == dominationship::Dominating&&vio_obj_result == dominationship::Equal)
				return dominationship::Dominating;
			if (nor_obj_result == dominationship::Equal&&vio_obj_result == dominationship::Dominating)
				return dominationship::Dominating;
			if (nor_obj_result == dominationship::Dominating&&vio_obj_result == dominationship::Dominating)
				return dominationship::Dominating;

			if (nor_obj_result == dominationship::Dominated&&vio_obj_result == dominationship::Equal)
				return dominationship::Dominated;
			if (nor_obj_result == dominationship::Equal&&vio_obj_result == dominationship::Dominated)
				return dominationship::Dominated;
			if (nor_obj_result == dominationship::Dominated&&vio_obj_result == dominationship::Dominated)
				return dominationship::Dominated;

			if (nor_obj_result == dominationship::Dominated&&vio_obj_result == dominationship::Dominating)
				return dominationship::Non_dominated;
			if (nor_obj_result == dominationship::Dominating&&vio_obj_result == dominationship::Dominated)
				return dominationship::Non_dominated;
			if (nor_obj_result == dominationship::Non_dominated || vio_obj_result == dominationship::Non_dominated)
				return dominationship::Non_dominated;

			if (nor_obj_result == dominationship::Equal&&vio_obj_result == dominationship::Equal)
				return dominationship::Equal;
		}
		//return objective_compare<Real>(s1->objective(), s2->objective(), CONTINUOUS_CAST->opt_mode());

	/* Both in-efeasible */
		else {
			if (s1->get_vio_obj()[0] < s2->get_vio_obj()[0])
				return dominationship::Dominating;
			else if (s1->get_vio_obj()[0] > s2->get_vio_obj()[0])
				return dominationship::Dominated;
			else
				return dominationship::Equal;
		}
	}

	dominationship MOLS_pop::Pareto_compare(Solution_type * const & s1, Solution_type * const & s2)
	{
		auto nor_obj_result = objective_compare<Real>(s1->objective(), s2->objective(), DVRP_CAST->opt_mode());
		auto vio_obj_result = objective_compare<Real>(s1->get_vio_obj(), s2->get_vio_obj(), optimization_mode::Minimization);

		if (nor_obj_result == dominationship::Dominating&&vio_obj_result == dominationship::Equal)
			return dominationship::Dominating;
		if (nor_obj_result == dominationship::Equal&&vio_obj_result == dominationship::Dominating)
			return dominationship::Dominating;
		if (nor_obj_result == dominationship::Dominating&&vio_obj_result == dominationship::Dominating)
			return dominationship::Dominating;

		if (nor_obj_result == dominationship::Dominated&&vio_obj_result == dominationship::Equal)
			return dominationship::Dominated;
		if (nor_obj_result == dominationship::Equal&&vio_obj_result == dominationship::Dominated)
			return dominationship::Dominated;
		if (nor_obj_result == dominationship::Dominated&&vio_obj_result == dominationship::Dominated)
			return dominationship::Dominated;

		if (nor_obj_result == dominationship::Dominated&&vio_obj_result == dominationship::Dominating)
			return dominationship::Non_dominated;
		if (nor_obj_result == dominationship::Dominating&&vio_obj_result == dominationship::Dominated)
			return dominationship::Non_dominated;
		if (nor_obj_result == dominationship::Non_dominated || vio_obj_result == dominationship::Non_dominated)
			return dominationship::Non_dominated;

		if (nor_obj_result == dominationship::Equal&&vio_obj_result == dominationship::Equal)
			return dominationship::Equal;
	}

	MOLS::MOLS(param_map & v) : algorithm(v.at("algorithm name")), m_pop(v.at("population size")) {
		DVRP_CAST->set_eval_monitor_flag(true);
	}

	void MOLS::initialize() {
		m_pop.initialize();
#ifdef OFEC_DEMO
		//std::cout << "pop[i]" << "\t" << "rank\troute_length\twait_time\tdelay_time\tover_load\tover_time\tnum_car" << std::endl;
		vector<vector<Solution<DVRP::routes, Real>*>> pops(1);
		for (size_t i = 0; i < m_pop.size(); ++i) {
			pops[0].emplace_back(&m_pop[i]);
			//rank_.push_back(m_pop[i].fitness());
		}
		dynamic_cast<ofec_demo::buffer_DVRP*>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops, 0);
#endif
	}

	void MOLS::run_() {
		DCMOEA_ind<Solution<DVRP::routes, Real>> best_ind = m_pop[0];
		DCMOEA_ind<Solution<DVRP::routes, Real>> ini_best_ind = m_pop[0];
		m_pop.sort();
		Real t_route_length = 0;
		Real t_wait_time = 0;
		Real t_delay_time = 0;
		Real t_over_load = 0;
		Real t_over_time = 0;
		std::string initialPOP_Fname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/InitialPOP" + to_string(m_pop.size()) + "_order" + to_string(DVRP_CAST->get_current_cus_order().size()) + "MOLS_MOVRPRTC_AHC_initial.txt";
		std::ofstream pop_(initialPOP_Fname);
		std::stringstream pop_str;
		for (size_t i = 0; i < m_pop.size(); i++) {
			t_route_length = m_pop[i].objective()[0];
			t_wait_time = m_pop[i].objective()[2];
			t_delay_time = m_pop[i].objective()[1];
			t_over_load = m_pop[i].constraint_value()[0];
			t_over_time = m_pop[i].constraint_value()[1];

			int t_num_car = m_pop[i].variable().m_cus_order.size();
			pop_str << to_string(i) << "\t" << to_string(m_pop[i].fitness()) << "\t" << t_route_length << "\t" << t_wait_time << "\t" << t_delay_time << "\t" << t_over_load << "\t" << t_over_time << "\t" << t_num_car << std::endl;
		}
		pop_ << pop_str.str();
		pop_.close();

		int iter = 0;
		while (!terminating()) {
			//auto domi_flag = m_pop.Pareto_compare(&best_ind, &ini_best_ind);
			if (iter >= 10000) {//10000  //&& domi_flag==dominationship::Dominating
				std::cout << "Online!" << std::endl;
				m_pop.Online();
				std::string Online_resFname;
				Online_resFname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/Online_result" + to_string(m_pop.size()) + "_order" + to_string(DVRP_CAST->get_current_cus_order().size()) + "MOLS_MOVRPRTC_AHC_initial.txt";
				std::ofstream Online_result(Online_resFname);
				std::stringstream Online_result_;
				if (DVRP_CAST->objective_size() == 3) {
					Online_result_ << m_pop.get_bestInd().objective()[0] << "\t" << m_pop.get_bestInd().objective()[2] << "\t" << m_pop.get_bestInd().objective()[1] << "\t" << m_pop.get_bestInd().constraint_value()[0] << "\t" << m_pop.get_bestInd().constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << m_pop.get_bestInd().variable().m_cus_order.size() << std::endl;
					std::cout << m_pop.get_bestInd().objective()[0] << "\t" << m_pop.get_bestInd().objective()[2] << "\t" << m_pop.get_bestInd().objective()[1] << "\t" << m_pop.get_bestInd().constraint_value()[0] << "\t" << m_pop.get_bestInd().constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << m_pop.get_bestInd().variable().m_cus_order.size() << std::endl;
				}
				if (DVRP_CAST->objective_size() == 1) {
					Online_result_ << m_pop.get_bestInd().objective()[0] << "\t" << m_pop.get_bestInd().constraint_value()[3] << "\t" << m_pop.get_bestInd().constraint_value()[2] << "\t" << m_pop.get_bestInd().constraint_value()[0] << "\t" << m_pop.get_bestInd().constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << m_pop.get_bestInd().variable().m_cus_order.size() << std::endl;
					std::cout << m_pop.get_bestInd().objective()[0] << "\t" << m_pop.get_bestInd().constraint_value()[3] << "\t" << m_pop.get_bestInd().constraint_value()[2] << "\t" << m_pop.get_bestInd().constraint_value()[0] << "\t" << m_pop.get_bestInd().constraint_value()[1] << "\t" << std::to_string(DVRP_CAST->get_future_cus_order().size()) << "\t" << m_pop.get_bestInd().variable().m_cus_order.size() << std::endl;
				}
				Online_result << Online_result_.str();
				Online_result.close();
				std::cout << "Online Complete!" << "\n";
				system("pause");
				//m_isOnline = true;
			}
			m_pop.evolve();
			m_pop.sort();
#ifdef OFEC_DEMO
			/*vector<vector<Solution<DVRP::routes, Real>*>> pops(1);
			for (size_t i = 0; i < m_pop.size(); ++i)
				pops[0].emplace_back(&m_pop[i]);
			dynamic_cast<ofec_demo::buffer_DVRP*>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops);*/

			/* Begin terminal window output */
			size_t evals = DVRP_CAST->evaluations();
			Real min_obj = (std::numeric_limits<Real>::max)();
			size_t idx_best;
			Real temp_obj;


			std::vector<std::pair<double, double>> objective_min_max(DVRP_CAST->objective_size());
			std::vector<std::pair<double, double>> constraint_min_max(DVRP_CAST->num_constraints());
			std::vector<std::pair<double, double>> vioObj_min_max(m_pop.get_num_vio_obj());

			for (int i = 0; i < DVRP_CAST->objective_size(); ++i) {
				std::pair<double, double> min_max;
				min_max.first = (std::numeric_limits<Real>::max)();
				min_max.second = (std::numeric_limits<Real>::min)();
				for (size_t j = 0; j < m_pop.size(); j++) {
					if (min_max.first >= m_pop[j].objective()[i])
						min_max.first = m_pop[j].objective()[i];
					if (min_max.second <= m_pop[j].objective()[i])
						min_max.second = m_pop[j].objective()[i];
				}
				objective_min_max[i] = min_max;
			}
			for (int i = 0; i < DVRP_CAST->num_constraints(); ++i) {
				std::pair<double, double> min_max;
				min_max.first = (std::numeric_limits<Real>::max)();
				min_max.second = (std::numeric_limits<Real>::min)();
				for (size_t j = 0; j < m_pop.size(); j++) {
					if (min_max.first >= m_pop[j].constraint_value()[i])
						min_max.first = m_pop[j].constraint_value()[i];
					if (min_max.second <= m_pop[j].constraint_value()[i])
						min_max.second = m_pop[j].constraint_value()[i];
				}
				constraint_min_max[i] = min_max;
			}
			for (int i = 0; i < m_pop.get_num_vio_obj(); ++i) {
				std::pair<double, double> min_max;
				min_max.first = (std::numeric_limits<Real>::max)();
				min_max.second = (std::numeric_limits<Real>::min)();
				for (size_t j = 0; j < m_pop.size(); j++) {
					if (min_max.first >= m_pop[j].get_vio_obj()[i])
						min_max.first = m_pop[j].get_vio_obj()[i];
					if (min_max.second <= m_pop[j].get_vio_obj()[i])
						min_max.second = m_pop[j].get_vio_obj()[i];
				}
				vioObj_min_max[i] = min_max;
			}

			Real route_length = 0;
			Real wait_time = 0;
			Real delay_time = 0;
			Real over_load = 0;
			Real over_time = 0;
			std::string result_Fname;
			result_Fname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/result_pop" + to_string(m_pop.pop_size()) + "_order" + to_string(DVRP_CAST->get_current_cus_order().size()) + "MOLS_MOVRPRTC_AHC_initial.txt";
			std::ofstream result(result_Fname, ios::app);
			std::stringstream result_;
			result_ << iter << "\n";

			for (size_t i = 0; i < m_pop.size(); i++) {
				route_length = m_pop[i].objective()[0];
				wait_time = m_pop[i].objective()[2];
				delay_time = m_pop[i].objective()[1];
				over_load = m_pop[i].constraint_value()[0];
				over_time = m_pop[i].constraint_value()[1];

				int num_car = m_pop[i].variable().m_cus_order.size();
				//std::cout << "pop[" + to_string(i) << "]\t" << to_string(m_pop[i].fitness()) << "\t" << route_length << "\t" << wait_time << "\t" << delay_time << "\t" << over_load << "\t" << over_time << "\t" << num_car << std::endl;
				result_ << to_string(i) << "\t" << to_string(m_pop[i].fitness()) << "\t" << route_length << "\t" << wait_time << "\t" << delay_time << "\t" << over_load << "\t" << over_time << "\t" << num_car << std::endl;



				if (m_pop[i].fitness() == 0) {
					temp_obj = 0;
					for (int j = 0; j < DVRP_CAST->objective_size(); ++j) {
						if (objective_min_max[j].second != 0 && objective_min_max[j].first != objective_min_max[j].second) {
							temp_obj += ((m_pop[i].objective()[j] - objective_min_max[j].first) / (objective_min_max[j].second - objective_min_max[j].first));
						}
						else
						{
							temp_obj += m_pop[i].objective()[j];
						}
					}
					for (int j = 0; j < DVRP_CAST->num_constraints(); ++j) {
						if (constraint_min_max[j].second != 0 && constraint_min_max[j].first != constraint_min_max[j].second)
							temp_obj += ((m_pop[i].constraint_value()[j] - constraint_min_max[j].first) / (constraint_min_max[j].second - constraint_min_max[j].first));
						else
							temp_obj += m_pop[i].constraint_value()[j];
					}
					for (int j = 0; j < m_pop.get_num_vio_obj(); ++j) {
						if (vioObj_min_max[j].second != 0 && vioObj_min_max[j].first != vioObj_min_max[j].second)
							temp_obj += ((m_pop[i].get_vio_obj()[j] - vioObj_min_max[j].first) / (vioObj_min_max[j].second - vioObj_min_max[j].first));
						else
							temp_obj += m_pop[i].get_vio_obj()[j];
					}
					if (min_obj >= temp_obj) {
						min_obj = temp_obj;
						idx_best = i;
					}
				}
			}
			result << result_.str();
			result.close();
			Real _route_length = 0;
			Real _wait_time = 0;
			Real _delay_time = 0;
			Real _over_load = 0;
			Real _over_time = 0;
			_route_length = m_pop[idx_best].objective()[0];
			_wait_time = m_pop[idx_best].objective()[2];
			_delay_time = m_pop[idx_best].objective()[1];
			_over_load = m_pop[idx_best].constraint_value()[0];
			_over_time = m_pop[idx_best].constraint_value()[1];

			std::cout << iter << "\t" << _route_length << "\t" << _wait_time << "\t" << _delay_time << "\t" << _over_load << "\t" << _over_time << std::endl;

			std::string progress_Fname;
			progress_Fname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/progress" + to_string(m_pop.pop_size()) + "_order" + to_string(DVRP_CAST->get_current_cus_order().size()) + "MOLS_MOVRPRTC_AHC_initial.txt";
			std::ofstream progress(progress_Fname, ios::app);
			std::stringstream progress_;
			progress_ << iter << "\t" << _route_length << "\t" << _wait_time << "\t" << _delay_time << "\t" << _over_load << "\t" << _over_time << std::endl;
			progress << progress_.str();
			progress.close();

			best_ind = m_pop[idx_best];
			if (iter == 0)
				ini_best_ind = best_ind;
			m_pop.set_bestInd(best_ind);

			vector<vector<Solution<DVRP::routes, Real>*>> pops(1);
			for (size_t i = 0; i < m_pop.size(); ++i)
				pops[0].emplace_back(&m_pop[i]);
			dynamic_cast<ofec_demo::buffer_DVRP*>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops, idx_best);
#endif
			iter++;
		}
		//m_pop.sort();
	}

	void MOLS::record() {
		m_pop.sort();
		size_t evals = DVRP_CAST->evaluations();
		Real min_obj = (std::numeric_limits<Real>::max)();
		size_t idx_best;
		Real temp_obj;
		for (size_t i = 0; i < m_pop.size(); i++) {
			if (m_pop[i].fitness() == 0) {
				temp_obj = 0;
				for (int j = 0; j < DVRP_CAST->objective_size(); ++j) {
					temp_obj += m_pop[i].objective()[j];
				}
				for (int j = 0; j < DVRP_CAST->num_constraints(); ++j) {
					temp_obj += m_pop[i].constraint_value()[j];
				}
				for (int j = 0; j < m_pop.get_num_vio_obj(); ++j) {
					temp_obj += m_pop[i].get_vio_obj()[j];
				}
				if (temp_obj < min_obj) {
					min_obj = temp_obj;
					idx_best = i;
				}
			}
		}
		Real route_length = 0;
		Real wait_time = 0;
		Real delay_time = 0;
		Real over_load = 0;
		Real over_time = 0;

		if (DVRP_CAST->objective_size() == 2) {
			route_length = m_pop[idx_best].objective()[0];
			wait_time = m_pop[idx_best].objective()[1];
			delay_time = m_pop[idx_best].constraint_value()[2];
			over_load = m_pop[idx_best].constraint_value()[0];
			over_time = m_pop[idx_best].constraint_value()[1];
		}
		if (DVRP_CAST->objective_size() == 1) {
			route_length = m_pop[idx_best].objective()[0];
			wait_time = m_pop[idx_best].constraint_value()[3];
			delay_time = m_pop[idx_best].constraint_value()[2];
			over_load = m_pop[idx_best].constraint_value()[0];
			over_time = m_pop[idx_best].constraint_value()[1];
		}
		measure::get_measure()->record(global::ms_global.get(), evals, route_length, wait_time, delay_time, over_load, over_time);
	}
}