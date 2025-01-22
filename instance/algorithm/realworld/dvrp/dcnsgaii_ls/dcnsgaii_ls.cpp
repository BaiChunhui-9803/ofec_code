#include "DCNSGAII_LS.h"
//#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../record/rcr_vec_real.h"
#include <numeric>

using namespace std::placeholders;

#ifdef OFEC_DEMO
#include "../../../../../../ui/buffer/realworld/DVRP/buffer_DVRP.h"
extern unique_ptr<ofec_demo::scene> ofec_demo::msp_buffer;
#endif

namespace ofec {
	void DCNSGAII_LS_pop::initialize(Problem *pro, Random *rnd) {
		//test_lkh();
		m_ahc_tsp.ahc_tsp_init_pop(m_individuals, pro, rnd);
		//m_ahc_tsp.randomly_init_pop(m_individuals);
		//evaluate();
		///****************conference paper experiment********************************/
		//记录订单，对比订单是否相同
		std::string current_order_Fname;
		current_order_Fname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/current_order" + to_string(DVRP_CAST(pro)->get_current_cus_order().size()) + "_" + m_algName_ + "_" + m_proName_ + "_" + m_initialName_ + ".txt";
		std::ofstream Corder_(current_order_Fname);
		std::stringstream Corder_str;
		for (auto &i : DVRP_CAST(pro)->get_current_cus_order()) {
			Corder_str << i << " ";
		}
		Corder_ << Corder_str.str();
		Corder_.close();
		///**************************************************************************/
	}

	void DCNSGAII_LS_pop::initializeAfterEvaluation() {
		calculate_initial_max_violation(m_individuals);
		m_e = m_max_G;
		calculate_violation_objective(m_individuals);
		mark_Solution_efeasible(m_individuals);
		m_score.clear();
		m_score.resize(m_num_operators, 1);
		m_p_select.clear();
		m_p_select.resize(m_num_operators);
		m_cp.clear();
		m_cp.resize(m_num_operators);
	}

	void DCNSGAII_LS_pop::sort() {
		std::vector<int> ranks_;
		std::vector<IndividualType *> pop_;
		std::function<Dominance(IndividualType *const &, IndividualType *const &)> comp = std::bind(&DCNSGAII_LS_pop::Pareto_compare, this, _1, _2);
		for (auto &i : m_individuals)
			pop_.emplace_back(i.get());
		nd_sort::fast_sort<IndividualType *>(pop_, ranks_, comp);
		for (size_t i = 0; i < m_individuals.size(); i++) {
			m_individuals[i]->setFitness(ranks_[i]);
		}
	}

	void DCNSGAII_LS_pop::Online(Problem *pro, Algorithm *alg, Random *rnd)
	{
		if (pro->numberObjectives() == 3)
			std::cout << "Offline WT: " << m_best_ind_Offline.objective()[2] << "\n";
		if (pro->numberObjectives() == 1)
			std::cout << "Offline WT: " << m_best_ind_Offline.constraint()[3] << "\n";
		DVRP_CAST(pro)->set_Online(true);

		m_best_ind_Online = m_best_ind_Offline;

		//以一个最优方案为基础，处理新的订单。最优订单的定义为个指标归一化之后

		for (int i = 0; i < m_best_ind_Online.variable().m_cus_order.size(); ++i) {
			double presentTime = DVRP_CAST(pro)->get_depart_time() + DVRP_CAST(pro)->get_cost_time(DVRP_CAST(pro)->get_depart_time(), m_best_ind_Online.variable().m_cus_order[i][0], m_best_ind_Online.variable().m_cus_order[i][1]);
			for (int j = 1; j < m_best_ind_Online.variable().m_cus_order[i].size() - 1; ++j) {
				if (DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][j]]->datum_t.is_served == true)
					continue;
				DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][j]]->datum_t.is_served = true;
				//updateSeqOnline(m_best_ind_Online, i);
				auto customer = m_best_ind_Online.variable().m_cus_order[i][j];
				auto presentTime_curOrder = presentTime;
				auto presentTime_newOrder = presentTime;
				std::pair<int, int> new_order;
				m_insert_success = false;
				for (int k = 0; k < DVRP_CAST(pro)->get_future_cus_order().size(); ++k) {
					if (presentTime >= DVRP_CAST(pro)->get_future_cus_order()[k].second) {
						new_order = DVRP_CAST(pro)->get_future_cus_order()[k];
						if (DVRP_CAST(pro)->get_net().get_road_net()[new_order.first]->type != DVRP::node_type::future_customer || new_order.first == 175) {
							std::cout << "error" << "\n";
						}
						size_t next_cur_order = m_best_ind_Online.variable().m_cus_order[i][j + 1];
						//DVRP_CAST(pro)->set_Online(false);
						presentTime = presentTime + DVRP_CAST(pro)->get_net().get_road_net()[customer]->server_time + DVRP_CAST(pro)->get_cost_time(presentTime, customer, new_order.first);
						if (presentTime >= DVRP_CAST(pro)->get_net().get_road_net()[new_order.first]->datum_t.begin_time && presentTime <= DVRP_CAST(pro)->get_net().get_road_net()[new_order.first]->datum_t.end_time) {

							m_best_ind_Online.variable().m_cus_order[i].insert(m_best_ind_Online.variable().m_cus_order[i].begin() + j + 1, new_order.first);
							DVRP_CAST(pro)->get_future_cus_order().erase(DVRP_CAST(pro)->get_future_cus_order().begin() + k);
							m_insert_success = true;
							if (m_best_ind_Online.variable().m_cus_order[i][0] != DVRP_CAST(pro)->get_depot()) {
								std::cout << "error" << "\n";
							}
							updateSeqOnline_inserted(m_best_ind_Online, i, new_order, pro, alg, rnd);

							presentTime_newOrder = presentTime;
							presentTime = presentTime + DVRP_CAST(pro)->get_net().get_road_net()[new_order.first]->server_time;
							presentTime = presentTime + DVRP_CAST(pro)->get_cost_time(presentTime, new_order.first, next_cur_order);

							/*double delay_time = 0;
							if (presentTime > DVRP_CAST(pro)->get_net().get_road_net()[next_cur_order]->datum_t.end_time) {
								delay_time = presentTime - DVRP_CAST(pro)->get_net().get_road_net()[next_cur_order]->datum_t.end_time;

								auto it= std::find(m_best_ind_Offline.variable().m_cus_order[i].begin(), m_best_ind_Offline.variable().m_cus_order[i].end(), next_cur_order);
								auto nextOrder_pos = std::distance(std::begin(m_best_ind_Offline.variable().m_cus_order[i]), it);
								if (delay_time > m_best_ind_Offline.variable().m_hi.delay_time[i][nextOrder_pos - 1]) {
									m_best_ind_Online.variable().m_cus_order[i].erase(m_best_ind_Online.variable().m_cus_order[i].begin() + j + 1);
									DVRP_CAST(pro)->get_future_cus_order().push_back(new_order);
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
										DVRP_CAST(pro)->get_future_cus_order().push_back(new_order);
										m_insert_success = false;
									}
								}
							}*/
							if (m_insert_success == true) {
								bool is_over_load = false;
								double deliveryWeight = 0;
								double over_load = 0;
								for (size_t m = 0; m < m_best_ind_Online.variable().m_cus_order[i].size(); ++m) {
									if (DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.is_served == false)
										deliveryWeight += DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.delivery_weight;
								}
								if (deliveryWeight > DVRP_CAST(pro)->get_vehicle_property()[2].capacity) {
									over_load += deliveryWeight - DVRP_CAST(pro)->get_vehicle_property()[2].capacity;
								}
								for (size_t m = 0; m < m_best_ind_Online.variable().m_cus_order[i].size(); ++m) {
									if (DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.is_served == false)
										deliveryWeight = deliveryWeight - DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.delivery_weight + DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.pick_weight;
									else
										deliveryWeight = deliveryWeight - 0 + DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.pick_weight;

									if (deliveryWeight > DVRP_CAST(pro)->get_vehicle_property()[2].capacity) {
										over_load += deliveryWeight - DVRP_CAST(pro)->get_vehicle_property()[2].capacity;
									}
								}
								if (over_load > 0) {
									m_best_ind_Online.variable().m_cus_order[i].erase(m_best_ind_Online.variable().m_cus_order[i].begin() + j + 1);
									DVRP_CAST(pro)->get_future_cus_order().push_back(new_order);
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
						//std::cout << m_best_ind_Online.objective()[0] << "\t" << m_best_ind_Online.constraint()[3] << "\t" << m_best_ind_Online.constraint()[2] << "\t" << m_best_ind_Online.constraint()[0] << "\t" << m_best_ind_Online.constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << "\t" << m_best_ind_Online.variable().m_cus_order.size() << std::endl;

					}

				}
				if (m_insert_success == false) {
					presentTime = presentTime_curOrder;
					presentTime = presentTime + DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][j]]->server_time + DVRP_CAST(pro)->get_cost_time(presentTime, m_best_ind_Online.variable().m_cus_order[i][j], m_best_ind_Online.variable().m_cus_order[i][j + 1]);
				}
			}
			//标记已被服务
			for (int m = 1; m < m_best_ind_Online.variable().m_cus_order[i].size() - 1; ++m) {
				DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.is_served = false;
			}
			for (int n = 1; n < m_best_ind_Online.variable().m_cus_order[i].size() - 1; ++n) {
				if (DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][n]]->type == DVRP::node_type::future_customer)
					break;
				DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][n]]->datum_t.is_served = true;
			}
			DVRP_CAST(pro)->set_Online(true);
			construct_path_Online(m_best_ind_Online, i, pro);
		}
		m_best_ind_Online.evaluate(pro, alg);
		if (DVRP_CAST(pro)->numberObjectives() == 1)
			std::cout << "Online WT: " << m_best_ind_Online.constraint()[3] << "\n";
		if (DVRP_CAST(pro)->numberObjectives() == 3)
			std::cout << "Online WT: " << m_best_ind_Online.objective()[2] << "\n";
	}

	void DCNSGAII_LS_pop::construct_path_Online(IndividualType &ind, size_t i, Problem *pro)
	{
		ind.variable().m_path[i].clear();
		ind.variable().m_path_node_type[i].clear();

		int cnt = 0;
		std::vector<std::vector<size_t>> path_temp(ind.variable().m_cus_order[i].size() - 1);
		std::vector<std::vector<DVRP::node_type>> path_node_type_temp(ind.variable().m_cus_order[i].size() - 1);
		//path_temp[0].push_back(m_depot);
		Real present_time = DVRP_CAST(pro)->get_depart_time();
		for (size_t j = 0; j < ind.variable().m_cus_order[i].size() - 1; j++) {
			std::vector<size_t> Shortest_path;

			DVRP_CAST(pro)->shortest_path(ind.variable().m_cus_order[i][j], ind.variable().m_cus_order[i][j + 1], Shortest_path, present_time);//时间要变


			if (DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::current_customer || DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::future_customer) {
				present_time = present_time + DVRP_CAST(pro)->get_cost_time(present_time, ind.variable().m_cus_order[i][j], ind.variable().m_cus_order[i][j + 1]) + DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->server_time;
			}
			if (DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::depot) {
				present_time = present_time + DVRP_CAST(pro)->get_cost_time(present_time, ind.variable().m_cus_order[i][j], ind.variable().m_cus_order[i][j + 1]);
			}
			auto tt = DVRP_CAST(pro)->get_cost_time(present_time, ind.variable().m_cus_order[i][j], ind.variable().m_cus_order[i][j + 1]);
			Shortest_path.pop_back();
			path_temp[j] = Shortest_path;
			for (auto &m : path_temp[j]) {
				ind.variable().m_path[i].push_back(m);
				path_node_type_temp[j].push_back(DVRP::node_type::general_node);
			}
			path_node_type_temp[j][0] = DVRP_CAST(pro)->get_net().get_road_net()[path_temp[j][0]]->type;
			for (auto &m : path_node_type_temp[j]) {
				ind.variable().m_path_node_type[i].push_back(m);
			}
		}
		ind.variable().m_path[i].push_back(DVRP_CAST(pro)->get_depot());
		ind.variable().m_path_node_type[i].push_back(DVRP_CAST(pro)->get_net().get_road_net()[DVRP_CAST(pro)->get_depot()]->type);
	}

	void DCNSGAII_LS_pop::updateSeqOnline(IndividualType &ind, size_t i, Problem *pro, Algorithm *alg, Random *rnd)
	{
		DVRP_CAST(pro)->set_Online(false);

		DVRP_CAST(pro)->construct_path(ind);
		ind.evaluate(pro, alg);

		m_score_Online.clear();
		m_score_Online.resize(m_num_operators_Online, 1);
		m_p_select_Online.clear();
		m_p_select_Online.resize(m_num_operators_Online);
		m_cp_Online.clear();
		m_cp_Online.resize(m_num_operators_Online);

		m_best_inds.clear();
		m_best_inds.emplace_back(new IndividualType(ind));
		calculate_initial_max_violation(m_best_inds);
		m_e = m_max_G;
		m_k = 0;
		std::cout << "updateSeqOnline() before: " << "\n";
		if (DVRP_CAST(pro)->numberObjectives() == 1)
			std::cout << ind.objective()[0] << "\t" << ind.constraint()[3] << "\t" << ind.constraint()[2] << "\t" << ind.constraint()[0] << "\t" << ind.constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		if (DVRP_CAST(pro)->numberObjectives() == 3)
			std::cout << ind.objective()[0] << "\t" << ind.objective()[2] << "\t" << ind.objective()[1] << "\t" << ind.constraint()[0] << "\t" << ind.constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;

		int iter_cnt = 0;
		while (iter_cnt < m_OnlineIters) {
			Online_ALS(ind, i, pro, alg, rnd);
			iter_cnt++;
		}
		std::cout << "updateSeqOnline() end: " << "\n";
		if (DVRP_CAST(pro)->numberObjectives() == 1)
			std::cout << ind.objective()[0] << "\t" << ind.constraint()[3] << "\t" << ind.constraint()[2] << "\t" << ind.constraint()[0] << "\t" << ind.constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		if (DVRP_CAST(pro)->numberObjectives() == 3)
			std::cout << ind.objective()[0] << "\t" << ind.objective()[2] << "\t" << ind.objective()[1] << "\t" << ind.constraint()[0] << "\t" << ind.constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;

		DVRP_CAST(pro)->set_Online(true);
	}

	void DCNSGAII_LS_pop::updateSeqOnline_inserted(IndividualType &ind, size_t i, std::pair<int, int> new_order, Problem *pro, Algorithm *alg, Random *rnd)
	{
		DVRP_CAST(pro)->set_Online(false);
		DVRP_CAST(pro)->construct_path(ind);
		ind.evaluate(pro, alg);

		m_score_Online.clear();
		m_score_Online.resize(m_num_operators_Online, 1);
		m_p_select_Online.clear();
		m_p_select_Online.resize(m_num_operators_Online);
		m_cp_Online.clear();
		m_cp_Online.resize(m_num_operators_Online);

		m_best_inds.clear();
		m_best_inds.emplace_back(new IndividualType(ind));
		calculate_initial_max_violation(m_best_inds);
		m_e = m_max_G;
		m_k = 0;

		int iter_cnt = 0;
		//auto is_dominated = Pareto_compare(&m_best_ind_Online, &best_ind_Online_tmp);//m_best_ind_Offline
		double max_delayTime = 0;
		auto averge_delayTime = (std::numeric_limits<Real>::max)();
		std::cout << "Before Online ALS" << "\n";
		if (DVRP_CAST(pro)->numberObjectives() == 3) {
			std::cout << ind.objective()[0] << "\t" << ind.objective()[2] << "\t" << ind.objective()[1] << "\t" << ind.constraint()[0] << "\t" << ind.constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		}
		if (DVRP_CAST(pro)->numberObjectives() == 1) {
			std::cout << ind.objective()[0] << "\t" << ind.constraint()[3] << "\t" << ind.constraint()[2] << "\t" << ind.constraint()[0] << "\t" << ind.constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		}
		while (iter_cnt < m_OnlineIters) {
			Online_ALS(ind, i, pro, alg, rnd);
			iter_cnt++;
			/*max_delayTime = 0;
			for (auto &i : m_best_ind_Online.variable().m_hi.delay_time) {
				if (max_delayTime <= *std::max_element(i.begin(), i.end()))
					max_delayTime = *std::max_element(i.begin(), i.end());
			}*/
			if (DVRP_CAST(pro)->numberObjectives() == 3) {
				averge_delayTime = ind.objective()[1] / ind.variable().m_cus_order.size();//(DVRP_CAST(pro)->get_current_cus_order().size() + (DVRP_CAST(pro)->get_num_futureOrder() - DVRP_CAST(pro)->get_future_cus_order().size()));;//m_best_ind_Online.constraint()[2] / (DVRP_CAST(pro)->get_current_cus_order().size() + (DVRP_CAST(pro)->get_num_futureOrder() - DVRP_CAST(pro)->get_future_cus_order().size()));
			}
			if (DVRP_CAST(pro)->numberObjectives() == 1) {
				averge_delayTime = ind.constraint()[2] / ind.variable().m_cus_order.size();//(DVRP_CAST(pro)->get_current_cus_order().size() + (DVRP_CAST(pro)->get_num_futureOrder() - DVRP_CAST(pro)->get_future_cus_order().size()));;//m_best_ind_Online.constraint()[2] / (DVRP_CAST(pro)->get_current_cus_order().size() + (DVRP_CAST(pro)->get_num_futureOrder() - DVRP_CAST(pro)->get_future_cus_order().size()));
			}
			if (averge_delayTime > 45 && iter_cnt > 0.9 * m_OnlineIters) {
				if (DVRP_CAST(pro)->get_num_futureOrder() - DVRP_CAST(pro)->get_future_cus_order().size() == 0)
					break;
				auto it = std::find(ind.variable().m_cus_order[i].begin(), ind.variable().m_cus_order[i].end(), new_order.first);
				auto newOrder_pos = std::distance(std::begin(ind.variable().m_cus_order[i]), it);
				ind.variable().m_cus_order[i].erase(ind.variable().m_cus_order[i].begin() + newOrder_pos);
				DVRP_CAST(pro)->get_future_cus_order().push_back(new_order);
				auto fSize = DVRP_CAST(pro)->get_future_cus_order().size();
				if (ind.variable().m_cus_order[i][0] == ind.variable().m_cus_order[i][1]) {
					std::cout << "error" << "\n";
				}
				DVRP_CAST(pro)->construct_path(ind);
				ind.evaluate(pro, alg);
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
		if (DVRP_CAST(pro)->numberObjectives() == 3) {
			std::cout << ind.objective()[0] << "\t" << ind.objective()[2] << "\t" << ind.objective()[1] << "\t" << ind.constraint()[0] << "\t" << ind.constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		}
		if (DVRP_CAST(pro)->numberObjectives() == 1) {
			std::cout << ind.objective()[0] << "\t" << ind.constraint()[3] << "\t" << ind.constraint()[2] << "\t" << ind.constraint()[0] << "\t" << ind.constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << "\t" << ind.variable().m_cus_order.size() << std::endl;
		}
		DVRP_CAST(pro)->set_Online(true);
	}

	void DCNSGAII_LS_pop::Online_ALS(IndividualType &ind, size_t i, Problem *pro, Algorithm *alg, Random *rnd)
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
		DCMOEA_ind<Solution<DVRP::routes>> best_ind_tmp = ind;

		Dominance is_dominated;
		size_t operator_id = 0;
		m_evo_success_Online = false;


		Real p = rnd->uniform.nextNonStd<Real>(0, 1);
		if (p <= m_cp_Online[0]) {
			ls_Rswap2_Online(ind, i, pro, rnd);
			operator_id = 0;
		}
		if (p > m_cp_Online[0] && p <= m_cp_Online[1]) {
			ls_maxDT_Online(ind, i, pro);
			operator_id = 1;
		}
		if (p > m_cp_Online[1] && p <= m_cp_Online[2]) {
			ls_2opt_Online(ind, i, pro, rnd);
			operator_id = 2;
		}
		DVRP_CAST(pro)->construct_path(ind);
		ind.evaluate(pro, alg);

		m_best_inds.clear();
		m_best_inds.emplace_back(new IndividualType(ind));
		calculate_violation_objective(m_best_inds);
		mark_Solution_efeasible(m_best_inds);
		ind = *m_best_inds[0];

		is_dominated = Pareto_compare(&ind, &best_ind_tmp);
		//if (m_convergence_Online) {
		//	if (is_dominated == Dominance::Dominated || is_dominated == Dominance::Equal || is_dominated == Dominance::Non_comparable) {
		//		ind = best_ind_tmp;
		//		//m_evo_success = false;
		//	}
		//	else {
		//		m_evo_success_Online = true;
		//		m_cnt_fail_Online = 0;
		//	}
		//	if (is_dominated == Dominance::Dominant)
		//		m_score_Online[operator_id] = m_score_Online[operator_id] + 1;
		//}
		//else {
		if (is_dominated != Dominance::Dominant) {
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
		if (is_dominated == Dominance::Equal)
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

	void DCNSGAII_LS_pop::ls_Rswap2_Online(IndividualType &ind, size_t i, Problem *pro, Random *rnd)
	{
		auto i_t = ind;
		int selected_car = i;
		int begin_pos = -1;
		for (int j = 1; j < ind.variable().m_cus_order[i].size(); ++j) {
			if (DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::depot)
				continue;
			if (DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->datum_t.is_served == false)
			{
				begin_pos = j + 1;
				break;
			}
		}
		if (begin_pos == -1 || begin_pos == 0 || begin_pos == 1 || begin_pos == ind.variable().m_cus_order[selected_car].size() - 1 || begin_pos == ind.variable().m_cus_order[selected_car].size() - 2)
			return;
		//随机选1个顾客
		int selected_pos1 = rnd->uniform.nextNonStd<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);

		int selected_pos2 = rnd->uniform.nextNonStd<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);
		while (selected_pos1 == selected_pos2)
			selected_pos2 = rnd->uniform.nextNonStd<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);

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

	void DCNSGAII_LS_pop::ls_maxDT_Online(IndividualType &ind, size_t i, Problem *pro)
	{
		auto i_t = ind;
		bool flag = true;
		if (ind.objective().size() <= 2) {
			if (ind.constraint()[2] == 0)
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
				if (DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::depot)
					continue;
				if (DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->datum_t.is_served == false)
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
			if (ind.variable().m_cus_order[i][0] != DVRP_CAST(pro)->get_depot()) {
				std::cout << "ls_maxDT_Online() error" << "\n";
			}
			if (ind.variable().m_cus_order[i][0] == ind.variable().m_cus_order[i][1]) {
				std::cout << "ls_maxDT_Online() error" << "\n";
			}
		}


	}

	void DCNSGAII_LS_pop::ls_2opt_Online(IndividualType &ind, size_t i, Problem *pro, Random *rnd)
	{
		auto i_t = ind;
		if (ind.variable().m_cus_order[i][0] != DVRP_CAST(pro)->get_depot()) {
			std::cout << "ls_2opt_Online() error" << "\n";
		}
		int selected_car = i;
		int begin_pos = -1;
		for (int j = 1; j < ind.variable().m_cus_order[i].size(); ++j) {
			if (DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->type == DVRP::node_type::depot)
				continue;
			if (DVRP_CAST(pro)->get_net().get_road_net()[ind.variable().m_cus_order[i][j]]->datum_t.is_served == false)
			{
				begin_pos = j + 1;
				break;
			}
		}
		if (begin_pos == -1 || begin_pos == 0 || begin_pos == 1 || begin_pos == ind.variable().m_cus_order[selected_car].size() - 1 || begin_pos == ind.variable().m_cus_order[selected_car].size() - 2)
			return;
		//随机选1个顾客
		int selected_pos1 = rnd->uniform.nextNonStd<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);

		int selected_pos2 = rnd->uniform.nextNonStd<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);
		while (selected_pos1 == selected_pos2)
			selected_pos2 = rnd->uniform.nextNonStd<int>(begin_pos, ind.variable().m_cus_order[selected_car].size() - 1);

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
		ind.variable().m_cus_order[selected_car].push_back(DVRP_CAST(pro)->get_depot());

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
		ind.variable().m_cus_order[selected_car].push_back(DVRP_CAST(pro)->get_depot());
		if (ind.variable().m_cus_order[i][0] != DVRP_CAST(pro)->get_depot()) {
			std::cout << "ls_2opt_Online() error" << "\n";
		}
		if (ind.variable().m_cus_order[selected_car][0] == ind.variable().m_cus_order[selected_car][1]) {
			std::cout << "ls_2opt_Online() error" << "\n";
		}
	}

	DCNSGAII_LS_pop::DCNSGAII_LS_pop(size_t size_pop, Problem *pro) :
		Population(size_pop, pro), 
		m_rand_seq(size_pop),
		m_best_ind_Online(pro->numberObjectives(), pro->numberConstraints()),
		m_best_ind_Offline(pro->numberObjectives(), pro->numberConstraints())
	{
		for (auto &i : m_individuals) {
			i->resize_vio_obj(m_num_vio_obj);
			m_offspring.emplace_back(new IndividualType(*i));
			m_temp_pop.emplace_back(new IndividualType(*i));
		}
		std::iota(m_rand_seq.begin(), m_rand_seq.end(), 0);
	}

	int DCNSGAII_LS_pop::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		int tag = kNormalEval;

		if (judge_population_efeasible(m_individuals)) {
			m_k++;
			if (m_k > m_max_K + 1) {
				return Terminate;
				//m_flag = true;
			}
			reduce_boundary();
			mark_Solution_efeasible(m_individuals);
		}

		//generate_offspring_LS();
		generate_offspring_ALS(pro, alg, rnd);

		//evaluate offspring pop
		std::vector<IndividualType *> pop;
		for (auto &i : m_offspring)
			i->evaluate(pro, alg);
		calculate_violation_objective(m_offspring);
		mark_Solution_efeasible(m_offspring);
		for (auto &i : m_individuals)
			pop.emplace_back(i.get());
		for (int i = 0; i < m_individuals.size(); ++i) {
			//添加子代是否和父代个体完全相同的判断 
			Dominance is_dominated = Pareto_compare(m_offspring[i].get(), m_individuals[i].get());
			if (is_dominated == Dominance::Dominant || is_dominated == Dominance::NonDominated) {
				pop.emplace_back(m_offspring[i].get());
			}
		}
		//select next pop
		select_next_parent_population(pop, pro);
		if (tag == kNormalEval)
			m_iteration++;
		return tag;
	}

	void DCNSGAII_LS_pop::select_next_parent_population(std::vector<IndividualType *> &pop, Problem *pro) {
		/* Nondominated sorting based on e-Pareto domination  */
		std::vector<int> ranks;
		std::function<ofec::Dominance(IndividualType *const &, IndividualType *const &)> comp = std::bind(&DCNSGAII_LS_pop::e_Pareto_compare, this, _1, _2);
		nd_sort::fast_sort<IndividualType *>(pop, ranks, comp);
		for (size_t i = 0; i < pop.size(); i++) {
			pop[i]->setFitness(ranks[i]);
		}

		size_t cur_rank = 0;
		size_t id_ind = 0;
		while (true) {
			int count = 0;
			for (size_t i = 0; i < pop.size(); i++)
				if (pop[i]->fitness() == cur_rank)
					count++;
			int size2 = id_ind + count;
			if (size2 > this->m_individuals.size()) {
				break;
			}
			for (size_t i = 0; i < pop.size(); i++)
				if (pop[i]->fitness() == cur_rank) {
					*m_temp_pop[id_ind] = *pop[i];
					++id_ind;
				}
			cur_rank++;
			if (id_ind >= this->m_individuals.size()) break;
		}
		if (id_ind < pop.size()) {
			std::vector<int> list;	// save the Solutions in the overflowed front
			for (size_t i = 0; i < pop.size(); i++)
				if (pop[i]->fitness() == cur_rank)
					list.push_back(i);
			int s2 = list.size();
			std::vector<Real> density(s2);
			std::vector<Real> obj(s2);
			std::vector<int> idx(s2);
			std::vector<int> idd(s2);
			for (size_t i = 0; i < s2; i++) {
				idx[i] = i;
				density[i] = 0;
			}
			for (size_t j = 0; j < 2; j++) {
				for (size_t i = 0; i < s2; i++) {
					idd[i] = i;
					obj[i] = pop[list[i]]->get_vio_obj()[0];
				}
				merge_sort(obj, s2, idd, true, 0, s2 - 1, s2);
				density[idd[0]] += -1.0e+30;
				density[idd[s2 - 1]] += -1.0e+30;
				for (int k = 1; k < s2 - 1; k++)
					density[idd[k]] += -(obj[idd[k]] - obj[idd[k - 1]] + obj[idd[k + 1]] - obj[idd[k]]);
			}
			for (size_t j = 0; j < DVRP_CAST(pro)->numberObjectives(); j++) {
				for (size_t i = 0; i < s2; i++) {
					idd[i] = i;
					obj[i] = pop[list[i]]->objective()[j];
				}
				merge_sort(obj, s2, idd, true, 0, s2 - 1, s2);
				density[idd[0]] += -1.0e+30;
				density[idd[s2 - 1]] += -1.0e+30;
				for (int k = 1; k < s2 - 1; k++)
					density[idd[k]] += -(obj[idd[k]] - obj[idd[k - 1]] + obj[idd[k + 1]] - obj[idd[k]]);
			}
			idd.clear();
			obj.clear();
			int s3 = this->m_individuals.size() - id_ind;
			merge_sort(density, s2, idx, true, 0, s2 - 1, s3);
			for (size_t i = 0; i < s3; i++) {
				*m_temp_pop[id_ind] = *pop[list[idx[i]]];
				++id_ind;
			}
			density.clear();
			idx.clear();
			list.clear();
		}
		for (size_t i = 0; i < m_temp_pop.size(); i++) {
			*m_individuals[i] = *m_temp_pop[i];
		}
	}

	Dominance DCNSGAII_LS_pop::e_Pareto_compare(IndividualType *const &s1, IndividualType *const &s2) {
		/* One efeasible one in-efeasible */
		if (s1->get_efeasible() != s2->get_efeasible()) {
			if (s1->get_efeasible())
				return Dominance::Dominant;
			else
				return Dominance::Dominated;
		}

		/* Both efeasible */
		else if (s1->get_efeasible() && s2->get_efeasible()) {
			auto nor_obj_result = objectiveCompare<Real>(s1->objective(), s2->objective(), OptimizeMode::kMinimize);
			auto vio_obj_result = objectiveCompare<Real>(s1->get_vio_obj(), s2->get_vio_obj(), OptimizeMode::kMinimize);

			if (nor_obj_result == Dominance::Dominant && vio_obj_result == Dominance::Equal)
				return Dominance::Dominant;
			if (nor_obj_result == Dominance::Equal && vio_obj_result == Dominance::Dominant)
				return Dominance::Dominant;
			if (nor_obj_result == Dominance::Dominant && vio_obj_result == Dominance::Dominant)
				return Dominance::Dominant;

			if (nor_obj_result == Dominance::Dominated && vio_obj_result == Dominance::Equal)
				return Dominance::Dominated;
			if (nor_obj_result == Dominance::Equal && vio_obj_result == Dominance::Dominated)
				return Dominance::Dominated;
			if (nor_obj_result == Dominance::Dominated && vio_obj_result == Dominance::Dominated)
				return Dominance::Dominated;

			if (nor_obj_result == Dominance::Dominated && vio_obj_result == Dominance::Dominant)
				return Dominance::NonDominated;
			if (nor_obj_result == Dominance::Dominant && vio_obj_result == Dominance::Dominated)
				return Dominance::NonDominated;
			if (nor_obj_result == Dominance::NonDominated || vio_obj_result == Dominance::NonDominated)
				return Dominance::NonDominated;

			if (nor_obj_result == Dominance::Equal && vio_obj_result == Dominance::Equal)
				return Dominance::Equal;
		}
		//return objectiveCompare<Real>(s1->objective(), s2->objective(), CONTINUOUS_CAST->opt_mode());

	/* Both in-efeasible */
		else {
			if (s1->get_vio_obj()[0] < s2->get_vio_obj()[0])
				return Dominance::Dominant;
			else if (s1->get_vio_obj()[0] > s2->get_vio_obj()[0])
				return Dominance::Dominated;
			else
				return Dominance::Equal;
		}
	}

	Dominance DCNSGAII_LS_pop::Pareto_compare(IndividualType *const &s1, IndividualType *const &s2)
	{
		auto nor_obj_result = objectiveCompare<Real>(s1->objective(), s2->objective(), OptimizeMode::kMinimize);
		auto vio_obj_result = objectiveCompare<Real>(s1->get_vio_obj(), s2->get_vio_obj(), OptimizeMode::kMinimize);

		if (nor_obj_result == Dominance::Dominant && vio_obj_result == Dominance::Equal)
			return Dominance::Dominant;
		if (nor_obj_result == Dominance::Equal && vio_obj_result == Dominance::Dominant)
			return Dominance::Dominant;
		if (nor_obj_result == Dominance::Dominant && vio_obj_result == Dominance::Dominant)
			return Dominance::Dominant;

		if (nor_obj_result == Dominance::Dominated && vio_obj_result == Dominance::Equal)
			return Dominance::Dominated;
		if (nor_obj_result == Dominance::Equal && vio_obj_result == Dominance::Dominated)
			return Dominance::Dominated;
		if (nor_obj_result == Dominance::Dominated && vio_obj_result == Dominance::Dominated)
			return Dominance::Dominated;

		if (nor_obj_result == Dominance::Dominated && vio_obj_result == Dominance::Dominant)
			return Dominance::NonDominated;
		if (nor_obj_result == Dominance::Dominant && vio_obj_result == Dominance::Dominated)
			return Dominance::NonDominated;
		if (nor_obj_result == Dominance::NonDominated || vio_obj_result == Dominance::NonDominated)
			return Dominance::NonDominated;

		if (nor_obj_result == Dominance::Equal && vio_obj_result == Dominance::Equal)
			return Dominance::Equal;
	}

	void DCNSGAII_LS_pop::generate_offspring_LS(Problem *pro, Algorithm *alg, Random *rnd) {
		//local search
		for (int i = 0; i < m_individuals.size(); ++i) {
			m_offspring[i] = *m_individuals[i];
			m_offspring[i]->set_vio_obj(0, m_individuals[i]->get_vio_obj()[0]);
			m_offspring[i]->set_efeasible(m_individuals[i]->get_efeasible());
		}
		Dominance is_dominated;
		//size_t operator_id = 0;
		m_evo_success = false;
		m_cnt_fail = 0;
		m_evo_success = false;
		for (size_t i = 0; i < m_offspring.size(); ++i) {
			calculate_cp();
			Real p = rnd->uniform.nextNonStd<Real>(0, 1);
			if (p <= m_cp[0]) {
				ls_Rswap2(*m_offspring[i], pro, rnd);
				//operator_id = 0;
			}
			if (p > m_cp[0] && p <= m_cp[1]) {
				ls_tsp(*m_offspring[i], pro, rnd);
				//operator_id = 1;
			}
			if (p > m_cp[1] && p <= m_cp[2]) {
				ls_maxDT(*m_offspring[i], pro);
				//operator_id = 2;
			}
			if (p > m_cp[2] && p <= m_cp[3]) {
				ls_r1(*m_offspring[i], pro, rnd);
				//operator_id = 3;
			}
			if (p > m_cp[3] && p <= m_cp[4]) {
				ls_rMulti(*m_offspring[i], pro, rnd);
				//operator_id = 4;
			}
			if (p > m_cp[4] && p <= m_cp[5]) {
				ls_2opt(*m_offspring[i], pro, rnd);
				//operator_id = 5;
			}
			if (p > m_cp[5] && p <= m_cp[6]) {
				ls_maxWTlength(*m_offspring[i], pro);
				//operator_id = 6;
			}
			if (p > m_cp[6] && p <= m_cp[7]) {
				ls_maxLength(*m_offspring[i], pro, rnd);
				//operator_id = 7;
			}
			check_order(*m_offspring[i], pro);
			m_offspring[i]->evaluate(pro, alg);
			calculate_violation_objective(m_offspring);
			mark_Solution_efeasible(m_offspring);

			is_dominated = Pareto_compare(m_offspring[i].get(), m_individuals[i].get());
			auto a = double(i);
			if (m_convergence) {//double(i) > (int(m_individuals.size()) / 1.5)
				if (is_dominated == Dominance::Dominated || is_dominated == Dominance::Equal || is_dominated == Dominance::NonComparable) {
					m_offspring[i] = *m_individuals[i];
					m_offspring[i]->set_vio_obj(0, m_individuals[i]->get_vio_obj()[0]);
					m_offspring[i]->set_efeasible(m_individuals[i]->get_efeasible());
					m_evo_success = false;
				}
				else {
					m_evo_success = true;
					//m_cnt_fail = 0;
				}
				/*if (is_dominated == Dominance::Dominant)
					m_score[operator_id] = m_score[operator_id] + 1;*/
			}
			else {
				if (is_dominated != Dominance::Dominant) {
					m_offspring[i] = *m_individuals[i];
					m_offspring[i]->set_vio_obj(0, m_individuals[i]->get_vio_obj()[0]);
					m_offspring[i]->set_efeasible(m_individuals[i]->get_efeasible());
					m_evo_success = false;
				}
				else {
					m_evo_success = true;
					//m_cnt_fail = 0;
					//m_score[operator_id] = m_score[operator_id] + 1;
				}

			}
			is_dominated = Pareto_compare(m_offspring[i].get(), m_individuals[i].get());
			if (is_dominated == Dominance::Equal)
				m_cnt_fail++;
			//auto decrease_level = int(m_individuals.size())*exp(-pow(int(m_cnt_fail), 0.5));//0.3*int(m_individuals.size()) + 1 * 
			auto decrease_level = int(m_individuals.size()) / 2.0;
			if (double(m_cnt_fail) > decrease_level)
				m_convergence = true;

			if (m_convergence) {
				//m_score.clear();
				//m_score.resize(m_num_operators, 1);
				if (m_evo_success) {
					//std::cout << "jump out local optimua successfully" << std::endl;
					m_evo_success = false;
					m_cnt_fail = 0;
					m_convergence = false;
				}
			}
		}

	}
	void DCNSGAII_LS_pop::generate_offspring_ALS(Problem *pro, Algorithm *alg, Random *rnd) {
		std::string scoreFname;
		scoreFname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/score_" + to_string(m_individuals.size()) + "_order" + to_string(DVRP_CAST(pro)->get_current_cus_order().size()) + "ALSDCMOEA_VRPRTC_AHC_initial.txt";
		std::ofstream score(scoreFname, ios::app);
		std::stringstream score_;
		score_ << "iteration: " << to_string(m_iteration) << "\n";

		//adaptive local search
		for (int i = 0; i < m_individuals.size(); ++i) {
			m_offspring[i] = *m_individuals[i];
			m_offspring[i]->set_vio_obj(0, m_individuals[i]->get_vio_obj()[0]);
			m_offspring[i]->set_efeasible(m_individuals[i]->get_efeasible());
		}
		Dominance is_dominated;
		size_t operator_id = 0;
		m_evo_success = false;
		m_cnt_fail = 0;
		m_evo_success = false;
		for (size_t i = 0; i < m_offspring.size(); ++i) {
			calculate_cp();
			Real p = rnd->uniform.nextNonStd<Real>(0, 1);
			if (p <= m_cp[0]) {
				ls_Rswap2(*m_offspring[i], pro, rnd);
				operator_id = 0;//1
			}
			if (p > m_cp[0] && p <= m_cp[1]) {
				ls_tsp(*m_offspring[i], pro, rnd);
				operator_id = 1;//3
			}
			if (p > m_cp[1] && p <= m_cp[2]) {
				ls_maxDT(*m_offspring[i], pro);
				operator_id = 2;//6
			}
			if (p > m_cp[2] && p <= m_cp[3]) {
				ls_r1(*m_offspring[i], pro, rnd);
				operator_id = 3;//4
			}
			if (p > m_cp[3] && p <= m_cp[4]) {
				ls_rMulti(*m_offspring[i], pro, rnd);
				operator_id = 4;//5
			}
			if (p > m_cp[4] && p <= m_cp[5]) {
				ls_2opt(*m_offspring[i], pro, rnd);
				operator_id = 5;//2
			}
			if (p > m_cp[5] && p <= m_cp[6]) {
				ls_maxWTlength(*m_offspring[i], pro);
				operator_id = 6;//7
			}
			if (p > m_cp[6] && p <= m_cp[7]) {
				ls_maxLength(*m_offspring[i], pro, rnd);
				operator_id = 7;//8
			}
			check_order(*m_offspring[i], pro);
			m_offspring[i]->evaluate(pro, alg);
			calculate_violation_objective(m_offspring);
			mark_Solution_efeasible(m_offspring);

			is_dominated = Pareto_compare(m_offspring[i].get(), m_individuals[i].get());
			if (m_convergence) {//double(i) > (int(m_individuals.size()) / 1.5)
				if (is_dominated == Dominance::Dominated || is_dominated == Dominance::Equal || is_dominated == Dominance::NonComparable) {
					m_offspring[i] = *m_individuals[i];
					m_offspring[i]->set_vio_obj(0, m_individuals[i]->get_vio_obj()[0]);
					m_offspring[i]->set_efeasible(m_individuals[i]->get_efeasible());
					m_evo_success = false;
				}
				else {
					m_evo_success = true;
					//m_cnt_fail = 0;
				}
				if (is_dominated == Dominance::Dominant)
					m_score[operator_id] = m_score[operator_id] + 1;
			}
			else {
				if (is_dominated != Dominance::Dominant) {
					m_offspring[i] = *m_individuals[i];
					m_offspring[i]->set_vio_obj(0, m_individuals[i]->get_vio_obj()[0]);
					m_offspring[i]->set_efeasible(m_individuals[i]->get_efeasible());
					m_evo_success = false;
				}
				else {
					m_evo_success = true;
					//m_cnt_fail = 0;
					m_score[operator_id] = m_score[operator_id] + 1;
				}

			}
			is_dominated = Pareto_compare(m_offspring[i].get(), m_individuals[i].get());
			if (is_dominated == Dominance::Equal)
				m_cnt_fail++;
			//auto decrease_level = int(m_individuals.size())*exp(-pow(int(m_cnt_fail), 0.5));//0.3*int(m_individuals.size()) + 1 * 
			auto decrease_level = int(m_individuals.size()) / 2.0;
			if (double(m_cnt_fail) > decrease_level)
				m_convergence = true;

			if (m_convergence) {
				if (!m_evo_success) {
					//score_ << "before reset:" << "\n";
					for (int ii = 0; ii < m_score.size(); ++ii) {
						score_ << m_score[ii] << "\t";
					}
					score_ << "\n";
				}
				m_score.clear();
				m_score.resize(m_num_operators, 1);
				/*if (!m_evo_success) {
					score_ << "after reset:" << "\n";
					for (int ii = 0; ii < m_score.size(); ++ii) {
						score_ << m_score[ii] << "\t";
					}
					score_ << "\n";
				}*/
				if (m_evo_success) {
					//std::cout << "jump out local optimua successfully" << std::endl;
					m_evo_success = false;
					m_cnt_fail = 0;
					m_convergence = false;
				}
			}
		}
		score << score_.str();
		score.close();
	}

	void DCNSGAII_LS_pop::calculate_cp()
	{
		Real sum_score = 0;
		for (auto &i : m_score)
			sum_score += i;
		for (int i = 0; i < m_num_operators; ++i) {
			m_p_select[i] = Real(Real(m_score[i]) / sum_score);
		}
		for (auto &i : m_cp)
			i = 0;
		m_cp[0] = m_p_select[0];
		for (int i = 1; i < m_num_operators; ++i) {
			m_cp[i] += (m_p_select[i] + m_cp[i - 1]);
		}
	}

	//void DCNSGAII_LS_pop::generate_offspring_ALS_Online()
	//{
	//	//adaptive local search
	//	Real sum_score = 0;
	//	for (auto &i : m_score_Online)
	//		sum_score += i;
	//	for (int i = 0; i < m_num_operators_Online; ++i) {
	//		m_p_select_Online[i] = Real(Real(m_score_Online[i]) / sum_score);
	//	}
	//	for (auto &i : m_cp_Online)
	//		i = 0;
	//	m_cp_Online[0] = m_p_select_Online[0];
	//	for (int i = 1; i < m_num_operators_Online; ++i) {
	//		m_cp_Online[i] += (m_p_select_Online[i] + m_cp_Online[i - 1]);
	//	}
	//	DCMOEA_ind<Solution<DVRP::routes>> best_ind_tmp = m_best_ind_Online;

	//	Dominance is_dominated;
	//	size_t operator_id = 0;
	//	m_evo_success_Online = false;


	//	Real p = rnd->uniform.nextNonStd<Real>(0, 1);
	//	if (p <= m_cp_Online[0]) {
	//		ls_Rswap2(m_best_ind_Online);
	//		operator_id = 0;
	//	}
	//	if (p > m_cp_Online[0] && p <= m_cp_Online[1]) {
	//		ls_tsp(m_best_ind_Online);
	//		operator_id = 1;
	//	}
	//	if (p > m_cp_Online[1] && p <= m_cp_Online[2]) {
	//		ls_maxDT(m_best_ind_Online);
	//		operator_id = 2;
	//	}
	//	if (p > m_cp_Online[2] && p <= m_cp_Online[3]) {
	//		ls_r1(m_best_ind_Online);
	//		operator_id = 3;
	//	}
	//	if (p > m_cp_Online[3] && p <= m_cp_Online[4]) {
	//		ls_rMulti(m_best_ind_Online);
	//		operator_id = 4;
	//	}
	//	if (p > m_cp_Online[4] && p <= m_cp_Online[5]) {
	//		ls_2opt(m_best_ind_Online);
	//		operator_id = 5;
	//	}
	//	if (p > m_cp_Online[5] && p <= m_cp_Online[6]) {
	//		ls_maxWTlength(m_best_ind_Online);
	//		operator_id = 6;
	//	}
	//	if (p > m_cp_Online[6] && p <= m_cp_Online[7]) {
	//		ls_maxLength(m_best_ind_Online);
	//		operator_id = 7;
	//	}
	//	/*if (p > m_cp_Online[7] && p <= m_cp_Online[8]) {
	//		insert_nOrder(m_best_ind_Online);
	//		operator_id = 8;
	//	}*/

	//	//标记已被服务
	//	for (int i = 0; i < m_best_ind_Online.variable().m_cus_order.size(); ++i) {
	//		for (int m = 1; m < m_best_ind_Online.variable().m_cus_order[i].size() - 1; ++m) {
	//			DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][m]]->datum_t.is_served = false;
	//		}
	//		for (int n = 1; n < m_best_ind_Online.variable().m_cus_order[i].size() - 1; ++n) {
	//			if (DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][n]]->type == DVRP::node_type::future_customer)
	//				break;
	//			DVRP_CAST(pro)->get_net().get_road_net()[m_best_ind_Online.variable().m_cus_order[i][n]]->datum_t.is_served = true;
	//		}
	//	}

	//	m_best_ind_Online.evaluate();

	//	m_best_inds.clear();
	//	m_best_inds.emplace_back(new IndividualType(m_best_ind_Online));
	//	calculate_violation_objective(m_best_inds);
	//	mark_Solution_efeasible(m_best_inds);
	//	m_best_ind_Online = *m_best_inds[0];

	//	is_dominated = e_Pareto_compare(&m_best_ind_Online, &best_ind_tmp);
	//	if (m_convergence_Online) {
	//		if (is_dominated == Dominance::Dominated || is_dominated == Dominance::Equal || is_dominated == Dominance::Non_comparable) {
	//			m_best_ind_Online = best_ind_tmp;
	//			//m_evo_success = false;
	//		}
	//		else {
	//			m_evo_success_Online = true;
	//			m_cnt_fail_Online = 0;
	//		}
	//		if (is_dominated == Dominance::Dominant)
	//			m_score_Online[operator_id] = m_score_Online[operator_id] + 1;
	//	}
	//	else {
	//		if (is_dominated != Dominance::Dominant) {
	//			m_best_ind_Online = best_ind_tmp;
	//			//m_evo_success = false;
	//		}
	//		else {
	//			m_evo_success_Online = true;
	//			m_cnt_fail_Online = 0;
	//			m_score_Online[operator_id] = m_score_Online[operator_id] + 1;
	//		}

	//	}
	//	is_dominated = e_Pareto_compare(&m_best_ind_Online, &best_ind_tmp);
	//	if (is_dominated == Dominance::Equal)
	//		m_cnt_fail_Online++;
	//	//auto b = exp(-pow(int(m_cnt_fail),0.5));
	//	auto decrease_level = 100;//50 * exp(-pow(int(m_cnt_fail_Online), 0.5));//0.3*int(m_individuals.size()) + 1 * 
	//	if (double(m_cnt_fail_Online) >= decrease_level)
	//		m_convergence_Online = true;

	//	if (m_convergence_Online) {
	//		m_score_Online.clear();
	//		m_score_Online.resize(m_num_operators_Online, 1);
	//		if (m_evo_success_Online) {
	//			std::cout << "jump out local optimua successfully" << std::endl;
	//			m_evo_success_Online = false;
	//			m_convergence_Online = false;
	//			m_cnt_fail_Online = 0;
	//		}
	//	}

	//}

	void DCNSGAII_LS_pop::tsp_nearest(std::vector<size_t> &member_index, Real presentTime, Problem *pro, Random *rnd)
	{
		std::vector<size_t> tsp;
		std::vector<std::tuple<size_t, Real, Real>> member_cost;
		size_t start = DVRP_CAST(pro)->get_depot();
		Real Time = presentTime;
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

				get_cost(start, member_index[i], present_time,pro,rnd);//此处耗时可能太长

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
			Real co_d = rnd->uniform.nextNonStd<Real>(0, 1);
			for (int i = 0; i < member_index.size(); ++i) {
				if (member_index.size() <= 1) {
					conflict_tw_min = 0;
					conflict_d_min = 0;
				}
				conflict[i] = co_d * ((conflict_d[i] - conflict_d_min) / (conflict_d_max - conflict_d_min)) + \
					(1 - co_d) * ((conflict_tw[i] - conflict_tw_min) / (conflict_tw_max - conflict_tw_min));
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
			auto it_check = std::find(member_index.begin(), member_index.end(), *(member_index.begin() + index));
			if (it_check == member_index.end()) {
				system("pause");
			}
			member_index.erase(member_index.begin() + index);
		}
		member_index = tsp;
	}
	void DCNSGAII_LS_pop::get_cost(size_t s, size_t e, Real &present_time, Problem *pro, Random *rnd)
	{
		//present_time = present_time + DVRP_CAST(pro)->get_cost_time(present_time, s, e);
		//present_time = present_time + DVRP_CAST(pro)->get_cost_time()[s][e];
		present_time = present_time + DVRP_CAST(pro)->get_cost_time(present_time, s, e);
	}
	void DCNSGAII_LS_pop::ls_Rswap2(IndividualType &offspring, Problem *pro, Random *rnd)
	{
		int selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		while (offspring.variable().m_cus_order[selected_car].size() <= 3) {
			selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		}
		//随机选1个顾客
		int selected_pos1 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);

		int selected_pos2 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);
		while (selected_pos1 == selected_pos2)
			selected_pos2 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);

		if (selected_pos1 > selected_pos2) {
			selected_pos1 = selected_pos1 + selected_pos2;
			selected_pos2 = selected_pos1 - selected_pos2;
			selected_pos1 = selected_pos1 - selected_pos2;
		}

		if (selected_pos2 < offspring.variable().m_cus_order[selected_car].size() - 1) {
			int customer = offspring.variable().m_cus_order[selected_car][selected_pos2];
			offspring.variable().m_cus_order[selected_car][selected_pos2] = offspring.variable().m_cus_order[selected_car][selected_pos1];
			offspring.variable().m_cus_order[selected_car][selected_pos1] = customer;
			DVRP_CAST(pro)->construct_path(offspring);
		}
		if (offspring.variable().m_cus_order[selected_car][0] != DVRP_CAST(pro)->get_depot()) {
			std::cout << "ls_Rswap2() error" << "\n";
		}
	}

	void DCNSGAII_LS_pop::ls_tsp(IndividualType &offspring, Problem *pro, Random *rnd)
	{
		int selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		while (offspring.variable().m_cus_order[selected_car].size() <= 3) {
			selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		}
		//随机选1个顾客
		int selected_pos1 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);

		int selected_pos2 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);
		while (selected_pos1 == selected_pos2 || abs(selected_pos1 - selected_pos2) > 4)
			selected_pos2 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);

		if (selected_pos1 > selected_pos2) {
			selected_pos1 = selected_pos1 + selected_pos2;
			selected_pos2 = selected_pos1 - selected_pos2;
			selected_pos1 = selected_pos1 - selected_pos2;
		}

		std::vector<size_t> seq_1, seq_2, seq_3;
		Real beginTime = 0;
		beginTime = DVRP_CAST(pro)->get_depart_time();
		/*s2_beginTime = offspring->variable().m_hi.arrive_time[selected_car][selected_pos1 - 1];
		s3_beginTime = offspring->variable().m_hi.arrive_time[selected_car][selected_pos2 - 1];*/
		for (int i = 1; i < selected_pos1; ++i) {
			seq_1.push_back(offspring.variable().m_cus_order[selected_car][i]);
		}
		for (int i = selected_pos1; i < selected_pos2; ++i) {
			seq_2.push_back(offspring.variable().m_cus_order[selected_car][i]);
		}
		for (int i = selected_pos2; i < offspring.variable().m_cus_order[selected_car].size() - 1; ++i) {
			seq_3.push_back(offspring.variable().m_cus_order[selected_car][i]);
		}
		offspring.variable().m_cus_order[selected_car].clear();
		offspring.variable().m_cus_order[selected_car].push_back(DVRP_CAST(pro)->get_depot());

		if (seq_1.size() != 0) {
			for (auto &i : seq_1)
				offspring.variable().m_cus_order[selected_car].push_back(i);
		}
		if (seq_3.size() != 0) {
			for (auto &i : seq_3)
				offspring.variable().m_cus_order[selected_car].push_back(i);
		}
		offspring.variable().m_cus_order[selected_car].push_back(DVRP_CAST(pro)->get_depot());
		int selected_car2 = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		while (selected_car2 == selected_car) {
			selected_car2 = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		}
		offspring.variable().m_cus_order[selected_car2].pop_back();
		offspring.variable().m_cus_order[selected_car2].erase(offspring.variable().m_cus_order[selected_car2].begin());
		if (seq_2.size() != 0) {
			for (auto &i : seq_2)
				offspring.variable().m_cus_order[selected_car2].push_back(i);
		}
		tsp_nearest(offspring.variable().m_cus_order[selected_car2], beginTime, pro, rnd);
		offspring.variable().m_cus_order[selected_car2].insert(offspring.variable().m_cus_order[selected_car2].begin(), DVRP_CAST(pro)->get_depot());
		offspring.variable().m_cus_order[selected_car2].push_back(DVRP_CAST(pro)->get_depot());
		DVRP_CAST(pro)->construct_path(offspring);
		if (offspring.variable().m_cus_order[selected_car][0] != DVRP_CAST(pro)->get_depot() || offspring.variable().m_cus_order[selected_car2][0] != DVRP_CAST(pro)->get_depot()) {
			std::cout << "ls_tsp() error" << "\n";
		}

	}

	void DCNSGAII_LS_pop::ls_maxDT(IndividualType &offspring, Problem *pro)
	{
		bool flag = true;
		if (offspring.objective().size() == 1) {
			if (offspring.constraint()[2] == 0)
				flag = false;
		}
		if (offspring.objective().size() == 3) {
			if (offspring.objective()[1] == 0)
				flag = false;
		}
		if (flag == true) {
			int max_delay_car = 0;
			Real max_delay_time = 0.0;
			int max_delay_pos = 0;
			for (int j = 0; j < offspring.variable().m_hi.delay_time.size(); ++j) {
				for (int k = 0; k < offspring.variable().m_hi.delay_time[j].size() - 1; ++k) {
					if (max_delay_time < offspring.variable().m_hi.delay_time[j][k]) {
						max_delay_time = offspring.variable().m_hi.delay_time[j][k];
						max_delay_car = j;
						max_delay_pos = k + 1;
					}
				}
			}
			int max_wait_car = 0;
			Real max_wait_time = 0.0;
			int max_wait_pos = 0;
			for (int j = 0; j < offspring.variable().m_hi.wait_time.size(); ++j) {
				for (int k = 0; k < offspring.variable().m_hi.wait_time[j].size() - 1; ++k) {
					if (max_wait_time < offspring.variable().m_hi.wait_time[j][k]) {
						max_wait_time = offspring.variable().m_hi.wait_time[j][k];
						max_wait_car = j;
						max_wait_pos = k + 1;
					}
				}
			}

			if (max_delay_time != 0) {
				int max_delay_cus = offspring.variable().m_cus_order[max_delay_car][max_delay_pos];
				offspring.variable().m_cus_order[max_delay_car].erase(offspring.variable().m_cus_order[max_delay_car].begin() + max_delay_pos);
				offspring.variable().m_cus_order[max_wait_car].insert(offspring.variable().m_cus_order[max_wait_car].begin() + max_wait_pos, max_delay_cus);
				DVRP_CAST(pro)->construct_path(offspring);
			}
			if (offspring.variable().m_cus_order[max_delay_car][0] != DVRP_CAST(pro)->get_depot() || offspring.variable().m_cus_order[max_wait_car][0] != DVRP_CAST(pro)->get_depot()) {
				std::cout << "ls_maxDT() error" << "\n";
			}
		}
	}

	void DCNSGAII_LS_pop::ls_r1(IndividualType &offspring, Problem *pro, Random *rnd)
	{
		int selected_car2 = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		int selected_pos2 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car2].size() - 1);
		int selected_car = selected_car2;
		while (selected_car == selected_car2)
			selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		int selected_pos = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);
		int customer = offspring.variable().m_cus_order[selected_car][selected_pos];
		offspring.variable().m_cus_order[selected_car][selected_pos] = offspring.variable().m_cus_order[selected_car2][selected_pos2];
		offspring.variable().m_cus_order[selected_car2][selected_pos2] = customer;

		DVRP_CAST(pro)->construct_path(offspring);
		if (offspring.variable().m_cus_order[selected_car][0] != DVRP_CAST(pro)->get_depot() || offspring.variable().m_cus_order[selected_car2][0] != DVRP_CAST(pro)->get_depot()) {
			std::cout << "ls_r1() error" << "\n";
		}
	}

	void DCNSGAII_LS_pop::ls_rMulti(IndividualType &offspring, Problem *pro, Random *rnd)
	{
		int cnt = rnd->uniform.nextNonStd<int>(1, 5);
		for (int i = 0; i < cnt; ++i) {
			int selected_car2 = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
			int selected_pos2 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car2].size() - 1);
			int selected_car = selected_car2;
			while (selected_car == selected_car2)
				selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
			int selected_pos = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);
			int customer = offspring.variable().m_cus_order[selected_car][selected_pos];
			offspring.variable().m_cus_order[selected_car][selected_pos] = offspring.variable().m_cus_order[selected_car2][selected_pos2];
			offspring.variable().m_cus_order[selected_car2][selected_pos2] = customer;
			if (offspring.variable().m_cus_order[selected_car][0] != DVRP_CAST(pro)->get_depot() || offspring.variable().m_cus_order[selected_car2][0] != DVRP_CAST(pro)->get_depot()) {
				std::cout << "ls_rMulti() error" << "\n";
			}
		}
		DVRP_CAST(pro)->construct_path(offspring);

	}

	void DCNSGAII_LS_pop::ls_2opt(IndividualType &offspring, Problem *pro, Random *rnd)
	{
		int selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		while (offspring.variable().m_cus_order[selected_car].size() <= 3) {
			selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		}
		//随机选1个顾客
		int selected_pos1 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);

		int selected_pos2 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);
		while (selected_pos1 == selected_pos2 || abs(selected_pos1 - selected_pos2) > 4)
			selected_pos2 = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);

		if (selected_pos1 > selected_pos2) {
			selected_pos1 = selected_pos1 + selected_pos2;
			selected_pos2 = selected_pos1 - selected_pos2;
			selected_pos1 = selected_pos1 - selected_pos2;
		}

		std::vector<size_t> seq_1, seq_2, seq_3;

		for (int i = 1; i < selected_pos1; ++i) {
			seq_1.push_back(offspring.variable().m_cus_order[selected_car][i]);
		}
		for (int i = selected_pos2; i >= selected_pos1; --i) {
			seq_2.push_back(offspring.variable().m_cus_order[selected_car][i]);
		}
		for (int i = selected_pos2 + 1; i < offspring.variable().m_cus_order[selected_car].size() - 1; ++i) {
			seq_3.push_back(offspring.variable().m_cus_order[selected_car][i]);
		}
		offspring.variable().m_cus_order[selected_car].clear();
		offspring.variable().m_cus_order[selected_car].push_back(DVRP_CAST(pro)->get_depot());

		if (seq_1.size() != 0) {
			for (auto &i : seq_1)
				offspring.variable().m_cus_order[selected_car].push_back(i);
		}
		if (seq_2.size() != 0) {
			for (auto &i : seq_2)
				offspring.variable().m_cus_order[selected_car].push_back(i);
		}
		if (seq_3.size() != 0) {
			for (auto &i : seq_3)
				offspring.variable().m_cus_order[selected_car].push_back(i);
		}
		offspring.variable().m_cus_order[selected_car].push_back(DVRP_CAST(pro)->get_depot());
		if (offspring.variable().m_cus_order[selected_car][0] != DVRP_CAST(pro)->get_depot()) {
			std::cout << "ls_2opt() error" << "\n";
		}
		DVRP_CAST(pro)->construct_path(offspring);
	}

	void DCNSGAII_LS_pop::ls_maxWTlength(IndividualType &offspring, Problem *pro)
	{
		if (offspring.objective().size() == 1) {
			if (offspring.constraint()[3] == 0) //wait time ==0
				return;
		}
		if (offspring.objective().size() == 3) {
			if (offspring.objective()[2] == 0) //wait time ==0
				return;
		}

		int max_wait_car = 0;
		Real max_wait_time = 0.0;
		int max_wait_pos = 0;
		for (int j = 0; j < offspring.variable().m_hi.wait_time.size(); ++j) {
			for (int k = 0; k < offspring.variable().m_hi.wait_time[j].size() - 1; ++k) {
				if (max_wait_time < offspring.variable().m_hi.wait_time[j][k]) {
					max_wait_time = offspring.variable().m_hi.wait_time[j][k];
					max_wait_car = j;
					max_wait_pos = k + 1;
				}
			}
		}

		int max_length_car = 0;
		Real max_length = 0.0;
		for (int j = 0; j < offspring.variable().m_hi.route_length.size(); ++j) {
			if (max_length <= offspring.variable().m_hi.route_length[j]) {
				max_length = offspring.variable().m_hi.route_length[j];
				max_length_car = j;
			}
		}
		Real max_cost_time = 0.0;
		int max_cost_pos = 0;
		for (int i = 0; i < offspring.variable().m_hi.arrive_time[max_length_car].size() - 1; ++i) {
			if (max_cost_time <= offspring.variable().m_hi.arrive_time[max_length_car][i + 1] - offspring.variable().m_hi.arrive_time[max_length_car][i]) {
				max_cost_time = offspring.variable().m_hi.arrive_time[max_length_car][i + 1] - offspring.variable().m_hi.arrive_time[max_length_car][i];
				max_cost_pos = i + 1;
			}
		}
		if (max_length_car == max_wait_car && max_cost_pos == max_wait_pos)
			return;

		if (max_wait_time != 0) {
			int max_cost_cus = offspring.variable().m_cus_order[max_length_car][max_cost_pos];
			offspring.variable().m_cus_order[max_length_car].erase(offspring.variable().m_cus_order[max_length_car].begin() + max_cost_pos);
			offspring.variable().m_cus_order[max_wait_car].insert(offspring.variable().m_cus_order[max_wait_car].begin() + max_wait_pos, max_cost_cus);
			DVRP_CAST(pro)->construct_path(offspring);
			if (offspring.variable().m_cus_order[max_length_car][0] != 175 || offspring.variable().m_cus_order[max_wait_car][0] != DVRP_CAST(pro)->get_depot()) {
				std::cout << "ls_maxWTlength() error" << "\n";
			}
		}
	}

	void DCNSGAII_LS_pop::ls_maxLength(IndividualType &offspring, Problem *pro, Random *rnd) {
		if (offspring.objective()[0] == 0) //wait time ==0
			return;


		int max_length_car = 0;
		Real max_length = 0.0;
		for (int j = 0; j < offspring.variable().m_hi.route_length.size(); ++j) {
			if (max_length <= offspring.variable().m_hi.route_length[j]) {
				max_length = offspring.variable().m_hi.route_length[j];
				max_length_car = j;
			}
		}
		Real max_cost_time = 0.0;
		int max_cost_pos = 0;
		for (int i = 0; i < offspring.variable().m_hi.arrive_time[max_length_car].size() - 1; ++i) {
			if (max_cost_time <= offspring.variable().m_hi.arrive_time[max_length_car][i + 1] - offspring.variable().m_hi.arrive_time[max_length_car][i]) {
				max_cost_time = offspring.variable().m_hi.arrive_time[max_length_car][i + 1] - offspring.variable().m_hi.arrive_time[max_length_car][i];
				max_cost_pos = i + 1;
			}
		}

		int selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		while (selected_car == max_length_car)
			selected_car = rnd->uniform.nextNonStd<int>(0, offspring.variable().m_cus_order.size());
		int selected_pos = rnd->uniform.nextNonStd<int>(1, offspring.variable().m_cus_order[selected_car].size() - 1);

		if (max_length != 0) {
			int max_length_cus = offspring.variable().m_cus_order[max_length_car][max_cost_pos];
			offspring.variable().m_cus_order[max_length_car].erase(offspring.variable().m_cus_order[max_length_car].begin() + max_cost_pos);
			offspring.variable().m_cus_order[selected_car].insert(offspring.variable().m_cus_order[selected_car].begin() + selected_pos, max_length_cus);
			DVRP_CAST(pro)->construct_path(offspring);
		}
	}

	bool DCNSGAII_LS_pop::check_order(IndividualType &offspring, Problem *pro)
	{
		int num_order = 0;
		for (auto &i : offspring.variable().m_cus_order) {
			num_order += (int(i.size()) - 2);
		}
		if (num_order == DVRP_CAST(pro)->get_current_cus_order().size())
			return true;
		else if (num_order < DVRP_CAST(pro)->get_current_cus_order().size()) {
			std::cout << "order is not enough" << "\n";
			return false;
		}
		else {
			std::cout << "order is repeat" << "\n";
			return false;
		}
	}

	//DCNSGAII_LS::DCNSGAII_LS(param_map & v) : algorithm(v.at("algorithm name")), m_pop(v.at("population size")) {
	//	DVRP_CAST(pro)->set_eval_monitor_flag(true);
	//}

	void DCNSGAII_LS::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;;
		m_pop.reset();
		m_pop_size = v.get<int>("population size");
		//#ifdef OFEC_DEMO
		//		std::cout << "pop[i]" << "\t" << "rank\troute_length\twait_time\tdelay_time\tover_load\tover_time\tnum_car" << std::endl;
		//#endif

#ifdef OFEC_DEMO
		//std::cout << "pop[i]" << "\t" << "rank\troute_length\twait_time\tdelay_time\tover_load\tover_time\tnum_car" << std::endl;
		vector<vector<Solution<DVRP::routes> *>> pops(1);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			pops[0].emplace_back(&m_pop->at(i));
			//rank_.push_back(m_pop->at(i).fitness());
		}
		dynamic_cast<ofec_demo::buffer_DVRP *>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops, 0);
		dynamic_cast<ofec_demo::buffer_DVRP *>(ofec_demo::msp_buffer.get())->get_initial_pop(pops);
#endif
	}

	void DCNSGAII_LS::run_() {
		m_pop.reset(new DCNSGAII_LS_pop(m_pop_size, m_problem.get()));
		m_pop->initialize(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
		m_pop->initializeAfterEvaluation();

		DCMOEA_ind<Solution<DVRP::routes>> best_ind = m_pop->at(0);
		DCMOEA_ind<Solution<DVRP::routes>> ini_best_ind = m_pop->at(0);
		m_pop->sort();
		Real t_route_length = 0;
		Real t_wait_time = 0;
		Real t_delay_time = 0;
		Real t_over_load = 0;
		Real t_over_time = 0;
		std::string initialPOP_Fname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/InitialPOP" + to_string(m_pop->size()) + "_order" + to_string(DVRP_CAST(m_problem.get())->get_current_cus_order().size()) + "_" + m_pop->m_algName_ + "_" + m_pop->m_proName_ + "_" + m_pop->m_initialName_ + ".txt";
		std::ofstream pop_(initialPOP_Fname);
		std::stringstream pop_str;
		for (size_t i = 0; i < m_pop->size(); i++) {
			if (DVRP_CAST(m_problem.get())->numberObjectives() == 3) {
				t_route_length = m_pop->at(i).objective()[0];
				t_delay_time = m_pop->at(i).objective()[1];
				t_wait_time = m_pop->at(i).objective()[2];
				t_over_load = m_pop->at(i).constraint()[0];
				t_over_time = m_pop->at(i).constraint()[1];
			}
			if (DVRP_CAST(m_problem.get())->numberObjectives() == 1) {
				t_route_length = m_pop->at(i).objective()[0];
				t_wait_time = m_pop->at(i).constraint()[3];
				t_delay_time = m_pop->at(i).constraint()[2];
				t_over_load = m_pop->at(i).constraint()[0];
				t_over_time = m_pop->at(i).constraint()[1];
			}
			int t_num_car = m_pop->at(i).variable().m_cus_order.size();
			pop_str << to_string(i) << "\t" << to_string(m_pop->at(i).fitness()) << "\t" << t_route_length << "\t" << t_wait_time << "\t" << t_delay_time << "\t" << t_over_load << "\t" << t_over_time << "\t" << t_num_car << std::endl;
		}
		pop_ << pop_str.str();
		pop_.close();

		int iter = 0;
		while (!terminating()) {
			//auto domi_flag = m_pop->Pareto_compare(&best_ind, &ini_best_ind);
			if (iter >= 10000) {//10000  //&& domi_flag==Dominance::Dominant
				std::cout << "Online!" << std::endl;
				m_pop->Online(m_problem.get(), this, m_random.get());
				std::string Online_resFname;
				Online_resFname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/Online_result" + to_string(m_pop->size()) + "_order" + to_string(DVRP_CAST(m_problem.get())->get_current_cus_order().size()) + "_" + m_pop->m_algName_ + "_" + m_pop->m_proName_ + "_" + m_pop->m_initialName_ + ".txt";
				std::ofstream Online_result(Online_resFname);
				std::stringstream Online_result_;
				if (DVRP_CAST(m_problem.get())->numberObjectives() == 3) {
					Online_result_ << m_pop->get_bestInd().objective()[0] << "\t" << m_pop->get_bestInd().objective()[2] << "\t" << m_pop->get_bestInd().objective()[1] << "\t" << m_pop->get_bestInd().constraint()[0] << "\t" << m_pop->get_bestInd().constraint()[1] << "\t" << std::to_string(DVRP_CAST(m_problem.get())->get_future_cus_order().size()) << "\t" << m_pop->get_bestInd().variable().m_cus_order.size() << std::endl;
					std::cout << m_pop->get_bestInd().objective()[0] << "\t" << m_pop->get_bestInd().objective()[2] << "\t" << m_pop->get_bestInd().objective()[1] << "\t" << m_pop->get_bestInd().constraint()[0] << "\t" << m_pop->get_bestInd().constraint()[1] << "\t" << std::to_string(DVRP_CAST(m_problem.get())->get_future_cus_order().size()) << "\t" << m_pop->get_bestInd().variable().m_cus_order.size() << std::endl;
				}
				if (DVRP_CAST(m_problem.get())->numberObjectives() == 1) {
					Online_result_ << m_pop->get_bestInd().objective()[0] << "\t" << m_pop->get_bestInd().constraint()[3] << "\t" << m_pop->get_bestInd().constraint()[2] << "\t" << m_pop->get_bestInd().constraint()[0] << "\t" << m_pop->get_bestInd().constraint()[1] << "\t" << std::to_string(DVRP_CAST(m_problem.get())->get_future_cus_order().size()) << "\t" << m_pop->get_bestInd().variable().m_cus_order.size() << std::endl;
					std::cout << m_pop->get_bestInd().objective()[0] << "\t" << m_pop->get_bestInd().constraint()[3] << "\t" << m_pop->get_bestInd().constraint()[2] << "\t" << m_pop->get_bestInd().constraint()[0] << "\t" << m_pop->get_bestInd().constraint()[1] << "\t" << std::to_string(DVRP_CAST(m_problem.get())->get_future_cus_order().size()) << "\t" << m_pop->get_bestInd().variable().m_cus_order.size() << std::endl;
				}
				Online_result << Online_result_.str();
				Online_result.close();
				std::cout << "Online Complete!" << "\n";
				system("pause");
				//m_isOnline = true;
			}
			m_pop->evolve(m_problem.get(), this, m_random.get());
			m_pop->sort();
#ifdef OFEC_DEMO
			/*vector<vector<Solution<DVRP::routes>*>> pops(1);
			for (size_t i = 0; i < m_pop->size(); ++i)
				pops[0].emplace_back(&m_pop->at(i));
			dynamic_cast<ofec_demo::buffer_DVRP*>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops);*/

			/* Begin terminal window output */
			size_t evals = DVRP_CAST(pro)->evaluations();
			Real min_obj = (std::numeric_limits<Real>::max)();
			size_t idx_best;
			Real temp_obj;


			std::vector<std::pair<double, double>> objective_min_max(DVRP_CAST(pro)->numberObjectives());
			std::vector<std::pair<double, double>> constraint_min_max(DVRP_CAST(pro)->num_constraints());
			std::vector<std::pair<double, double>> vioObj_min_max(m_pop->get_num_vio_obj());

			for (int i = 0; i < DVRP_CAST(pro)->numberObjectives(); ++i) {
				std::pair<double, double> min_max;
				min_max.first = (std::numeric_limits<Real>::max)();
				min_max.second = (std::numeric_limits<Real>::min)();
				for (size_t j = 0; j < m_pop->size(); j++) {
					if (min_max.first >= m_pop->at(j].objective()[i])
						min_max.first = m_pop->at(j].objective()[i];
					if (min_max.second <= m_pop->at(j].objective()[i])
						min_max.second = m_pop->at(j].objective()[i];
				}
				objective_min_max[i] = min_max;
			}
			for (int i = 0; i < DVRP_CAST(pro)->num_constraints(); ++i) {
				std::pair<double, double> min_max;
				min_max.first = (std::numeric_limits<Real>::max)();
				min_max.second = (std::numeric_limits<Real>::min)();
				for (size_t j = 0; j < m_pop->size(); j++) {
					if (min_max.first >= m_pop->at(j].constraint()[i])
						min_max.first = m_pop->at(j].constraint()[i];
					if (min_max.second <= m_pop->at(j].constraint()[i])
						min_max.second = m_pop->at(j].constraint()[i];
				}
				constraint_min_max[i] = min_max;
			}
			for (int i = 0; i < m_pop->get_num_vio_obj(); ++i) {
				std::pair<double, double> min_max;
				min_max.first = (std::numeric_limits<Real>::max)();
				min_max.second = (std::numeric_limits<Real>::min)();
				for (size_t j = 0; j < m_pop->size(); j++) {
					if (min_max.first >= m_pop->at(j].get_vio_obj()[i])
						min_max.first = m_pop->at(j].get_vio_obj()[i];
					if (min_max.second <= m_pop->at(j].get_vio_obj()[i])
						min_max.second = m_pop->at(j].get_vio_obj()[i];
				}
				vioObj_min_max[i] = min_max;
			}

			Real route_length = 0;
			Real wait_time = 0;
			Real delay_time = 0;
			Real over_load = 0;
			Real over_time = 0;
			std::string result_Fname;
			result_Fname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/result_pop" + to_string(m_pop->size()) + "_order" + to_string(DVRP_CAST(pro)->get_current_cus_order().size()) + "_" + m_pop->m_algName_ + "_" + m_pop->m_proName_ + "_" + m_pop->m_initialName_ + ".txt";
			std::ofstream result(result_Fname, ios::app);
			std::stringstream result_;
			result_ << iter << "\n";

			Real delay_timeTmp = (std::numeric_limits<Real>::max)();
			for (size_t i = 0; i < m_pop->size(); i++) {
				if (DVRP_CAST(pro)->numberObjectives() == 3) {
					route_length = m_pop->at(i).objective()[0];
					delay_time = m_pop->at(i).objective()[1];
					wait_time = m_pop->at(i).objective()[2];
					over_load = m_pop->at(i).constraint()[0];
					over_time = m_pop->at(i).constraint()[1];
				}
				if (DVRP_CAST(pro)->numberObjectives() == 1) {
					route_length = m_pop->at(i).objective()[0];
					wait_time = m_pop->at(i).constraint()[3];
					delay_time = m_pop->at(i).constraint()[2];
					over_load = m_pop->at(i).constraint()[0];
					over_time = m_pop->at(i).constraint()[1];
				}
				int num_car = m_pop->at(i).variable().m_cus_order.size();
				//std::cout << "pop[" + to_string(i) << "]\t" << to_string(m_pop->at(i).fitness()) << "\t" << route_length << "\t" << wait_time << "\t" << delay_time << "\t" << over_load << "\t" << over_time << "\t" << num_car << std::endl;
				result_ << to_string(i) << "\t" << to_string(m_pop->at(i).fitness()) << "\t" << route_length << "\t" << wait_time << "\t" << delay_time << "\t" << over_load << "\t" << over_time << "\t" << num_car << std::endl;



				if (m_pop->at(i).fitness() == 0) {
					if (delay_timeTmp >= delay_time) {
						delay_timeTmp = delay_time;
						idx_best = i;
					}
					/*temp_obj = 0;
					for (int j = 0; j < DVRP_CAST(pro)->numberObjectives(); ++j) {
						if (objective_min_max[j].second != 0 && objective_min_max[j].first != objective_min_max[j].second) {
							temp_obj += ((m_pop->at(i).objective()[j] - objective_min_max[j].first) / (objective_min_max[j].second - objective_min_max[j].first));
						}
						else
						{
							temp_obj += m_pop->at(i).objective()[j];
						}
					}
					for (int j = 0; j < DVRP_CAST(pro)->num_constraints(); ++j) {
						if (constraint_min_max[j].second != 0 && constraint_min_max[j].first != constraint_min_max[j].second)
							temp_obj += ((m_pop->at(i).constraint()[j] - constraint_min_max[j].first) / (constraint_min_max[j].second - constraint_min_max[j].first));
						else
							temp_obj += m_pop->at(i).constraint()[j];
					}
					for (int j = 0; j < m_pop->get_num_vio_obj(); ++j) {
						if (vioObj_min_max[j].second != 0 && vioObj_min_max[j].first != vioObj_min_max[j].second)
							temp_obj += ((m_pop->at(i).get_vio_obj()[j] - vioObj_min_max[j].first) / (vioObj_min_max[j].second - vioObj_min_max[j].first));
						else
							temp_obj += m_pop->at(i).get_vio_obj()[j];
					}
					if (min_obj >= temp_obj) {
						min_obj = temp_obj;
						idx_best = i;
					}*/
				}
			}
			result << result_.str();
			result.close();
			Real _route_length = 0;
			Real _wait_time = 0;
			Real _delay_time = 0;
			Real _over_load = 0;
			Real _over_time = 0;

			if (DVRP_CAST(pro)->numberObjectives() == 3) {
				_route_length = m_pop->at(idx_best].objective()[0];
				_delay_time = m_pop->at(idx_best].objective()[1];
				_wait_time = m_pop->at(idx_best].objective()[2];
				_over_load = m_pop->at(idx_best].constraint()[0];
				_over_time = m_pop->at(idx_best].constraint()[1];
			}
			if (DVRP_CAST(pro)->numberObjectives() == 1) {
				_route_length = m_pop->at(idx_best].objective()[0];
				_wait_time = m_pop->at(idx_best].constraint()[3];
				_delay_time = m_pop->at(idx_best].constraint()[2];
				_over_load = m_pop->at(idx_best].constraint()[0];
				_over_time = m_pop->at(idx_best].constraint()[1];
			}
			std::cout << iter << "\t" << _route_length << "\t" << _wait_time << "\t" << _delay_time << "\t" << _over_load << "\t" << _over_time << std::endl;

			std::string progress_Fname;
			progress_Fname = static_cast<std::string>(g_working_dir) + "/instance/algorithm/realworld/DVRP/Compare_data/progress" + to_string(m_pop->size()) + "_order" + to_string(DVRP_CAST(pro)->get_current_cus_order().size()) + "_" + m_pop->m_algName_ + "_" + m_pop->m_proName_ + "_" + m_pop->m_initialName_ + ".txt";
			std::ofstream progress(progress_Fname, ios::app);
			std::stringstream progress_;
			progress_ << iter << "\t" << _route_length << "\t" << _wait_time << "\t" << _delay_time << "\t" << _over_load << "\t" << _over_time << std::endl;
			progress << progress_.str();
			progress.close();

			best_ind = m_pop->at(idx_best];
			if (iter == 0)
				ini_best_ind = best_ind;
			m_pop->set_bestInd(best_ind);
			/*	std::cout << "Online" << std::endl;
				m_pop->Online();
				std::cout << m_pop->get_bestInd().objective()[0] << "\t" << m_pop->get_bestInd().constraint()[3] << "\t" << m_pop->get_bestInd().constraint()[2] << "\t" << m_pop->get_bestInd().constraint()[0] << "\t" << m_pop->get_bestInd().constraint()[1] << "\t" << std::to_string(DVRP_CAST(pro)->get_future_cus_order().size()) << std::endl;
	*/
	/* End terminal window output */
	//system("pause");
			vector<vector<Solution<DVRP::routes> *>> pops(1);
			for (size_t i = 0; i < m_pop->size(); ++i)
				pops[0].emplace_back(&m_pop->at(i));
			dynamic_cast<ofec_demo::buffer_DVRP *>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops, idx_best);
#endif
			iter++;
		}
		//m_pop->sort();
	}

	void DCNSGAII_LS::record() {
		m_pop->sort();
		size_t evals = m_evaluations;
		Real min_obj = (std::numeric_limits<Real>::max)();
		size_t idx_best;
		Real temp_obj;
		for (size_t i = 0; i < m_pop->size(); i++) {
			if (m_pop->at(i).fitness() == 0) {
				temp_obj = 0;
				for (int j = 0; j < DVRP_CAST(m_problem.get())->numberObjectives(); ++j) {
					temp_obj += m_pop->at(i).objective()[j];
				}
				for (int j = 0; j < DVRP_CAST(m_problem.get())->numberConstraints(); ++j) {
					temp_obj += m_pop->at(i).constraint()[j];
				}
				for (int j = 0; j < m_pop->get_num_vio_obj(); ++j) {
					temp_obj += m_pop->at(i).get_vio_obj()[j];
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

		if (DVRP_CAST(m_problem.get())->numberObjectives() == 3) {
			route_length = m_pop->at(idx_best).objective()[0];
			delay_time = m_pop->at(idx_best).objective()[1];
			wait_time = m_pop->at(idx_best).objective()[2];
			over_load = m_pop->at(idx_best).constraint()[0];
			over_time = m_pop->at(idx_best).constraint()[1];
		}
		if (DVRP_CAST(m_problem.get())->numberObjectives() == 1) {
			route_length = m_pop->at(idx_best).objective()[0];
			wait_time = m_pop->at(idx_best).constraint()[3];
			delay_time = m_pop->at(idx_best).constraint()[2];
			over_load = m_pop->at(idx_best).constraint()[0];
			over_time = m_pop->at(idx_best).constraint()[1];
		}
		std::vector<Real> entry = { (Real)evals, route_length, wait_time, delay_time, over_load, over_time };
		dynamic_cast<RecordVectorReal*>(m_record.get())->record(this, entry);
	}
}