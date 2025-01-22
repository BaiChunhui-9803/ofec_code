#include "dynamic_vrp.h"
#include "../../../../core/algorithm/solution.h"
#include <functional>
using namespace std::placeholders;

namespace ofec {
	namespace DVRP {
		/*std::vector<std::vector<std::multimap<int, std::vector<size_t>>>> dynamic_vrp::m_shortest_path_all_day;
		std::vector<std::vector<std::multimap<int, Real>>> dynamic_vrp::m_cost_time_all_day;
		std::vector<std::vector<std::multimap<int, Real>>> dynamic_vrp::m_shortest_dis_all_day;*/

		//dynamic_vrp::dynamic_vrp(const ParameterMap &v) : dynamic_vrp(v.at("problem name"), v.at("numCus"), v.at("numNewCus"), v.at("numObj"), v.at("dataFile1"), v.at("dataFile2"),v.at("dataFile3")) {

		//}
		//dynamic_vrp::dynamic_vrp(const std::string& pro_name, size_t num_cus, size_t num_new_cus, size_t size_obj, const std::string &road_file_name, const std::string &velocity_file_name_Offline, const std::string &velocity_file_name_Online) :
		//	problem(pro_name, 0, size_obj), m_current_order_num(num_cus), m_future_order_num(num_new_cus), m_road_file_name(road_file_name), m_velocity_file_name_Offline(velocity_file_name_Offline), m_velocity_file_name_Online(velocity_file_name_Online) {
		//}

		void dynamic_vrp::initialize_() {
			Problem::initialize_();
			m_road_net.read_data(m_road_file_name, m_velocity_file_name_Offline, m_velocity_file_name_Online);
			m_road_net.set_velocity(false);
			m_road_net.constrct_net(m_depot);
			//readOnlineData();
			m_road_net.node_distance();
			m_road_net.select_customers(m_depot, m_current_order_num, m_future_order_num, m_depart_time, m_dynamic_order_endtime, m_random.get());//(depot, static_order, dynamic_order, depart_time, dynamic_endtime)
			m_nodes_size = m_road_net.get_road_net().size();

			m_vehicle_property.push_back({ 1 });//set vehicle type, capacity = 1
			m_vehicle_property.push_back({ 2 });
			m_vehicle_property.push_back({ 3 });
			m_constraint.assign(4, Constraint::Equality);

			//get customer sequence
			for (size_t i = 0; i < m_road_net.get_customers().size(); i++) {
				if (m_road_net.get_road_net()[m_road_net.get_customers()[i].first]->type == node_type::current_customer) {
					m_current_cus_order.push_back(m_road_net.get_customers()[i].first);
				}
			}
			for (size_t i = 0; i < m_road_net.get_customers().size(); i++) {
				if (m_road_net.get_road_net()[m_road_net.get_customers()[i].first]->type == node_type::future_customer) {
					m_future_cus_order.push_back(m_road_net.get_customers()[i]);
				}
			}

			//set_online();

			//std::function<void(const routes &, Real, size_t, size_t)> obj0 = std::bind(&dynamic_vrp::used_car, this, _1, _2, _3, _4);
			std::function<void(routes &, Real, size_t, size_t)> obj1 = std::bind(&dynamic_vrp::length, this, _1, _2, _3, _4);
			std::function<void(routes &, Real, size_t, size_t)> obj2 = std::bind(&dynamic_vrp::delay_time, this, _1, _2, _3, _4);
			std::function<void(routes &, Real, size_t, size_t)> obj3 = std::bind(&dynamic_vrp::arrive_time, this, _1, _2, _3, _4);
			std::function<void(routes &, Real, size_t, size_t)> obj4 = std::bind(&dynamic_vrp::wait_time, this, _1, _2, _3, _4);
			//std::function<void(const routes &, Real, size_t, size_t)> constr0 = std::bind(&dynamic_vrp::over_load, this, _1, _2, _3, _4);
			//std::function<void(const routes &, Real, size_t, size_t)> constr1 = std::bind(&dynamic_vrp::over_depot_end_time, this, _1, _2, _3, _4);

			m_objective.push_back(obj1);
			m_objective.push_back(obj2);
			m_objective.push_back(obj3);
			m_objective.push_back(obj4);

			m_cost_time_all_day.clear();
			m_shortest_dis_all_day.clear();
			m_shortest_path_all_day.clear();

			m_cost_time_all_day.resize(24 * 60);
			for(auto&i:m_cost_time_all_day){
				i.resize(m_nodes_size);
			}

			m_shortest_dis_all_day.resize(24 * 60);
			for (auto&i : m_shortest_dis_all_day) {
				i.resize(m_nodes_size);
			}

			m_shortest_path_all_day.resize(24 * 60);
			for (auto&i : m_shortest_path_all_day) {
				i.resize(m_nodes_size);
			}

			bool calculated = false;
			calcul_shorest_dis_24h(calculated);
		}

		void dynamic_vrp::set_Online(bool is_online) {
			road_net road_net_t;
			road_net_t.read_data(m_road_file_name, m_velocity_file_name_Offline, m_velocity_file_name_Online);
			road_net_t.set_velocity(is_online);
			road_net_t.constrct_net(m_depot);
			for (int i = 0; i < m_road_net.get_road_net().size(); ++i) {
				m_road_net.get_road_net()[i]->datum_t.v = road_net_t.get_road_net()[i]->datum_t.v;
			}
			road_net_t.~road_net();

			m_road_net.node_distance();
		}

		std::vector<bool> dynamic_vrp::is_valid(SolutionBase &s) const {
			const routes &x = dynamic_cast<Solution<routes> &>(s).variable();
			std::vector<bool> m_is_initial_valid(x.m_cus_order.size());
			Real over_load = 0;
			double deliveryWeight = 0;
			for (size_t i = 0; i < x.m_cus_order.size(); ++i) {
				double deliveryWeight = 0;
				for (size_t j = 0; j < x.m_cus_order[i].size(); ++j) {
					deliveryWeight = deliveryWeight + m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.delivery_weight;
				}
				if (deliveryWeight > m_vehicle_property[x.m_vehicle_type[i]].capacity) {
					m_is_initial_valid[i] = true;
				}
				else m_is_initial_valid[i] = false;
			}
			return m_is_initial_valid;

			/*for (size_t i = 0; i < x.m_cus_order.size(); ++i) {
				over_load = 0;
				deliveryWeight = 0;
				for (size_t j = 1; j < x.m_cus_order[i].size(); ++j) {
					if (m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.is_served == false)
						deliveryWeight += m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.delivery_weight;
				}
				if (deliveryWeight > m_vehicle_property[x.m_vehicle_type[i]].capacity) {
					over_load += deliveryWeight - m_vehicle_property[x.m_vehicle_type[i]].capacity;
				}
				for (size_t j = 1; j < x.m_cus_order[i].size(); ++j) {
					if (m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.is_served == false)
						deliveryWeight = deliveryWeight - m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.delivery_weight + m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.pick_weight;
					else
						deliveryWeight = deliveryWeight - 0 + m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.pick_weight;

					if (deliveryWeight > m_vehicle_property[x.m_vehicle_type[i]].capacity) {
						over_load += deliveryWeight - m_vehicle_property[x.m_vehicle_type[i]].capacity;
					}
				}
				if (over_load > 0)
					m_is_initial_valid[i] = true;
				else
					m_is_initial_valid[i] = false;
			}
			return m_is_initial_valid;*/
		}

		void dynamic_vrp::initializeSolution(SolutionBase &s, Random *rnd) const{
			routes &x = dynamic_cast<Solution<routes> &>(s).variable();
			auto cus_order_temp = m_current_cus_order;
			rnd->uniform.shuffle(cus_order_temp.begin(), cus_order_temp.end());
			x.m_cus_order.push_back(cus_order_temp);
			cus_order_temp.clear();
			x.m_vehicle_type.clear();
			x.m_vehicle_type.resize(x.m_cus_order.size());
			for (size_t i = 0; i < x.m_vehicle_type.size(); i++) {
				//size_t type = global::ms_global->m_uniform[caller::Problem]->next_non_standard<size_t>(0, m_vehicle_property.size());
				x.m_vehicle_type[i] = 2;
			}
			for (size_t i = 0; i < x.m_cus_order.size(); i++) {
				while (is_valid(s)[i]) {
					cus_order_temp.push_back(x.m_cus_order[i][x.m_cus_order[i].size() - 1]);
					x.m_cus_order[i].pop_back();
				}
				if (!cus_order_temp.empty()) {
					x.m_cus_order.push_back(cus_order_temp);
					cus_order_temp.clear();
					//size_t type = global::ms_global->m_uniform[caller::Problem]->next_non_standard<size_t>(0, m_vehicle_property.size());
					x.m_vehicle_type.push_back(2);
					/*if (m_road_net.get_road_net()[x.m_cus_order[x.m_cus_order.size() - 1][0]]->datum_t.delivery_weight > 2) {
						x.m_vehicle_type.push_back(2);
					}
					else
						x.m_vehicle_type.push_back(1);*/
				}
			}
			//to add sort for customers
			/*for (size_t i = 0; i < x.m_cus_order.size(); ++i) {
				sort_customers(x.m_cus_order[i]);
			}*/
			std::vector<size_t>::iterator it;
			for (auto &i : x.m_cus_order) {
				it = i.begin();
				i.insert(it, m_depot);
				it = i.end();
				i.insert(it, m_depot);
			}
			//get m_path
			construct_path(s);
		}
		void dynamic_vrp::sort_customers(std::vector<size_t> &customer_order) const
		{
			std::vector<Real> dis_from_nearest;
			std::vector<Real> dis_from_nearest_t;
			std::vector<size_t> cus_o_t;
			int nearest = m_depot;
			while (!customer_order.empty()) {
				for (size_t j = 0; j < customer_order.size(); ++j) {
					if (nearest == customer_order[j])
						continue;
					std::vector<size_t> path_from_nearest;
					auto it_path = m_shortest_path_all_day[static_cast<int>(m_depart_time) % 1439][nearest].find(customer_order[j]);
					if (it_path != m_shortest_path_all_day[static_cast<int>(m_depart_time) % 1439][nearest].end()) {
						path_from_nearest = it_path->second;
					}
					else
					{
						shortest_path(nearest, customer_order[j], path_from_nearest, m_depart_time);
					}
					
					Real dis = 0;
					for (size_t z = 0; z < path_from_nearest.size() - 1; ++z) {
						dis += m_road_net.get_node_dis()[path_from_nearest[z]][path_from_nearest[z + 1]];
					}
					dis_from_nearest.push_back(dis);
				}
				dis_from_nearest_t = dis_from_nearest;
				sort(dis_from_nearest_t.begin(), dis_from_nearest_t.end());
				auto it = find(dis_from_nearest.begin(), dis_from_nearest.end(), dis_from_nearest_t[0]);
				auto pos = std::distance(dis_from_nearest.begin(), it);
				nearest = customer_order[pos];
				cus_o_t.push_back(nearest);
				dis_from_nearest.clear();
				customer_order.erase(customer_order.begin() + pos);
			}
			customer_order = cus_o_t;
		}
		void dynamic_vrp::construct_path(SolutionBase & s) {
			routes &x = dynamic_cast<Solution<routes> &>(s).variable();
			x.m_path.clear();
			x.m_path_node_type.clear();
			x.m_path.resize(x.m_cus_order.size());
			x.m_path_node_type.resize(x.m_cus_order.size());
			for (size_t i = 0; i < x.m_cus_order.size(); i++) {
				int cnt = 0;
				std::vector<std::vector<size_t>> path_temp(x.m_cus_order[i].size() - 1);
				std::vector<std::vector<node_type>> path_node_type_temp(x.m_cus_order[i].size() - 1);
				//path_temp[0].push_back(m_depot);
				Real present_time = m_depart_time;
				for (size_t j = 0; j < x.m_cus_order[i].size() - 1; j++) {
					std::vector<size_t> Shortest_path;
					auto it_path = m_shortest_path_all_day[static_cast<int>(present_time) % 1439][x.m_cus_order[i][j]].find(x.m_cus_order[i][j + 1]);
					if (it_path != m_shortest_path_all_day[static_cast<int>(present_time) % 1439][x.m_cus_order[i][j]].end()) {
						Shortest_path = it_path->second;
					}
					else
					{
						shortest_path(x.m_cus_order[i][j], x.m_cus_order[i][j + 1], Shortest_path, present_time);//时间要变
					}
					
					if (m_road_net.get_road_net()[x.m_cus_order[i][j]]->type == node_type::current_customer || m_road_net.get_road_net()[x.m_cus_order[i][j]]->type == node_type::future_customer) {
						present_time = present_time + get_cost_time(present_time, x.m_cus_order[i][j], x.m_cus_order[i][j + 1]) + m_road_net.get_road_net()[x.m_cus_order[i][j]]->server_time;
					}
					if (m_road_net.get_road_net()[x.m_cus_order[i][j]]->type == node_type::depot) {
						present_time = present_time + get_cost_time(present_time, x.m_cus_order[i][j], x.m_cus_order[i][j + 1]);
					}
					auto tt = get_cost_time(present_time, x.m_cus_order[i][j], x.m_cus_order[i][j + 1]);
					Shortest_path.pop_back();
					path_temp[j] = Shortest_path;
					for (auto &m : path_temp[j]) {
						x.m_path[i].push_back(m);
						path_node_type_temp[j].push_back(node_type::general_node);
					}
					path_node_type_temp[j][0] = m_road_net.get_road_net()[path_temp[j][0]]->type;
					for (auto &m : path_node_type_temp[j]) {
						x.m_path_node_type[i].push_back(m);
					}
				}
				x.m_path[i].push_back(m_depot);
				x.m_path_node_type[i].push_back(m_road_net.get_road_net()[m_depot]->type);

			}
			/*std::ofstream f_path("./instance/problem/realworld/DVRP/data/OffLine/m_path.txt");
			for (size_t i = 0; i < x.m_path.size(); i++)
			{
				for (auto j : x.m_path[i])
					f_path << j << " ";
			}
			f_path.close();*/
		}

		void dynamic_vrp::construct_path(SolutionBase & s) const{
			routes &x = dynamic_cast<Solution<routes> &>(s).variable();
			x.m_path.clear();
			x.m_path_node_type.clear();
			x.m_path.resize(x.m_cus_order.size());
			x.m_path_node_type.resize(x.m_cus_order.size());
			for (size_t i = 0; i < x.m_cus_order.size(); i++) {
				int cnt = 0;
				std::vector<std::vector<size_t>> path_temp(x.m_cus_order[i].size() - 1);
				std::vector<std::vector<node_type>> path_node_type_temp(x.m_cus_order[i].size() - 1);
				//path_temp[0].push_back(m_depot);
				Real present_time = m_depart_time;
				for (size_t j = 0; j < x.m_cus_order[i].size() - 1; j++) {
					std::vector<size_t> Shortest_path;
					auto it_path = m_shortest_path_all_day[static_cast<int>(present_time) % 1439][x.m_cus_order[i][j]].find(x.m_cus_order[i][j + 1]);
					if (it_path != m_shortest_path_all_day[static_cast<int>(present_time) % 1439][x.m_cus_order[i][j]].end()) {
						Shortest_path = it_path->second;
					}
					else
					{
						shortest_path(x.m_cus_order[i][j], x.m_cus_order[i][j + 1], Shortest_path, present_time);//时间要变
					}
					Shortest_path.pop_back();
					path_temp[j] = Shortest_path;
					for (auto &m : path_temp[j]) {
						x.m_path[i].push_back(m);
						path_node_type_temp[j].push_back(node_type::general_node);
					}
					path_node_type_temp[j][0] = m_road_net.get_road_net()[path_temp[j][0]]->type;
					for (auto &m : path_node_type_temp[j]) {
						x.m_path_node_type[i].push_back(m);
					}
				}
				x.m_path[i].push_back(m_depot);
				x.m_path_node_type[i].push_back(m_road_net.get_road_net()[m_depot]->type);

			}
			/*std::ofstream f_path("./instance/problem/realworld/DVRP/data/m_path.txt");
			for (size_t i = 0; i < x.m_path.size(); i++)
			{
				for (auto j : x.m_path[i])
					f_path << j << " ";
			}
			f_path.close();*/
		}

		bool dynamic_vrp::same(const SolutionBase &s1, const SolutionBase &s2) const {
			const routes &x1 = dynamic_cast<const Solution<routes> &>(s1).variable();
			const routes &x2 = dynamic_cast<const Solution<routes> &>(s2).variable();
			if (x1.m_cus_order.size() != x2.m_cus_order.size()) return false;
			if (x1.m_cus_order != x2.m_cus_order) return false;
			if (x1.m_path != x2.m_path) return false;
			if (x1.m_vehicle_type != x2.m_vehicle_type) return false;
			else return true;
		}

		void dynamic_vrp::evaluate_(SolutionBase &s, bool effective) {
			routes &x = dynamic_cast<Solution<routes> &>(s).variable();
			std::vector<Real> &obj = dynamic_cast<Solution<routes>&>(s).objective();
			std::vector<Real> &con = dynamic_cast<Solution<routes>&>(s).constraint();

			x.m_hi.route_length.clear();
			x.m_hi.delay_time.clear();
			x.m_hi.arrive_time.clear();
			x.m_hi.wait_time.clear();

			x.m_hi.over_time.clear();
			x.m_hi.over_load.clear();
			x.m_hi.total_delivery_weight.clear();


			/*m_route_length.clear();
			m_delay_time.clear();
			m_travel_time.clear();
			m_wait_time.clear();
			m_over_time.clear();


			m_total_delivery_weight.clear();
			m_over_load.clear();*/

			x.m_hi.route_length.resize(x.m_cus_order.size());
			x.m_hi.delay_time.resize(x.m_cus_order.size());
			x.m_hi.arrive_time.resize(x.m_cus_order.size());
			x.m_hi.wait_time.resize(x.m_cus_order.size());

			x.m_hi.over_time.resize(x.m_cus_order.size());
			x.m_hi.total_delivery_weight.resize(x.m_cus_order.size());
			x.m_hi.over_load.resize(x.m_cus_order.size());

			visit_routes(x);

			if (obj.size() == 5) {
				used_car(x);
				obj[0] = x.m_hi.used_car;

				obj[1] = 0;
				for (auto &row : x.m_hi.route_length) {
					obj[1] += row;
				}

				obj[2] = 0;
				for (auto &row : x.m_hi.delay_time) {
					for (auto &col : row) {
						obj[2] += col;
					}
				}

				obj[3] = 0;
				for (auto &row : x.m_hi.arrive_time) { 
					obj[3] += (row[row.size() - 1] - m_depart_time);
				}

				obj[4] = 0;
				for (auto &row : x.m_hi.wait_time) {
					for (auto &col : row) {
						obj[4] += col;
					}
				}
			}

			if (obj.size() == 4) {
				used_car(x);
				obj[0] = x.m_hi.used_car;

				obj[1] = 0;
				for (auto &row : x.m_hi.route_length) {
					obj[1] += row;
				}

				obj[2] = 0;
				for (auto &row : x.m_hi.delay_time) {
					for (auto &col : row) {
						obj[2] += col;
					}
				}

				obj[3] = 0;
				for (auto &row : x.m_hi.wait_time) {
					for (auto &col : row) {
						obj[3] += col;
					}
				}
			}

			if (obj.size() == 3) {
				used_car(x);
				obj[0] = 0;
				for (auto &row : x.m_hi.route_length) {
					obj[0] += row;
				}

				obj[1] = 0;
				for (auto &row : x.m_hi.delay_time) {
					for (auto &col : row) {
						obj[1] += col;
					}
				}

				obj[2] = 0;
				for (auto &row : x.m_hi.wait_time) {
					for (auto &col : row) {
						obj[2] += col;
					}
				}
			}
			if (obj.size() == 2) {
				used_car(x);
				obj[0] = 0;
				for (auto &row : x.m_hi.route_length) {
					obj[0] += row;
				}

				obj[1] = 0;
				for (auto &row : x.m_hi.wait_time) {
					for (auto &col : row) {
						obj[1] += col;
					}
				}

				con[2] = 0;
				for (auto &row : x.m_hi.delay_time) {
					for (auto &col : row) {
						con[2] += col;
					}
				}
			}
			if (obj.size() == 1) {
				used_car(x);
				obj[0] = 0;
				for (auto &row : x.m_hi.route_length) {
					obj[0] += row;
				}

				con[3] = 0;
				for (auto &row : x.m_hi.wait_time) {
					for (auto &col : row) {
						con[3] += col;
					}
				}

				con[2] = 0;
				for (auto &row : x.m_hi.delay_time) {
					for (auto &col : row) {
						con[2] += col;
					}
				}
			}
			
			//constraints
			over_load(x);
			con[0] = 0;
			for (auto &row : x.m_hi.over_load) {
				con[0] += row;
			}
			over_depot_end_time(x);
			con[1] = 0;
			for (auto &row : x.m_hi.over_time) {
				con[1] += row;
			}
		}

		Real dynamic_vrp::variableDistance(const SolutionBase &s1, const SolutionBase &s2) const {
			std::vector<std::vector<bool>> s1_edge(m_nodes_size);
			const routes &x1 = dynamic_cast<const Solution<routes> &>(s1).variable();
			const routes &x2 = dynamic_cast<const Solution<routes> &>(s2).variable();
			for (auto &row : s1_edge)
				row.assign(m_nodes_size, false);
			for (int i = 0; i < m_nodes_size; ++i) {
				for (int m = 0; m < x1.m_cus_order.size(); m++) {
					for (int n = 0; n < x1.m_cus_order[m].size() - 1; n++) {
						s1_edge[x1.m_cus_order[m][n]][x1.m_cus_order[m][n + 1]] = true;
					}
				}
			}
			Real dis = 0;
			for (int m = 0; m < x2.m_cus_order.size(); m++) {
				for (int n = 0; n < x2.m_cus_order[m].size() - 1; n++) {
					if (!s1_edge[x2.m_cus_order[m][n]][x2.m_cus_order[m][n + 1]]) dis += 1;
				}
			}
			return dis;
		}
		Real dynamic_vrp::variableDistance(const VariableBase &s1, const VariableBase &s2) const {
			std::vector<std::vector<bool>> s1_edge(m_nodes_size);
			const routes &x1 = dynamic_cast<const routes&>(s1);
			const routes &x2 = dynamic_cast<const routes&>(s2);
			for (auto &row : s1_edge)
				row.assign(m_nodes_size, false);
			for (int i = 0; i < m_nodes_size; ++i) {
				for (int m = 0; m < x1.m_cus_order.size(); m++) {
					for (int n = 0; n < x1.m_cus_order[m].size() - 1; n++) {
						s1_edge[x1.m_cus_order[m][n]][x1.m_cus_order[m][n + 1]] = true;
					}
				}
			}
			Real dis = 0;
			for (int m = 0; m < x2.m_cus_order.size(); m++) {
				for (int n = 0; n < x2.m_cus_order[m].size() - 1; n++) {
					if (!s1_edge[x2.m_cus_order[m][n]][x2.m_cus_order[m][n + 1]]) dis += 1;
				}
			}
			return dis;
		}

		void dynamic_vrp::shortest_path(size_t start, size_t end, std::vector<size_t> &Shortest_path, Real present_time){
			//A star
			if (start == end) {
				std::cout << "The start and end of shotest_path is equal" << std::endl;
				system("pause");
			}
			else {
				std::vector<size_t> open_list;
				std::vector<size_t> closeed_list;
				open_list.push_back(start);
				std::vector<int> came_from(m_nodes_size);
				/*for (size_t i = 0; i < m_road_net.get_road_net()[start]->next.size(); ++i) {
					open_list.push_back(m_road_net.get_road_net()[start]->next[i]->id);
				}*/
				std::vector<double> fScore(m_nodes_size);
				std::vector<double> gScore(m_nodes_size);
				for (auto &i : gScore)
					i = std::numeric_limits<double>::max();//km(Infinity)
				gScore[start] = 0;
				for (auto &i : fScore)
					i = std::numeric_limits<double>::max();// Infinity
				fScore[start] = m_road_net.get_node_dis()[start][end] / m_road_net.get_road_net()[start]->datum_t.v[static_cast<int>(present_time) % 1439];
				while (!open_list.empty()) {
					auto min_f = std::min_element(std::begin(fScore), std::end(fScore));
					size_t min_f_p = std::distance(std::begin(fScore), min_f);
					//find(open_list.begin(), open_list.end(), min_f_p);
					size_t current_node = min_f_p;
					if (current_node == end) {//回溯
						Shortest_path.push_back(current_node);
						do {
							current_node = came_from[current_node];
							Shortest_path.push_back(current_node);
						} while (current_node != start);
						std::reverse(Shortest_path.begin(), Shortest_path.end());
						m_shortest_path_all_day[static_cast<int>(present_time) % 1439][start].insert(std::make_pair(end, Shortest_path));
						break;
					}
					if (find(open_list.begin(), open_list.end(), current_node) == open_list.end()) {
						system("pause");
					}
					open_list.erase(find(open_list.begin(), open_list.end(), current_node));
					closeed_list.push_back(current_node);
					for (auto &i : m_road_net.get_road_net()[current_node]->next) {
						if (find(closeed_list.begin(), closeed_list.end(), i->id) != closeed_list.end()) {
							continue;
						}
						//时间要动态变化
						double temp_gScore = gScore[current_node] + m_road_net.get_node_dis()[current_node][i->id] / m_road_net.get_road_net()[i->id]->datum_t.v[static_cast<int>(present_time) % 1439];
						if (find(closeed_list.begin(), closeed_list.end(), i->id) == closeed_list.end()) {
							open_list.push_back(i->id);
						}
						else if (temp_gScore >= gScore[i->id])
							continue;
						came_from[i->id] = current_node;
						gScore[i->id] = temp_gScore;
						fScore[i->id] = gScore[i->id] + m_road_net.get_node_dis()[i->id][end] / m_road_net.get_road_net()[i->id]->datum_t.v[static_cast<int>(present_time) % 1439];
					}
					gScore[current_node] = std::numeric_limits<double>::max();
					fScore[current_node] = std::numeric_limits<double>::max();
				}
				if (open_list.empty())
					std::cout << std::endl << "无法找到点: " << start << "-->" << end << "之间的路径！" << std::endl;
			}
		}

		void dynamic_vrp::shortest_path(size_t start, size_t end, std::vector<size_t> &Shortest_path, Real present_time) const{
			//A star
			if (start == end) {
				std::cout << "The start and end of shotest_path is equal" << std::endl;
				system("pause");
			}
			else {
				std::vector<size_t> open_list;
				std::vector<size_t> closeed_list;
				open_list.push_back(start);
				std::vector<int> came_from(m_nodes_size);
				/*for (size_t i = 0; i < m_road_net.get_road_net()[start]->next.size(); ++i) {
					open_list.push_back(m_road_net.get_road_net()[start]->next[i]->id);
				}*/
				std::vector<double> fScore(m_nodes_size);
				std::vector<double> gScore(m_nodes_size);
				for (auto &i : gScore)
					i = std::numeric_limits<double>::max();//km(Infinity)
				gScore[start] = 0;
				for (auto &i : fScore)
					i = std::numeric_limits<double>::max();// Infinity
				fScore[start] = m_road_net.get_node_dis()[start][end] / m_road_net.get_road_net()[start]->datum_t.v[static_cast<int>(present_time) % 1439];
				while (!open_list.empty()) {
					auto min_f = std::min_element(std::begin(fScore), std::end(fScore));
					size_t min_f_p = std::distance(std::begin(fScore), min_f);
					//find(open_list.begin(), open_list.end(), min_f_p);
					size_t current_node = min_f_p;
					if (current_node == end) {//回溯
						Shortest_path.push_back(current_node);
						do {
							current_node = came_from[current_node];
							Shortest_path.push_back(current_node);
						} while (current_node != start);
						std::reverse(Shortest_path.begin(), Shortest_path.end());
						break;
					}
					if (find(open_list.begin(), open_list.end(), current_node) == open_list.end()) {
						system("pause");
					}
					open_list.erase(find(open_list.begin(), open_list.end(), current_node));
					closeed_list.push_back(current_node);
					for (auto &i : m_road_net.get_road_net()[current_node]->next) {
						if (find(closeed_list.begin(), closeed_list.end(), i->id) != closeed_list.end()) {
							continue;
						}
						//时间要动态变化
						double temp_gScore = gScore[current_node] + m_road_net.get_node_dis()[current_node][i->id] / m_road_net.get_road_net()[i->id]->datum_t.v[static_cast<int>(present_time) % 1439];
						if (find(closeed_list.begin(), closeed_list.end(), i->id) == closeed_list.end()) {
							open_list.push_back(i->id);
						}
						else if (temp_gScore >= gScore[i->id])
							continue;
						came_from[i->id] = current_node;
						gScore[i->id] = temp_gScore;
						fScore[i->id] = gScore[i->id] + m_road_net.get_node_dis()[i->id][end] / m_road_net.get_road_net()[i->id]->datum_t.v[static_cast<int>(present_time) % 1439];
					}
					gScore[current_node] = std::numeric_limits<double>::max();
					fScore[current_node] = std::numeric_limits<double>::max();
				}
				if (open_list.empty())
					std::cout << std::endl << "无法找到点: " << start << "-->" << end << "之间的路径！" << std::endl;
			}
		}

		void dynamic_vrp::calcul_shorest_dis_24h(bool calculated)
		{
			/*DWORD start_time = GetTickCount();
			auto t = get_cost_time(m_depart_time, 0, 880);
			auto t1 = get_cost_time(m_depart_time, 0, 96);
			auto t2 = get_cost_time(m_depart_time, 1, 880);
			auto d = get_shortest_dis(m_depart_time, 0, 880);
			auto d1 = get_shortest_dis(m_depart_time, 0, 880);
			DWORD end_time = GetTickCount();
			std::cout << std::endl << "The run time is:" << (end_time - start_time) << "ms!" << std::endl;*/
		}

		void dynamic_vrp::calcul_shorest_dis_1min(int i)
		{
			//std::string shorest_path_fileName;
			std::string cost_time_fileName;
			std::string shortest_dis_fileName;

			//shorest_path_fileName = "./instance/problem/realworld/DVRP/data/shorest_path/shorest_path" + std::to_string(i) + ".txt";
			cost_time_fileName = "./instance/problem/realworld/DVRP/data/OffLine/cost_time/cost_time" + std::to_string(i) + ".txt";
			shortest_dis_fileName = "./instance/problem/realworld/DVRP/data/OffLine/shortest_dis/shortest_dis" + std::to_string(i) + ".txt";
			//std::ofstream shortest_path_t(shorest_path_fileName);
			std::ofstream cost_time(cost_time_fileName);
			std::ofstream shortest_dis(shortest_dis_fileName);

			//std::stringstream shortest_path_str;
			std::stringstream cost_time_str;
			std::stringstream shortest_dis_str;

			//std::vector<std::vector<size_t>> all_path;
			std::vector<std::vector<Real>> cost_time_t(m_nodes_size);
			std::vector<std::vector<Real>> shortest_dis_t(m_nodes_size);

			for (size_t ii = 0; ii < m_road_net.get_road_net().size(); ++ii) {
				for (size_t jj = 0; jj < m_road_net.get_road_net().size(); ++jj) {
					if (ii == jj) {
						cost_time_t[ii].push_back(std::numeric_limits<double>::max());
						//cost_time << std::numeric_limits<double>::max() << " ";
						cost_time_str << std::numeric_limits<double>::max() << " ";

						shortest_dis_t[ii].push_back(std::numeric_limits<double>::max());
						//shortest_dis << std::numeric_limits<double>::max() << " ";
						shortest_dis_str << std::numeric_limits<double>::max() << " ";

						//all_path.push_back({ std::numeric_limits<int>::max() });
						//shortest_path_t << std::numeric_limits<int>::max() << "\n";
						//shortest_path_str << std::numeric_limits<int>::max() << "\n";
					}
					else
					{
						Real sum = 0.0;
						Real time = i;

						std::vector<size_t> Shortest_path;
						//std::cout << ii << " " << jj << std::endl;
						auto it_path = m_shortest_path_all_day[i % 1439][ii].find(jj);
						if (it_path != m_shortest_path_all_day[i % 1439][ii].end()) {
							Shortest_path = it_path->second;
						}
						else
						{
							shortest_path(ii, jj, Shortest_path, i);
						}
						//all_path.push_back(Shortest_path);

						for (int k = 0; k < Shortest_path.size() - 1; ++k) {
							sum += m_road_net.get_node_dis()[Shortest_path[k]][Shortest_path[k + 1]];
							time += (m_road_net.get_node_dis()[Shortest_path[k]][Shortest_path[k + 1]] / m_road_net.get_road_net()[Shortest_path[k + 1]]->datum_t.v[static_cast<int>(time)%1439]) * 60;
							//shortest_path_t << Shortest_path[k] << " ";
							//shortest_path_str << Shortest_path[k] << " ";
						}
						//shortest_path_t << Shortest_path[Shortest_path.size() - 1];
						//shortest_path_t << "\n";
						//shortest_path_str << Shortest_path[Shortest_path.size() - 1];
						//shortest_path_str << "\n";

						cost_time_t[ii].push_back(time - i);//减去出发时刻i
						shortest_dis_t[ii].push_back(sum);

						cost_time_str << time - i << " ";

						shortest_dis_str << sum << " ";
					}

				}
				cost_time_str << "\n";
				shortest_dis_str << "\n";
			}
			//auto mm= shortest_path_str.str();
			//shortest_path_t << shortest_path_str.str();
			cost_time << cost_time_str.str();
			shortest_dis << shortest_dis_str.str();

			//shortest_path_t.close();
			cost_time.close();
			shortest_dis.close();

			/*auto it_path = shorest_path_all_day.begin();
			for (int j = 0; j < i; ++j) {
				it_path++;
			}
			*it_path = all_path;*/

			/*auto it_time = cost_time_all_day.begin();
			for (int j = 0; j < i; ++j) {
				it_time++;
			}
			*it_time = cost_time_t;

			auto it_dis = shortest_dis_all_day.begin();
			for (int j = 0; j < i; ++j) {
				it_dis++;
			}
			*it_dis = shortest_dis_t;*/

			/*cost_time_all_day.push_back(cost_time_t);
			shortest_dis_all_day.push_back(shortest_dis_t);*/
		}

		Real dynamic_vrp::get_cost_time(Real present_time, size_t start, size_t end) 
		{
			int present_time_ = static_cast<int>(present_time) % 1439;
			Real cost_time = 0.0;
			bool is_calcu = true;
			if (m_cost_time_all_day[present_time_][start].find(end)!= m_cost_time_all_day[present_time_][start].end()) {
				auto it = m_cost_time_all_day[present_time_][start].find(end);
				is_calcu = false;
				cost_time = it->second;
				return it->second;
			}
			if(is_calcu==true) {
				Real time = present_time;
				std::vector<size_t> Shortest_path;
				auto it_path = m_shortest_path_all_day[present_time_][start].find(end);
				if (it_path != m_shortest_path_all_day[present_time_][start].end()) {
					Shortest_path = it_path->second;
				}
				else
				{
					shortest_path(start, end, Shortest_path, present_time);
				}
				for (int k = 0; k < Shortest_path.size() - 1; ++k) {
					time += (m_road_net.get_node_dis()[Shortest_path[k]][Shortest_path[k + 1]] / m_road_net.get_road_net()[Shortest_path[k + 1]]->datum_t.v[static_cast<int>(time) % 1439]) * 60;
				}
				m_cost_time_all_day[present_time_][start].insert(std::make_pair(end, time - present_time));
				cost_time = time - present_time;
				return time - present_time;
			}
		}

		Real dynamic_vrp::get_shortest_dis(Real present_time, size_t start, size_t end)
		{
			Real dis = 0.0;
			int present_time_ = static_cast<int>(present_time) % 1439;
			bool is_calcu = true;
			if (m_shortest_dis_all_day[present_time_][start].find(end) != m_shortest_dis_all_day[present_time_][start].end()) {
				auto it = m_shortest_dis_all_day[present_time_][start].find(end);
				is_calcu = false;
				dis = it->second;
				return it->second;
			}
			if (is_calcu == true) {
				std::vector<size_t> Shortest_path;
				auto it_path = m_shortest_path_all_day[present_time_][start].find(end);
				if (it_path != m_shortest_path_all_day[present_time_][start].end()) {
					Shortest_path = it_path->second;
				}
				else
				{
					shortest_path(start, end, Shortest_path, present_time);
				}
				for (int k = 0; k < Shortest_path.size() - 1; ++k) {
					dis += m_road_net.get_node_dis()[Shortest_path[k]][Shortest_path[k + 1]];
				}
				m_shortest_dis_all_day[present_time_][start].insert(std::make_pair(end, dis));
				return dis;
			}
		}

		

		void dynamic_vrp::visit_routes(routes &x) const {
			for (size_t i = 0; i < x.m_path.size(); ++i) {
				double present_time = m_depart_time;
				for (size_t j = 0; j < x.m_path[i].size() - 1; ++j) {
					Real node_distance = m_road_net.get_node_dis()[x.m_path[i][j]][x.m_path[i][j + 1]];
					Real arriveTime = 0;
					Real travelTime = 0;
					size_t serverTime = 0;
					//if (present_time >= m_final_back_time) present_time = m_final_back_time;//会出现超出仓库关闭时间的情况
					Real velocity = m_road_net.get_road_net()[x.m_path[i][j + 1]]->datum_t.v[static_cast<int>(present_time)%1439];
					travelTime = node_distance / velocity * 60;
					/*for (auto &row : m_objective) {
						row(x, present_time, i, j);
					}*/
					if (x.m_path_node_type[i][j] == node_type::current_customer || x.m_path_node_type[i][j] == node_type::future_customer) {
						serverTime = m_road_net.get_road_net()[x.m_path[i][j]]->server_time;
					}
					arriveTime = present_time + serverTime + travelTime;
					present_time = arriveTime;
					for (auto &row : m_objective) {
						row(x, present_time, i, j);
					}
					//auto it2 = find(x.m_cus_order[i].begin(), x.m_cus_order[i].end(), x.m_path[i][j + 1]);
				}
			}
		}

		void dynamic_vrp::over_load(routes &x) {
			for (size_t i = 0; i < x.m_cus_order.size(); ++i) {
				x.m_hi.over_load[i] = 0;
				for (size_t j = 1; j < x.m_cus_order[i].size(); ++j) {
					if(m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.is_served==false)
						x.m_hi.total_delivery_weight[i] += m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.delivery_weight;
				}
				if (x.m_hi.total_delivery_weight[i] > m_vehicle_property[x.m_vehicle_type[i]].capacity) {
					x.m_hi.over_load[i] += x.m_hi.total_delivery_weight[i] - m_vehicle_property[x.m_vehicle_type[i]].capacity;
				}
				for (size_t j = 1; j < x.m_cus_order[i].size(); ++j) {
					if (m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.is_served == false)
						x.m_hi.total_delivery_weight[i] = x.m_hi.total_delivery_weight[i] - m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.delivery_weight + m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.pick_weight;
					else
						x.m_hi.total_delivery_weight[i] = x.m_hi.total_delivery_weight[i] - 0 + m_road_net.get_road_net()[x.m_cus_order[i][j]]->datum_t.pick_weight;

					if (x.m_hi.total_delivery_weight[i] > m_vehicle_property[x.m_vehicle_type[i]].capacity) {
						x.m_hi.over_load[i] += x.m_hi.total_delivery_weight[i] - m_vehicle_property[x.m_vehicle_type[i]].capacity;
					}
				}
			}
		}

		void dynamic_vrp::over_depot_end_time(routes &x) {
			for (size_t i = 0; i < x.m_hi.arrive_time.size(); i++) {
				if (x.m_hi.arrive_time[i][x.m_hi.arrive_time[i].size() - 1] > m_final_back_time)
					x.m_hi.over_time[i] = x.m_hi.arrive_time[i][x.m_hi.arrive_time[i].size() - 1] - m_final_back_time;
				else x.m_hi.over_time[i] = 0;
			}
		}

	}
}
