#include "road_net.h"
#include "../../../../core/global.h"
#include <numeric>

namespace ofec {
	namespace DVRP {
		//split string
		std::vector<std::string> split(std::string str, std::string pattern)
		{
			std::string::size_type pos;
			std::vector<std::string> routes;
			str += pattern;
			size_t size = str.size();

			for (size_t i = 0; i < size; i++)
			{
				pos = str.find(pattern, i);
				if (pos < size)
				{
					std::string s = str.substr(i, pos - i);
					routes.push_back(s);
					i = pos + pattern.size() - 1;
				}
			}
			return routes;
		}

		void road_net::read_data(const std::string &road_file_name, const std::string &velocity_file_name_Offline, const std::string &velocity_file_name_Online) {

			std::string road_file = static_cast<std::string>(g_working_directory) + "/instance/problem/realworld/dvrp/" + road_file_name + ".txt";
			std::string velocity_file_Offline = static_cast<std::string>(g_working_directory) + "/instance/problem/realworld/dvrp/OffLine/" + velocity_file_name_Offline + ".txt";
			std::string velocity_file_Online = static_cast<std::string>(g_working_directory) + "/instance/problem/realworld/dvrp/OnLine/" + velocity_file_name_Online + ".txt";

			m_roads.clear();
			std::fstream read_velocity_Offline(velocity_file_Offline);
			std::fstream read_velocity_Online(velocity_file_Online);
			std::fstream read_coordinate(road_file);
			if (read_coordinate && read_velocity_Offline && read_velocity_Online) {
				size_t roadNum(0);
				std::string line_coordinate, temp, line_velocity_Offline, line_velocity_Online;
				while (getline(read_coordinate, line_coordinate) && getline(read_velocity_Offline, line_velocity_Offline)&& getline(read_velocity_Online, line_velocity_Online)) {
					std::stringstream iss(line_coordinate);
					iss >> temp;
					if (roadNum != stoi(temp)) {
						m_roads.push_back({});
						roadNum = stoi(temp);
					}
					iss >> temp;
					iss >> temp;
					iss >> temp;
					double x, y;
					iss >> temp;
					x = stod(temp);
					iss >> temp;
					y = stod(temp);

					double pickWeight, deliveryWeight;
					iss >> temp;
					deliveryWeight = stod(temp);
					iss >> temp;
					pickWeight = stod(temp);

					int beginTime, endTime;
					std::vector<std::string> begin_time_temp;
					std::vector<std::string> end_time_temp;
					iss >> temp;
					begin_time_temp = split(temp, ":");
					beginTime = stoi(begin_time_temp[0]) * 60 + stoi(begin_time_temp[1]);
					iss >> temp;
					end_time_temp = split(temp, ":");
					endTime = stoi(end_time_temp[0]) * 60 + stoi(end_time_temp[1]);

					m_roads.back().push_back({ roadNum,-1,x,y,pickWeight,deliveryWeight,beginTime,endTime });
					std::stringstream row_Offline(line_velocity_Offline);
					std::stringstream row_Online(line_velocity_Online);

					m_roads.back().back().v_Offline.clear();
					m_roads.back().back().v_Offline.resize(1440);
					for (size_t i = 0; i < 1440; i++) {
						row_Offline >> m_roads.back().back().v_Offline[i];
					}

					m_roads.back().back().v_Online.clear();
					m_roads.back().back().v_Online.resize(1440);
					for (size_t i = 0; i < 1440; i++) {
						row_Online >> m_roads.back().back().v_Online[i];
					}
					
				}
			}
			else {
				std::cout << "Open file failed" << std::endl;
			}
			read_coordinate.close();
			read_velocity_Offline.close();
			read_velocity_Online.close();
		}
		void road_net::set_velocity(bool is_Online)
		{
			if (is_Online) {
				for (auto &i : m_roads) {
					for (auto &j : i) {
						j.v = j.v_Online;
					}
				}
			}
			else {
				for (auto &i : m_roads) {
					for (auto &j : i) {
						j.v = j.v_Offline;
					}
				}
			}
		}
		void road_net::constrct_net(size_t depot) {
			m_road_nodes.clear();
			size_t 	cnt = 0;
			for (size_t i = 0; i < m_roads.size(); i++)
			{
				for (auto& t_node : m_roads[i]) {
					bool flag = false;
					for (size_t j = 0; j < cnt; j++)
					{
						if (t_node.x == m_road_nodes[j]->datum_t.x&&t_node.y == m_road_nodes[j]->datum_t.y)
						{
							flag = true;
							if (t_node.pick_weight != m_road_nodes[j]->datum_t.pick_weight || t_node.delivery_weight != m_road_nodes[j]->datum_t.delivery_weight || t_node.begin_time != m_road_nodes[j]->datum_t.begin_time || t_node.end_time != m_road_nodes[j]->datum_t.end_time) {
								std::cout << t_node.road_num << ' ' << t_node.x << ' ' << t_node.y << std::endl;
							}
							t_node.id = j;
							break;
						}
					}
					if (!flag) {
						m_road_nodes.push_back(new road_node({ t_node }));
						m_road_nodes.back()->id = m_road_nodes.size() - 1;
						t_node.id = m_road_nodes.size() - 1;
					}
				}
				cnt = m_road_nodes.size();
			}

			for (size_t i = 0; i < m_roads.size(); i++) {
				for (size_t k = 0; k < m_roads[i].size() - 1; ++k) {
					m_road_nodes[m_roads[i][k].id]->next.push_back(m_road_nodes[m_roads[i][k + 1].id]);
					m_road_nodes[m_roads[i][k + 1].id]->prev.push_back(m_road_nodes[m_roads[i][k].id]);
				}
			}
			m_road_nodes[depot]->datum_t.begin_time = 0;
			m_road_nodes[depot]->datum_t.end_time = 1440;
			m_road_nodes[depot]->server_time = 0;
			m_road_nodes[depot]->datum_t.pick_weight = 0;
			m_road_nodes[depot]->datum_t.delivery_weight = 0;
			m_road_nodes[depot]->type = node_type::depot;

			/*All customers print */
			std::ofstream toal_cus_nodes("./instance/problem/realworld/DVRP/data/all_customers.txt");
			for (size_t i = 0; i < m_road_nodes.size(); i++)
			{
				toal_cus_nodes << std::fixed << std::setprecision(7) << m_road_nodes[i]->id << " " << m_road_nodes[i]->datum_t.x << " "<< m_road_nodes[i]->datum_t.y;
				toal_cus_nodes << std::endl;
			}
			toal_cus_nodes.close();

			//std::cout << "Verify if there is a breakpoint " << std::endl;
			for (size_t i = 0; i < m_road_nodes.size(); i++)
			{
				if (m_road_nodes[i]->next.size() < 1) {
					std::cout << "Breakpoint is£º" << std::fixed << std::setprecision(7) << m_road_nodes[i]->datum_t.road_num << " " << m_road_nodes[i]->id << " " << m_road_nodes[i]->datum_t.x << " " << m_road_nodes[i]->datum_t.y << std::endl;
				}
			}
		}

		double road_net::node_distance(road_node *n1, road_node *n2) const{
			const double pi = 3.14159265358979323846264338328;
			const double earth_radius = 6378.137;//km
			double dis_x = (n1->datum_t.x - n2->datum_t.x)*pi / 180.0;
			double dis_y = (n1->datum_t.y - n2->datum_t.y)*pi / 180.0;
			double dis;
			dis = earth_radius * 2 * asin(sqrt(pow(sin(dis_y / 2), 2) + cos(n1->datum_t.y*pi / 180.0)*cos(n2->datum_t.y*pi / 180.0)*pow(sin(dis_x / 2), 2)));
			return dis;
		}
		void road_net::node_distance() {
			static double pi = 3.14159265358979323846264338328;
			static double earth_radius = 6378.137;//km
			m_node_dis.clear();
			m_node_dis.resize(m_road_nodes.size());
			for (size_t i = 0; i < m_road_nodes.size(); i++) {
				for (size_t j = 0; j < m_road_nodes.size(); j++)
				{
					double dis_x = (m_road_nodes[i]->datum_t.x - m_road_nodes[j]->datum_t.x)*pi / 180.0;
					double dis_y = (m_road_nodes[i]->datum_t.y - m_road_nodes[j]->datum_t.y)*pi / 180.0;
					double dis = earth_radius * 2 * asin(sqrt(pow(sin(dis_y / 2), 2) + cos(m_road_nodes[i]->datum_t.y*pi / 180.0)*cos(m_road_nodes[j]->datum_t.y*pi / 180.0)*pow(sin(dis_x / 2), 2)));
					m_node_dis[i].push_back(dis);
				}
			}
		}

		void road_net::select_customers(const size_t depot, const size_t static_order_num, const size_t dynamic_order_num, const double depart_time, const double dynamic_endtime, Random *rnd) {

			m_customers.clear();
			/**********traversing graph************************************************/
			std::deque<int> q;
			int present_node = depot, last_node = depot;
			std::vector<size_t> visited;
			std::vector<int> customer_temp;
			int static_cnt = 0;
			int dynamic_cnt = 0;
			visited.resize(m_road_nodes.size());
			for (size_t i = 0; i < m_road_nodes.size(); i++) {
				visited[i] = false;
			}
			q.push_front(depot);
			//set<int> all_custom;
			while (!q.empty()) {
				int node = q.front();
				q.pop_front();
				if (!visited[node]) {
					visited[node] = true;
					/*present_node = node;
					if (m_Distance(m_road_net[present_node], m_road_net[last_node]) >= customer_dis) {
						++cnt;
						std::cout << node << "-";
						m_Customers.push_back(present_node);
						last_node = present_node;
					}*/
					customer_temp.push_back(node);
				}
				//if (cnt == static_order_num) break;
				for (size_t i = 0; i < m_road_nodes[node]->next.size(); ++i) {
					if (!visited[m_road_nodes[node]->next[i]->id]) {
						q.push_back(m_road_nodes[node]->next[i]->id);
					}
				}
				if (q.empty()) {
					for (size_t i = 0; i < m_road_nodes.size(); i++) {
						if (!visited[i]) {
							q.push_back(i);
							std::cout << "There is node not accessed" << std::endl;
							break;
						}
					}
				}
				//last_node = node;
			}
			/***Choose a certain number of customers£¨argument£ºnum£©*****************************************/

			/*Select customers by traversing the graph*/

			/*for (size_t i = 0; i < customer_temp.size(); ++i) {
				present_node = customer_temp[i];
				if (node_distance(m_road_net[present_node], m_road_net[last_node]) >= m_customer_dis) {
					++cnt;
					std::cout << customer_temp[i] << "-";
					m_customers.push_back(present_node);
					m_road_net[present_node]->type = 2;
					last_node = present_node;
				}
				if (cnt == static_order_num) break;
			}*/

			/*Select customer by text file node number*/
			/*for (size_t i = 0; i < m_road_net.size(); ++i) {
				present_node = i;
				if (node_distance(m_road_net[present_node], m_road_net[last_node]) >= m_customer_dis) {
					++cnt;
					std::cout << m_road_net[i]->id << "-";
					m_customers.push_back(present_node);
					m_road_net[present_node]->type = 2;
					last_node = present_node;
				}
				if (cnt == static_order_num) break;
			}*/
		
			/*Randomly choose customers*/
			std::vector<int> random_sequence(m_road_nodes.size());
			std::iota(random_sequence.begin(), random_sequence.end(), 0);
			rnd->uniform.shuffle(random_sequence.begin(), random_sequence.end());

			std::vector<int> random_disclosure_time1(depart_time - 0 + 1);
			std::iota(random_disclosure_time1.begin(), random_disclosure_time1.end(), 0);
			rnd->uniform.shuffle(random_disclosure_time1.begin(), random_disclosure_time1.end());

			//static order
			std::cout << "static order: ";
			for (size_t i = 0; i < random_sequence.size(); ++i) {
				if (random_sequence[i] == depot) continue;
				std::cout << random_sequence[i] << ", ";
				std::pair<int, int> static_temp;
				static_temp.first = random_sequence[i];
				static_temp.second = -1;
				m_customers.push_back(static_temp);
				m_road_nodes[random_sequence[i]]->type = node_type::current_customer;
				//m_road_net[random_sequence[i]]->is_dynamic = false;
				static_cnt++;
				if (static_cnt == static_order_num) break;
			}
			std::cout << std::endl;

			//dynamic order
			/*std::vector<int> random_disclosure_time2(dynamic_endtime - depart_time + 1);
			std::iota(random_disclosure_time2.begin(), random_disclosure_time2.end(), depart_time);
			global::ms_global->m_uniform[caller::Problem]->shuffle(random_disclosure_time2.begin(), random_disclosure_time2.end());*/
			std::cout << "dynamic order: ";
			for (size_t i = 0; i < random_sequence.size(); ++i) {
				if (random_sequence[i] == depot) continue;
				if (m_road_nodes[random_sequence[i]]->type == node_type::current_customer) continue;
				std::cout << random_sequence[i] << ", ";
				std::pair<int, int> dynamic_temp;
				dynamic_temp.first = random_sequence[i];

				int beginTime = m_road_nodes[dynamic_temp.first]->datum_t.begin_time;
				std::vector<int> random_disclosure_time2(beginTime - 60 - depart_time + 1);
				std::iota(random_disclosure_time2.begin(), random_disclosure_time2.end(), depart_time);
				rnd->uniform.shuffle(random_disclosure_time2.begin(), random_disclosure_time2.end());
				dynamic_temp.second = random_disclosure_time2[i % random_disclosure_time2.size()];

				m_customers.push_back(dynamic_temp);
				m_road_nodes[random_sequence[i]]->type = node_type::future_customer;
				m_road_nodes[random_sequence[i]]->datum_t.delivery_weight = 0;//new customers can only pick
				dynamic_cnt++;
				if (dynamic_cnt == dynamic_order_num) break;
			}
			std::cout << "\n";
			//std::ofstream customers("./instance/problem/realworld/DVRP/data/customers.txt");
			//for (size_t i = 0; i < m_current_customers.size(); i++)
			//{
			//	customers << m_current_customers[i] << " ";
			//}
			//customers.close();
			//return  m_current_customers;
		}

		road_net::~road_net() {
			if (!m_road_nodes.empty()) {
				for (size_t i = 0; i < m_road_nodes.size(); i++) {
					delete m_road_nodes[i];
				}
			}
		}

	}
}