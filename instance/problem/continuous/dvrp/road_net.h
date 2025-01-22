#ifndef OFEC_ROAD_NET_H
#define OFEC_ROAD_NET_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <deque>
#include <cstdlib>

namespace ofec {
	namespace DVRP {
		struct datum {
			size_t road_num = 0;
			int id = -1;//mapping flag for building a road network
			double x, y;
			double pick_weight, delivery_weight;
			int begin_time, end_time;
			bool is_served = false;
			std::vector<double> v;//store 24h speed value, interval 1min unit:km/h
			std::vector<double> v_Offline;//store 24h speed value, interval 1min unit:km/h
			std::vector<double> v_Online;//store 24h speed value, interval 1min unit:km/h
		};

		enum class node_type { general_node, depot, current_customer, future_customer };

		struct road_node {
			datum datum_t;
			std::vector<road_node *> prev;
			std::vector<road_node *> next;
			int id = -1;             //node number, for select customers
			size_t server_time = 30; //unit:min
			node_type type = node_type::general_node;

		};

		class road_net {
		private:
			std::vector<std::vector<datum>> m_roads;
			std::vector<road_node *> m_road_nodes;  //road net
			std::vector<std::pair<int,int>> m_customers;
			std::vector<std::vector<double>> m_node_dis;
			//double m_average_dis;
		public:
			void read_data(const std::string &, const std::string &, const std::string &velocity_file_name_Online);
			void set_velocity(bool is_Online);
			void constrct_net(size_t);
			
			double node_distance(road_node *, road_node *) const;
			void node_distance();
			void select_customers(const size_t, const size_t, const size_t, const double, const double, Random *rnd);
			const std::vector<std::vector<datum>>& get_roads() const { return m_roads; }
			const std::vector<road_node *>& get_road_net() const { return m_road_nodes; }
			const std::vector<std::pair<int, int>>& get_customers() const { return m_customers; } //{ {node_index,disclosure_time},{}}
			const std::vector<std::vector<double>>& get_node_dis() const { return m_node_dis; }
			~road_net();
		};	
	}
}
#endif // !OFEC_ROAD_NET
