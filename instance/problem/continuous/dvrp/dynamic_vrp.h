/********* Begin Register Information **********
{
	"name": "DVRP",
	"identifier": "DVRP",
	"problem tags": [ "DVRP" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Long Xiao
* Email: 917003976@qq.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://gitlab.com/OFEC/console for more information
*
*********************************************************************************/
// updated Feb 8, 2019 by Long Xiao

#ifndef OEFC_DYNAMIC_VRP_H
#define OEFC_DYNAMIC_VRP_H
#include "../../../../core/problem/problem.h"
#include "road_net.h"
#include <map>

namespace ofec {

#define DVRP_CAST(pro) dynamic_cast<DVRP::dynamic_vrp&>(pro)

	namespace DVRP {
		struct vehicle {//vehicle property
			size_t capacity;
		};

		struct HI {//HI(heuristic information)
			int used_car;
			std::vector<std::vector<Real>> delay_time;
			std::vector<std::vector<Real>> wait_time;
			std::vector<std::vector<Real>> arrive_time;
			std::vector<Real> route_length;
			std::vector<Real> over_time;
			std::vector<Real> over_load;
			std::vector<Real> total_delivery_weight;
		};

		class routes : public VariableBase{//encode of solution
		public:
			std::vector<size_t> m_vehicle_type;       //e.g.:{1,2} 
			std::vector<std::vector<size_t>> m_cus_order;//e.g.:{{0,1,2,3,4,0},{0,11,45,8,0}}
			std::vector<std::vector<size_t>> m_path;     //e.g.:{{0,13,15,17,1,7,6,2,55,67,78,3,35,65,23,4,7,10,0},{0,76,34,2,11,13,54,45,90,110,8,6,5,0}}
			std::vector<std::vector<node_type>> m_path_node_type;

			HI m_hi;
		};

		class dynamic_vrp : public Problem {
		protected:
			void initialize_() override;
			void evaluate_(SolutionBase &s, bool effective) override;

		public:
			void initializeSolution(SolutionBase &s, Random *rnd) const override;
			bool same(const SolutionBase &s1, const SolutionBase &s2) const override;
			Real variableDistance(const SolutionBase &s1, const SolutionBase &s2) const override;
			Real variableDistance(const VariableBase &s1, const VariableBase &s2) const override;

			void set_Online(bool is_online);
			void construct_path(SolutionBase &s);
			std::vector<bool> is_valid(SolutionBase &s) const;
			
			void shortest_path(size_t, size_t, std::vector<size_t>&, Real)const;
			void shortest_path(size_t, size_t, std::vector<size_t>&, Real);
			
			Real get_cost_time(Real present_time, size_t start, size_t end);
			Real get_shortest_dis(Real present_time, size_t start, size_t end);
			const road_net &get_net() const { return m_road_net; }
			const size_t get_depot() const { return m_depot; }
			const Real get_depart_time() const { return m_depart_time; }
			const Real get_final_back_time() const { return m_final_back_time; }
			const std::vector<size_t> &get_current_cus_order() const { return m_current_cus_order; }
			const std::vector<std::pair<int, int>> &get_future_cus_order() const { return m_future_cus_order; }
			std::vector<std::pair<int, int>> &get_future_cus_order(){ return m_future_cus_order; }
			const std::vector<vehicle> &get_vehicle_property() const { return m_vehicle_property; }
			const int &get_num_futureOrder() const { return m_future_order_num; }

		private:
			void sort_customers(std::vector<size_t> &) const;
			void construct_path(SolutionBase &s)const;
			void calcul_shorest_dis_24h(bool calculated);
			void calcul_shorest_dis_1min(int i);
			//caculate objectives and constraints
			void visit_routes(routes &) const;

			void used_car(routes &x) {
				x.m_hi.used_car = x.m_cus_order.size();
			}

			void length(routes &x, Real presentTime, size_t i, size_t j) {
				Real node_distance = m_road_net.get_node_dis()[x.m_path[i][j]][x.m_path[i][j + 1]];
				x.m_hi.route_length[i] += node_distance;
			}

			void delay_time(routes &x, Real presentTime, size_t i, size_t j) {
				Real delayTime = 0;
				Real endTime = 0;
				//auto it2 = find(x.m_cus_order[i].begin(), x.m_cus_order[i].end(), x.m_path[i][j + 1]);
				if (x.m_path_node_type[i][j + 1] == node_type::current_customer || x.m_path_node_type[i][j + 1] == node_type::future_customer) {
					endTime = m_road_net.get_road_net()[x.m_path[i][j + 1]]->datum_t.end_time;
					if (presentTime > endTime)
						delayTime = presentTime - endTime;
					x.m_hi.delay_time[i].push_back(delayTime);
				}
			}

			void arrive_time(routes &x, Real presentTime, size_t i, size_t j) {
				if (x.m_path_node_type[i][j + 1] != node_type::general_node) {
					x.m_hi.arrive_time[i].push_back(presentTime);
				}
			}

			void wait_time(routes &x, Real presentTime, size_t i, size_t j) {
				Real waitTime = 0;
				Real beginTime = 0;
				//auto it2 = find(x.m_cus_order[i].begin(), x.m_cus_order[i].end(), x.m_path[i][j + 1]);
				if (x.m_path_node_type[i][j + 1] == node_type::current_customer || x.m_path_node_type[i][j + 1] == node_type::future_customer) {
					beginTime = m_road_net.get_road_net()[x.m_path[i][j + 1]]->datum_t.begin_time;
					if (presentTime < beginTime)
						waitTime = beginTime - presentTime;
					x.m_hi.wait_time[i].push_back(waitTime);//�ȴ�ʱ����Ϊ���Ա�Ƿ����?
				}
			}
			void over_load(routes &);
			void over_depot_end_time(routes &);

		protected:
			road_net m_road_net;
			size_t m_depot = 175;
			size_t m_nodes_size = 0;
			size_t m_current_order_num = 10;
			size_t m_future_order_num = 5;
			std::string m_road_file_name;
			std::string m_velocity_file_name_Offline;
			std::string m_velocity_file_name_Online;
			std::vector<vehicle> m_vehicle_property;    //index represents vehicle type
			bool m_is_vehilce_compatible = false;		//if false, different vehicle can not interaction

			std::vector<size_t> m_current_cus_order;
			std::vector<std::pair<int, int>> m_future_cus_order;//<id,time>
			//std::vector<size_t> m_future_cus_order;

			//Real m_average_cost;
			Real m_depart_time = 8 * 60 + 0;     //480 min
			Real m_final_back_time = 24 * 60 + 0;//1440 min
			Real m_dynamic_order_endtime = 14 * 60 + 0;//840 min
			std::vector<std::function <void(routes &, Real, size_t, size_t)>> m_objective;

		private:
			std::vector<std::vector<std::multimap<int, std::vector<size_t>>>>m_shortest_path_all_day;
			std::vector<std::vector<std::multimap<int, Real>>>m_cost_time_all_day;
			std::vector<std::vector<std::multimap<int, Real>>>m_shortest_dis_all_day;
			bool m_is_Online = false;

		};
	}
	using dynamic_vrp = DVRP::dynamic_vrp;
}
#endif // !OEFC_DYNAMIC_VRP
