#include "water_resource_plan.h"
#include <fstream>
#include "../../../../../core/global.h"
#include"../../../../../utility/nondominated_sorting/fast_sort.h"
#include"../../../../../utility/nondominated_sorting/filter_sort.h"
#include <algorithm>
#include<cmath>

namespace ofec {
	/*
	2 objectives, each objective have only one peak
	*/
	void Water_resource_plan::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeObjective(6);
		for (size_t i = 0; i < m_number_objectives; ++i) {
			m_optimize_mode[i] = OptimizeMode::kMinimize;
		}
		resizeVariable(3);//n=3
		for (size_t i = 0; i < 3; ++i) {
			if (i == 0) {
				m_domain.setRange(0.01,0.45,i);
			}
			else {
				m_domain.setRange(0.01,0.1,i);
			}
		}
		
		generateAdLoadPF();
	}

	void Water_resource_plan::generateAdLoadPF() {
		m_optima.reset(new Optima<>());

		std::string file_path = g_working_dir + "/instance/problem/continuous/multi_objective/real_world/data/wrp_var_obj.dat";
		std::ifstream in_file(file_path);
		if (!in_file.fail()) {
			std::string line;
			while (std::getline(in_file, line)) {
				if (line.empty())
					continue;
				std::vector<Real> temp_var(m_number_variables);
				std::vector<Real> temp_obj(m_number_objectives);
				std::stringstream ss(line);
				for (size_t j = 0; j < m_number_variables; ++j) {
					ss >> temp_var[j];
				}
				for (size_t j = 0; j < m_number_objectives; ++j) {
					ss >> temp_obj[j];
				}
				dynamic_cast<Optima<>&>(*m_optima).appendVar(temp_var);
				dynamic_cast<Optima<>&>(*m_optima).appendObj(temp_obj);
			}
		}
		else {
			//网格采样，再评价排序
			size_t div_num = 101;
			//先二维
			std::vector<VariableVector<Real>> sample_2d;
			std::vector<std::vector<Real>> all_objs;
			for (size_t i = 0; i < div_num; i++) {
				VariableVector<Real> temp_var(2);
				temp_var[0] = m_domain.range(0).limit.first + i * (m_domain.range(0).limit.second - m_domain.range(0).limit.first) / (div_num - 1.);
				for (size_t j = 0; j < div_num; j++) {
					temp_var[1] = m_domain.range(1).limit.first + j * (m_domain.range(1).limit.second - m_domain.range(1).limit.first) / (div_num - 1.);
					sample_2d.emplace_back(temp_var);
				}
			}
			//再3维
			size_t div_num3 = 11;
			std::vector<VariableVector<Real>> sample_3d;
			for (size_t i = 0; i < sample_2d.size(); i++) {
				VariableVector<Real> temp_var(3);
				temp_var[0] = sample_2d[i][0];
				temp_var[1] = sample_2d[i][1];
				for (size_t j = 0; j < div_num3; j++) {
					temp_var[2] = m_domain.range(2).limit.first + j * (m_domain.range(2).limit.second - m_domain.range(2).limit.first) / (div_num - 1.);
					sample_3d.emplace_back(temp_var);
				}
			}

			//评价目标值
			for (size_t i = 0; i < sample_3d.size(); ++i) {
				std::vector<Real> temp_obj(m_number_objectives);
				evaluateObjective(sample_3d[i].data(), temp_obj);
				all_objs.emplace_back(temp_obj);
			}

			//先排序，再加入
			std::vector<std::vector<Real>*> objs;
			for (size_t i = 0; i < all_objs.size(); ++i) {
				objs.emplace_back(&all_objs[i]);
			}
			std::vector<int> rank;
			ofec::nd_sort::filterSortP<Real>(objs, rank, m_optimize_mode);
			for (size_t i = 0; i < rank.size(); ++i) {
				if (rank[i] == 0) {
					dynamic_cast<Optima<>&>(*m_optima).appendVar(sample_3d[i]);
					m_optima->appendObj(all_objs[i]);
				}
			}

			//std::string file_path = g_working_dir + "/instance/problem/continuous/multi_objective/real_world/data/vrw_var_obj.dat";
			std::stringstream out_stream;
			int flag_blank_line = 0;
			out_stream << std::fixed << std::setprecision(10);
			for (size_t i = 0; i < rank.size(); ++i) {
				if (rank[i] == 0) {
					for (size_t j = 0; j < m_number_variables; ++j) {
						out_stream << std::setw(10) << sample_3d[i][j] << ' ';
					}

					for (size_t j = 0; j < m_number_objectives; ++j) {
						out_stream << std::setw(10) << all_objs[i][j] << ' ';
					}
					out_stream << '\n';
				}
			}
			std::ofstream out_file(file_path);
			out_file << out_stream.str();
		}

		m_optima->setObjectiveGiven(true);
		m_optima->setVariableGiven(true);

	}

	void Water_resource_plan::evaluateObjective(Real* x, std::vector<Real>& obj) {

		obj[0] = 106780.37 * (x[1] + x[2]) + 61704.67;
		obj[1] = 3000 * x[0];
		obj[2] = 305700. * 2289. * x[1] / std::pow(0.06 * 2289, 0.65);
		obj[3] = 250 * 2289 * std::exp(-39.75 * x[1] + 9.9 * x[2] + 2.74);
		obj[4] = 25 * (1.39 / (x[0]*x[1]) + 4940 * x[2] - 80);
		
		Real g1 = 1-      0.00139  / x[0] / x[1] - 4.94     * x[2] + 0.08;
		Real g2 = 1 -     0.000306 / x[0] / x[1] - 1.082    * x[2] + 0.0986;
		Real g3 = 50000 - 12.307   / x[0] / x[1] - 49408.24 * x[2] - 4051.02;
		Real g4 = 16000 - 2.098    / x[0] / x[1] - 8046.33  * x[2] + 696.71;
		Real g5 = 10000 - 2.138    / x[0] / x[1] - 7883.39  * x[2] + 705.04;
		Real g6 = 2000 - 0.417 * x[0] * x[1] - 1721.26 * x[2] + 136.54;
		Real g7 = 550 - 0.164 / x[0] / x[1] - 631.13 * x[2] + 54.48;

		g1 = g1 >= 0 ? 0 : -1 * g1;
		g2 = g2 >= 0 ? 0 : -1 * g2;
		g3 = g3 >= 0 ? 0 : -1 * g3;
		g4 = g4 >= 0 ? 0 : -1 * g4;
		g5 = g5 >= 0 ? 0 : -1 * g5;
		g6 = g6 >= 0 ? 0 : -1 * g6;
		g7 = g7 >= 0 ? 0 : -1 * g7;

		/*g1 = g1 >= 0 ? g1 : 0;
		g2 = g2 >= 0 ? g2 : 0;
		g3 = g3 >= 0 ? g3 : 0;
		g4 = g4 >= 0 ? g4 : 0;
		g5 = g5 >= 0 ? g5 : 0;
		g6 = g6 >= 0 ? g6 : 0;
		g7 = g7 >= 0 ? g7 : 0;*/

		obj[5] = g1 + g2 + g3 + g4 + g5 + g6 + g7;
	}
}