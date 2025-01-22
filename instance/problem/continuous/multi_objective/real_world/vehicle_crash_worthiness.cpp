#include "vehicle_crash_worthiness.h"
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
	void Vehicle_crash_worthiness::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeObjective(3);
		
		for (size_t i = 0; i < m_number_objectives; ++i) {
			m_optimize_mode[i] = OptimizeMode::kMinimize;
		}
		resizeVariable(5);//n=5
		
		for (size_t i = 0; i < 5; ++i) {
			m_domain.setRange(1., 3., i);
		}

		generateAdLoadPF();
	}

	void Vehicle_crash_worthiness::generateAdLoadPF() {
		m_optima.reset(new Optima<>());

		std::string file_path = g_working_dir + "/instance/problem/continuous/multi_objective/real_world/data/vcw_cluster_obj.dat";
		std::ifstream in_file(file_path);
		if (!in_file.fail()) {
			std::string line;
			while (std::getline(in_file, line)) {
				if (line.empty())
					continue;
				std::vector<Real> temp_obj(m_number_objectives);
				std::stringstream ss(line);
				
				for (size_t j = 0; j < m_number_objectives; ++j) {
					ss >> temp_obj[j];
				}
				m_optima->appendObj(temp_obj);
			}
		}

		std::string file_path2 = g_working_dir + "/instance/problem/continuous/multi_objective/real_world/data/vcw_cluster_var.dat";
		std::ifstream in_file2(file_path2);
		if (!in_file2.fail()) {
			std::string line;
			while (std::getline(in_file2, line)) {
				if (line.empty())
					continue;
				std::vector<Real> temp_var(m_number_variables);
				std::stringstream ss(line);
				for (size_t j = 0; j < m_number_variables; ++j) {
					ss >> temp_var[j];
				}
				dynamic_cast<Optima<>&>(*m_optima).appendVar(temp_var);
			}
		}

		


		//std::string file_path = g_working_dir + "/instance/problem/continuous/multi_objective/real_world/data/vcw_var_obj.dat";
		//std::ifstream in_file(file_path);
		//if (!in_file.fail()) {
		//	std::string line;
		//	while (std::getline(in_file, line)) {
		//		if (line.empty())
		//			continue;
		//		std::vector<Real> temp_var(m_number_variables);
		//		std::vector<Real> temp_obj(m_number_objectives);
		//		std::stringstream ss(line);
		//		for (size_t j = 0; j < m_number_variables; ++j) {
		//			ss >> temp_var[j];
		//		}
		//		for (size_t j = 0; j < m_number_objectives; ++j) {
		//			ss >> temp_obj[j];
		//		}
		//		dynamic_cast<Optima<>&>(*m_optima).appendVar(temp_var);
		//		dynamic_cast<Optima<>&>(*m_optima).appendObj(temp_obj);
		//	}
		//}
		//else {
		//	//网格采样，再评价排序
		//	size_t div_num = 11;
		//	//先二维
		//	std::vector<VariableVector<Real>> sample_2d;
		//	std::vector<std::vector<Real>> all_objs;
		//	for (size_t i = 0; i < div_num; i++) {
		//		VariableVector<Real> temp_var(2);
		//		temp_var[0] = m_domain.range(0).limit.first + i * (m_domain.range(0).limit.second- m_domain.range(0).limit.first) / (div_num-1.);
		//		for (size_t j = 0; j < div_num; j++) {
		//			temp_var[1] = m_domain.range(1).limit.first + j * (m_domain.range(1).limit.second - m_domain.range(1).limit.first) / (div_num - 1.);
		//			sample_2d.emplace_back(temp_var);
		//		}
		//	}
		//	//再4维
		//	std::vector<VariableVector<Real>> sample_4d;
		//	for (size_t i = 0; i < sample_2d.size(); i++) {
		//		VariableVector<Real> temp_var(4);
		//		temp_var[0] = sample_2d[i][0];
		//		temp_var[1] = sample_2d[i][1];
		//		for (size_t j = 0; j < sample_2d.size(); j++) {
		//			temp_var[2] = sample_2d[j][0];
		//			temp_var[3] = sample_2d[j][1];
		//			sample_4d.emplace_back(temp_var);
		//		}
		//	}
		//	//最后得到5d
		//	std::vector<VariableVector<Real>> sample_5d;
		//	for (size_t i = 0; i < div_num; i++) {
		//		VariableVector<Real> temp_var(5);
		//		for (size_t j = 0; j < sample_4d.size(); j++) {
		//			temp_var[0] = sample_4d[j][0];
		//			temp_var[1] = sample_4d[j][1];
		//			temp_var[2] = sample_4d[j][2];
		//			temp_var[3] = sample_4d[j][3];
		//			temp_var[4] = m_domain.range(4).limit.first + i * (m_domain.range(4).limit.second - m_domain.range(4).limit.first) / (div_num - 1.);
		//			sample_5d.emplace_back(temp_var);
		//		}
		//	}

		//	//评价目标值
		//	for (size_t i = 0; i < sample_5d.size(); ++i) {
		//		std::vector<Real> temp_obj(m_number_objectives);
		//		evaluateObjective(sample_5d[i].data(), temp_obj);
		//		all_objs.emplace_back(temp_obj);
		//	}

		//	//先排序，再加入
		//	std::vector<std::vector<Real>*> objs;
		//	for (size_t i = 0; i < all_objs.size(); ++i) {
		//		objs.emplace_back(&all_objs[i]);
		//	}
		//	std::vector<int> rank;
		//	ofec::nd_sort::filterSortP<Real>(objs, rank, m_optimize_mode);
		//	for (size_t i = 0; i < rank.size(); ++i) {
		//		if (rank[i] == 0) {
		//			dynamic_cast<Optima<>&>(*m_optima).appendVar(sample_5d[i]);
		//			m_optima->appendObj(all_objs[i]);
		//		}
		//	}

		//	//std::string file_path = g_working_dir + "/instance/problem/continuous/multi_objective/real_world/data/vrw_var_obj.dat";
		//	std::stringstream out_stream;
		//	int flag_blank_line = 0;
		//	out_stream << std::fixed << std::setprecision(3);
		//	for (size_t i = 0; i < rank.size(); ++i) {
		//		if (rank[i] == 0) {
		//			for (size_t j = 0; j < m_number_variables; ++j) {
		//				out_stream << std::setw(10) << sample_5d[i][j] << ' ';
		//			}

		//			for (size_t j = 0; j < m_number_objectives; ++j) {
		//				out_stream << std::setw(10) << all_objs[i][j] << ' ';
		//			}
		//			out_stream << '\n';
		//		}
		//	}
		//	std::ofstream out_file(file_path);
		//	out_file << out_stream.str();
		//}

		m_optima->setObjectiveGiven(true);
		m_optima->setVariableGiven(true);

	}

	void Vehicle_crash_worthiness::evaluateObjective(Real* x, std::vector<Real>& obj) {

		obj[0] = 1640.2823 + 2.3573285 * x[0] + 2.3220035 * x[1] + 4.5688768 * x[2] + \
			7.7213633 * x[3] + 4.4559504 * x[4];
		obj[1] = 6.5856 + 1.15 * x[0] - 1.0427 * x[1] + 0.9738 * x[2] + 0.8364 * x[3] - \
			0.3695 * x[0]*x[3] + 0.0861 * x[0]*x[4] + 0.3628 * x[1] * x[3] - \
			0.1106 * x[0]* x[0] - 0.3437 * x[2]* x[2] + 0.1764 * x[3]* x[3];
		obj[2] = -0.0551 + 0.0181 * x[0] + 0.1024 * x[1] + 0.0421 * x[2] - 0.0073 * x[0] * x[1] + \
			0.024 * x[1] * x[2] - 0.0118 * x[1] * x[3] - 0.0204 * x[2] * x[3] - 0.008 * x[2] * x[4] - \
			0.0241 * x[1]* x[1] + 0.0109 * x[3]* x[3];
	}
}