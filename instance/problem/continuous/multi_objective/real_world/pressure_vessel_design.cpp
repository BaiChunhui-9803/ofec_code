#include "pressure_vessel_design.h"
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
	void Pressure_vessel_design::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeObjective(2);

		for (size_t i = 0; i < m_number_objectives; ++i) {
			m_optimize_mode[i] = OptimizeMode::kMinimize;
		}
		resizeVariable(4);//n=4

		m_domain.setRange(1., 100., 0);
		m_domain.setRange(1., 100., 1);
		m_domain.setRange(10., 200., 2);
		m_domain.setRange(10., 240., 2);

		generateAdLoadPF();
	}

	void Pressure_vessel_design::generateAdLoadPF() {
		m_optima.reset(new Optima<>());

		std::string file_path = g_working_dir + "/instance/problem/continuous/multi_objective/real_world/data/pvd_var_obj.dat";
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
			size_t div_num1 = 21;
			size_t div_num2 = 21;
			size_t div_num3 = 21;
			size_t div_num4 = 21;

			//先二维
			std::vector<VariableVector<Real>> sample_2d;
			std::vector<std::vector<Real>> all_objs;
			for (size_t i = 0; i < div_num1; i++) {
				VariableVector<Real> temp_var(2);
				temp_var[0] = m_domain.range(0).limit.first + i * (m_domain.range(0).limit.second - m_domain.range(0).limit.first) / (div_num1 - 1.);
				for (size_t j = 0; j < div_num2; j++) {
					temp_var[1] = m_domain.range(1).limit.first + j * (m_domain.range(1).limit.second - m_domain.range(1).limit.first) / (div_num2 - 1.);
					sample_2d.emplace_back(temp_var);
				}
			}

			//再3维
			std::vector<VariableVector<Real>> sample_3d;
			for (size_t i = 0; i < sample_2d.size(); i++) {
				VariableVector<Real> temp_var(3);
				temp_var[0] = sample_2d[i][0];
				temp_var[1] = sample_2d[i][1];
				for (size_t j = 0; j < div_num3; j++) {
					temp_var[2] = m_domain.range(2).limit.first + j * (m_domain.range(2).limit.second - m_domain.range(2).limit.first) / (div_num3 - 1.);
					sample_3d.emplace_back(temp_var);
				}
			}

			//再4维
			std::vector<VariableVector<Real>> sample_4d;
			for (size_t i = 0; i < sample_3d.size(); i++) {
				VariableVector<Real> temp_var(4);
				temp_var[0] = sample_3d[i][0];
				temp_var[1] = sample_3d[i][1];
				temp_var[2] = sample_3d[i][2];
				for (size_t j = 0; j < div_num4; j++) {
					temp_var[3] = m_domain.range(3).limit.first + j * (m_domain.range(3).limit.second - m_domain.range(3).limit.first) / (div_num4 - 1.);
					sample_4d.emplace_back(temp_var);
				}
			}

			//评价目标值
			for (size_t i = 0; i < sample_4d.size(); ++i) {
				std::vector<Real> temp_obj(m_number_objectives);
				evaluateObjective(sample_4d[i].data(), temp_obj);
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
					dynamic_cast<Optima<>&>(*m_optima).appendVar(sample_4d[i]);
					m_optima->appendObj(all_objs[i]);
				}
			}

			std::stringstream out_stream;
			int flag_blank_line = 0;
			out_stream << std::fixed << std::setprecision(6);
			for (size_t i = 0; i < rank.size(); ++i) {
				if (rank[i] == 0) {
					for (size_t j = 0; j < m_number_variables; ++j) {
						out_stream << std::setw(6) << sample_4d[i][j] << ' ';
					}

					for (size_t j = 0; j < m_number_objectives; ++j) {
						out_stream << std::setw(6) << all_objs[i][j] << ' ';
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

	void Pressure_vessel_design::evaluateObjective(Real* x, std::vector<Real>& obj) {

		obj[0] = 0.6224*x[0]*x[2]*x[3]+1.7781*x[1]*std::pow(x[2],2)+3.1661*std::pow(x[0],2)*x[3]+19.84*std::pow(x[0],2) * x[2];


		Real g1 = x[0]-0.0193*x[2];
		Real g2 = x[1]-0.00954*x[2];
		Real g3 = OFEC_PI * std::pow(x[2], 2) * x[3] + 4 * OFEC_PI / 3 * std::pow(x[2], 3) - 1296000;

		g1 = g1 >= 0 ? 0 : -1 * g1;
		g2 = g2 >= 0 ? 0 : -1 * g2;
		g3 = g3 >= 0 ? 0 : -1 * g3;


		/*g1 = g1 >= 0 ? g1 : 0;
		g2 = g2 >= 0 ? g2 : 0;
		g3 = g3 >= 0 ? g3 : 0;*/

		obj[1] = g1 + g2 +g3;
	}
}