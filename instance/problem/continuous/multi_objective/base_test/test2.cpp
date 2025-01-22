#include "test2.h"
#include <fstream>
#include "../../../../../core/instance_manager.h"
#include "../../../../../core/global.h"
#include"../../../../../utility/nondominated_sorting/fast_sort.h"
#include <algorithm>

namespace ofec {
	/*
	2 objectives, each objective have only one peak
	*/
	void TEST2::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeObjective(std::get<int>(v.at("number of objectives")));
		if (m_number_objectives != 2)
			throw MyExcept("The number of objectives must be equal to 2");
		for (size_t i = 0; i < m_number_objectives; ++i) {
			m_optimize_mode[i] = OptimizeMode::kMinimize;
		}

		resizeVariable(v.get<int>("number of variables"));//recomend n=10
		if (m_number_variables < 2)
			throw MyExcept("The number of variables must be no less than 2");
		setDomain(-1., 1.);

		//initial the location of the peaks for each objective
		m_num_PS = std::get<int>(v.at("number of PS shapes"));
		m_radius.push_back(std::get<Real>(v.at("radiusPS")));
		for (size_t i = 0; i < m_num_PS-1; ++i) {
			m_radius.push_back(m_radius[0]*(1+ GET_RND(m_random.get()).uniform.next()));
		}
		m_rand_flag = false;
		if (m_rand_flag) {
			for (size_t i = 0; i < m_num_PS; ++i) {
				std::vector<Real> temp_center;
				for (size_t j = 0; j < m_number_variables; ++j) {
					temp_center.push_back((domain().range().limit.second - m_radius[i]) *(2* GET_RND(m_random.get()).uniform.next() - 1));
				}
				m_center_pos.push_back(temp_center);
			}
		}
		else {
			//每个象限初始化一个PS
			std::vector<std::vector<Real>> all_centers;
			std::vector<Real> temp_location = {-0.5,0.5};
			std::vector<Real> optima_location;
			for (size_t j = 0; j < m_number_variables; ++j) {
				optima_location.push_back(0);
			}
			all_centers.emplace_back(optima_location);
			while (all_centers.size()<std::pow(2,(m_number_variables>8?8:m_number_variables))) {
				std::vector<Real> temp;
				for (size_t j = 0; j < m_number_variables; ++j) {
					auto inx = (size_t)floor(temp_location.size() * GET_RND(m_random.get()).uniform.next());
					temp.push_back(temp_location[inx]);
				}
				bool flag = false;
				for (size_t i = 0; i < all_centers.size(); ++i) {
					auto tt = all_centers[i];
					bool gg = false;
					Real diff_sum = 0;
					for (size_t j = 0; j < tt.size(); ++j) {
						diff_sum += (std::fabs(temp[j]-tt[j]));
					}
					if (diff_sum < 0.5) {
						flag = true;
						break;
					}
				}
				if (!flag) {
					all_centers.emplace_back(temp);
				}
			}
			for (size_t i = 0; i < m_num_PS; ++i) {
				m_center_pos.emplace_back(all_centers[i]);
			}
			
			/*std::vector<Real> center1 = {3,0};
			std::vector<Real> center2 = { -2,3 };
			std::vector<Real> center3 = { -2,-3 };
			m_center_pos.push_back(center1);
			m_center_pos.push_back(center2);
			m_center_pos.push_back(center3);*/
		}
		for (size_t i = 0; i < m_number_objectives; ++i) {
			std::vector<Real> temp_angles;
			if (m_rand_flag) {
				for (size_t j = 0; j < m_number_variables - 1; ++j) {
					temp_angles.push_back(2 * OFEC_PI * GET_RND(m_random.get()).uniform.next());
				}
			}
			else {
				for (size_t j = 0; j < m_number_variables - 1; ++j) {
					temp_angles.push_back(i * 2 * OFEC_PI / m_number_objectives + OFEC_PI / 3);
				}
			}
			m_angles.emplace_back(temp_angles);
		}
		std::vector<std::vector<Real>> obj_optimal_peak_pos;//global ps location
		for (size_t i = 0; i < m_number_objectives; ++i) {
			std::vector<Real> optimal_peak_pos;
			for (size_t j = 0; j < m_number_variables; ++j) {
				Real temp = 1;
				for (size_t k = 0; k < m_angles[i].size() - j; ++k) {
					temp *= std::cos(m_angles[i][k]);
				}
				if (j == 0) {
					optimal_peak_pos.push_back(m_center_pos[0][j] + m_radius[0]*temp);
				}
				else {
					optimal_peak_pos.push_back(m_center_pos[0][j] + m_radius[0] * temp * std::sin(m_angles[i][m_angles[i].size() - j]));
				}
			}
			obj_optimal_peak_pos.emplace_back(optimal_peak_pos);
		}
		m_peaks_pos.emplace_back(obj_optimal_peak_pos);
		for (size_t i = 0; i < m_num_PS - 1; ++i) {
			std::vector<std::vector<Real>> peak_pos = obj_optimal_peak_pos;
			for (size_t j = 0; j < peak_pos.size(); ++j) {
				for (size_t k = 0; k < m_number_variables; ++k) {
					peak_pos[j][k] = m_center_pos[0][k] + (obj_optimal_peak_pos[j][k] - m_center_pos[0][k]) * (m_radius[i+1] / m_radius[0]);
				}
			}
			for (size_t j = 0; j < m_number_variables; ++j) {
				for (size_t k = 0; k < peak_pos.size(); ++k) {
					peak_pos[k][j] = peak_pos[k][j] + m_center_pos[i+1][j] - m_center_pos[0][j];
				}
			}
			m_peaks_pos.emplace_back(peak_pos);
		}
		
		//initial heights and slopes
		m_heights = { 10,10 };
		m_slopes = { 5,5 };

		generateAdLoadPF();
	}

	void TEST2::generateAdLoadPF() {
		size_t num = 2000;
		m_optima.reset(new Optima<>());
		for (size_t i = 0; i < num; i++) {
			VariableVector<Real> temp_var(m_number_variables);
			std::vector<Real> temp_obj(m_number_objectives);
			for (size_t j = 0; j < m_number_variables; j++) {
				temp_var[j] = m_peaks_pos[0][0][j] + (m_peaks_pos[0][1][j] - m_peaks_pos[0][0][j]) * i / (num - 1.);
			}
			evaluateObjective(temp_var.data(), temp_obj);
			dynamic_cast<Optima<>&>(*m_optima).appendVar(temp_var);
			m_optima->appendObj(temp_obj);
		}
		m_optima->setObjectiveGiven(true);
		m_optima->setVariableGiven(true);
	}

	void TEST2::evaluateObjective(Real* x, std::vector<Real>& obj) {
		std::vector<Real> temp_sol;
		for (size_t i = 0; i < m_number_variables; ++i) {
			temp_sol.push_back(x[i]);
		}
		std::vector<std::vector<Real>> temp_dist;
		for (size_t j = 0; j < m_number_objectives; ++j) {
			std::vector<Real> temp_ps_dist;
			for (size_t i = 0; i < m_num_PS; ++i) {
				temp_ps_dist.push_back(euclideanDistance(m_peaks_pos[i][j].begin(), m_peaks_pos[i][j].end(), temp_sol.begin()));
			}
			temp_dist.push_back(temp_ps_dist);
		}
		
		std::vector<Real> min_dist;
		for (size_t i = 0; i < temp_dist.size(); ++i) {
			min_dist.push_back(*std::min_element(temp_dist[i].begin(), temp_dist[i].end()));
		}

		for (size_t i = 0; i < m_number_objectives; ++i) {
			//obj[i] = m_heights[i] * (m_slopes[i] * min_dist[i]) / (1 + m_slopes[i] * min_dist[i]);
			obj[i] = m_slopes[i] * min_dist[i];
		}
	}
}