#include "test3.h"
#include <fstream>
#include "../../../../../core/instance_manager.h"
#include "../../../../../core/global.h"
#include"../../../../../utility/nondominated_sorting/fast_sort.h"
#include <algorithm>
#include<cmath>

namespace ofec {
	/*
	2 objectives, each objective have only one peak
	*/
	void TEST3::initialize_() {
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

		m_peak_num = std::get<int>(v.at("numPeak"));
		m_max_height = 5;
		for (size_t i = 0; i < m_number_variables; ++i) {
			m_peak_norm.push_back(2);
		}
		for (size_t i = 0; i < m_peak_num; ++i) {
			m_peak_width.push_back(3);
			//m_peak_angle.push_back(OFEC_PI/3);
			m_peak_angle.push_back(0.);
		}

		bool if_multimodal = std::get<bool>(v.at("MOPmodality"));
		bool if_uniform = std::get<bool>(v.at("MOPuniform"));
		if (if_uniform) {
			size_t peak_num = std::pow(2, m_number_variables) + 1;
			std::vector<std::vector<Real>> peak_pos;
			for (size_t i = 0; i < peak_num - 1; ++i) {
				auto bin_v = number2bin(i, m_number_variables);
				std::vector<Real> temp_x;
				for (size_t j = 0; j < m_number_variables; ++j) {
					if (bin_v[j] == 1) {
						temp_x.push_back(m_domain[j].limit.first / 4 + 3 * m_domain[j].limit.second / 4);
					}
					else {
						temp_x.push_back(3 * m_domain[j].limit.first / 4 + m_domain[j].limit.second / 4);
					}
				}
				peak_pos.emplace_back(temp_x);
			}
			std::vector<Real> middle_point;
			for (size_t j = 0; j < m_number_variables; ++j) {
				middle_point.push_back(m_domain[j].limit.first / 2 + m_domain[j].limit.second / 2);
			}
			peak_pos.emplace_back(middle_point);
			m_total_peaks_pos.emplace_back(peak_pos);
		}
		else {
			size_t num = 9;
			for (size_t i = 0; i < num; ++i) {
				std::vector<std::vector<Real>> temp_peaks;
				std::vector<Real> temp_p(m_number_variables, 0.);
				switch (i) {
				case 0:
					temp_peaks.emplace_back(temp_p);
					break;
				case 1:
					temp_p[0] = 0.5; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					break;
				case 2:
					temp_p[0] = 0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					break;
				case 3:
					temp_p[0] = 0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.5; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					break;
				case 4:
					temp_p[0] = 0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.5; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0; temp_p[1] = 0;
					temp_peaks.emplace_back(temp_p);
					break;
				case 5:
					temp_p[0] = 0.5; temp_p[1] = 0;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.25; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.25; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = 0;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.25; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.25; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					break;
				case 6:
					temp_p[0] = 0.5; temp_p[1] = 0;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.25; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.25; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = 0;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.25; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.25; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.; temp_p[1] = 0.;
					temp_peaks.emplace_back(temp_p);
					break;
				case 7:
					temp_p[0] = 0.5; temp_p[1] = 0;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.5; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = 0;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					break;
				case 8:
					temp_p[0] = 0.5; temp_p[1] = 0;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.5; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = 0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = 0;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = -0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.5; temp_p[1] = -0.5;
					temp_peaks.emplace_back(temp_p);
					temp_p[0] = 0.; temp_p[1] = 0.;
					temp_peaks.emplace_back(temp_p);
					break;
				default:
					break;
				}
				m_total_peaks_pos.emplace_back(temp_peaks);
			}
		}
		for (size_t i = 0; i < m_number_objectives; ++i) {
			std::vector<std::shared_ptr<Peak>> obj_peaks;
			std::vector<std::vector<Real>> temp_peak_pos;
			for (size_t j = 0; j < m_peak_num; ++j) {
				Real h =0.;
				Real w = 0.;
				if (if_multimodal) {
					h = m_max_height;
				}
				else {
					if (j == 0) {
						h = m_max_height;
					}
					else {
						h = m_max_height*(0.7+0.3*GET_RND(m_random.get()).uniform.next());
					}
				}
				std::vector<Real> peak_x;
				if (if_uniform) {
					peak_x = m_total_peaks_pos.back()[j];
				}
				else {
					peak_x = m_total_peaks_pos[m_peak_num-1][j];
				}
				
				for (size_t k = 0; k < peak_x.size(); ++k) {
					peak_x[k] = peak_x[k] - 0.1 * i;
				}
				temp_peak_pos.emplace_back(peak_x);
				w = m_peak_width[j];
				obj_peaks.emplace_back(std::make_shared<Nonliear_peak>(peak_x,h,w,m_peak_norm,2,m_peak_angle[j]));
			}
			m_obj_mpb.emplace_back(std::make_shared<Mpb_class>(obj_peaks));
			m_peaks_pos.emplace_back(temp_peak_pos);
		}
		generateAdLoadPF();
	}

	std::vector<size_t> TEST3::number2bin(size_t v,size_t bits) {
		std::vector<size_t> flag(bits, 0);
		std::vector<size_t> bin_flag(bits, 0);
		if (v == 1) {
			bin_flag.back() = 1;
		}
		else if(v>1){
			size_t stop = 0;
			for (size_t i = 0; stop != 1; ++i) {
				flag[i] = v % 2;
				if (v == 1) {
					stop = 1;
				}
				else {
					v = v / 2;
				}
			}
			//µπ÷√
			for (size_t i = 0; i < bits; ++i) {
				bin_flag[i] = flag[bits - i - 1];
			}
		}
		return bin_flag;
	}

	void TEST3::generateAdLoadPF() {
		m_optima.reset(new Optima<>());
		size_t num = 100;
		for (size_t k = 0; k < m_peak_num; ++k) {
			auto first_peak_pos = m_peaks_pos[0][k];
			auto second_peak_pos = m_peaks_pos[1][k];
			for (size_t i = 0; i < num; i++) {
				VariableVector<Real> temp_var(m_number_variables);
				std::vector<Real> temp_obj(m_number_objectives);
				for (size_t j = 0; j < m_number_variables; j++) {
					temp_var[j] = first_peak_pos[j] + (second_peak_pos[j] - first_peak_pos[j]) * i / (num - 1.);
				}
				evaluateObjective(temp_var.data(), temp_obj);
				dynamic_cast<Optima<>&>(*m_optima).appendVar(temp_var);
				m_optima->appendObj(temp_obj);
			}
		}
		
		m_optima->setObjectiveGiven(true);
		m_optima->setVariableGiven(true);

	}

	void TEST3::evaluateObjective(Real* x, std::vector<Real>& obj) {
		std::vector<Real> temp_sol;
		for (size_t i = 0; i < m_number_variables; ++i) {
			temp_sol.push_back(x[i]);
		}

		for (size_t i = 0; i < m_number_objectives; ++i) {
			obj[i] = m_max_height-m_obj_mpb[i]->calMpbValue(temp_sol);
		}
	}
}