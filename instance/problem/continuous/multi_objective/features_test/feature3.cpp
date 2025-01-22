#include "feature3.h"
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
	void Feature3::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;
		resizeObjective(std::get<int>(v.at("number of objectives")));
		if (m_number_objectives != 2)
			throw MyExcept("The number of objectives must be equal to 2");
		for (size_t i = 0; i < m_number_objectives; ++i) {
			m_optimize_mode[i] = OptimizeMode::kMinimize;
		}
		resizeVariable(std::get<int>(v.at("number of variables")));//recomend n=10
		if (m_number_variables < 2)
			throw MyExcept("The number of variables must be no less than 2");
		setDomain(-1., 1.);

		m_peak_num = 2;
		m_max_height = 5;
		for (size_t i = 0; i < m_number_variables; ++i) {
			m_peak_norm.push_back(2);
		}
		for (size_t i = 0; i < m_peak_num; ++i) {
			m_peak_width.push_back(3);
			//m_peak_angle.push_back(OFEC_PI/3);
			m_peak_angle.push_back(0.);
		}

		size_t num = m_peak_num;
		for (size_t i = 0; i < num; ++i) {
			std::vector<std::vector<Real>> temp_peaks;//第一个目标峰的位置
			std::vector<Real> temp_p(m_number_variables, 0.);
			if (i == 0) {
				for (size_t j = 0; j < m_number_variables; ++j) {
					if (j == m_number_variables - 2) {
						temp_p[j] = 1.;
					}
					else {
						temp_p[j] = 1. * GET_RND(m_random.get()).uniform.next() - 0.5;//位置随机
					}
				}
			}
			else {
				for (size_t j = 0; j < m_number_variables; ++j) {
					if (j == m_number_variables - 1) {
						temp_p[j] = 1.;
					}
					else {
						temp_p[j] = 1. * GET_RND(m_random.get()).uniform.next() - 0.5;//位置随机
					}
				}
			}
			temp_peaks.emplace_back(temp_p);
			m_total_peaks_pos.emplace_back(temp_peaks);
		}

		
		for (size_t i = 0; i < m_number_objectives; ++i) {
			std::vector<std::shared_ptr<Peak>> obj_peaks;
			std::vector<std::vector<Real>> temp_peak_pos;
			for (size_t j = 0; j < m_peak_num; ++j) {
				Real h = m_max_height * (0.7 + 0.3 * GET_RND(m_random.get()).uniform.next());
				Real w = 0.;
				std::vector<Real> peak_x;
				peak_x = m_total_peaks_pos[j][0];
				
				if (i > 0) {
					if (j == 0) {
						for (size_t k = 1; k < peak_x.size(); ++k) {
							peak_x[k] = peak_x[k] - 0.1 * i;
						}
					}
					else {
						for (size_t k = 0; k < peak_x.size()-1; ++k) {
							peak_x[k] = peak_x[k] - 0.1 * i;
						}
					}
				}

				temp_peak_pos.emplace_back(peak_x);
				w = m_peak_width[j];
				obj_peaks.emplace_back(std::make_shared<Nonliear_peak>(peak_x, h, w, m_peak_norm, 2, m_peak_angle[j]));
			}
			m_obj_mpb.emplace_back(std::make_shared<Mpb_class>(obj_peaks));
			m_peaks_pos.emplace_back(temp_peak_pos);
		}
		generateAdLoadPF();
	}

	void Feature3::generateAdLoadPF() {
		size_t num = 500;
		m_optima.reset(new Optima<>());
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

	void Feature3::evaluateObjective(Real* x, std::vector<Real>& obj) {
		std::vector<Real> temp_sol;
		for (size_t i = 0; i < m_number_variables; ++i) {
			temp_sol.push_back(x[i]);
		}

		for (size_t i = 0; i < m_number_objectives; ++i) {
			obj[i] = m_max_height - m_obj_mpb[i]->calMpbValue(temp_sol);
		}
	}
}