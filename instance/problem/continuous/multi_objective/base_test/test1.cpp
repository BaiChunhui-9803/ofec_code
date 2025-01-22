#include "test1.h"
#include <fstream>
#include "../../../../../core/instance_manager.h"
#include "../../../../../core/global.h"
#include"../../../../../utility/nondominated_sorting/fast_sort.h"

namespace ofec {
	/*
	2 objectives, each objective have only one peak
	*/
	void TEST1::initialize_() {
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
		m_radius = 0.4;
		m_rand_flag = false;
		if (m_rand_flag) {
			for (size_t j = 0; j < m_number_variables; ++j) {
				m_center_pos.push_back((2 * GET_RND(m_random.get()).uniform.next() - 1)*(1-m_radius));
			}
		}
		else{
			for (size_t j = 0; j < m_number_variables; ++j) {
				m_center_pos.push_back(0.1);
			}
		}
		for (size_t i = 0; i < m_number_objectives; ++i) {
			std::vector<Real> temp_angles;
			if (m_rand_flag) {
				for (size_t j = 0; j < m_number_variables - 1; ++j) {
					temp_angles.push_back(2 * OFEC_PI *GET_RND(m_random.get()).uniform.next());
				}
			}
			else {
				for (size_t j = 0; j < m_number_variables-1; ++j) {
					temp_angles.push_back(i*2*OFEC_PI/m_number_objectives + OFEC_PI / 3);
				}
			}
			m_angles.emplace_back(temp_angles);
		}
		for (size_t i = 0; i < m_number_objectives; ++i) {
			std::vector<Real> peak_pos;
			for (size_t j = 0; j < m_number_variables; ++j) {
				Real temp = 1;
				for (size_t k = 0; k < m_angles[i].size() - j; ++k) {
					temp *= std::cos(m_angles[i][k]);
				}
				if (j == 0) {
					peak_pos.push_back(m_center_pos[j] + m_radius*temp);
				}
				else {
					peak_pos.push_back(m_center_pos[j] + m_radius*temp * std::sin(m_angles[i][m_angles[i].size() - j]));
				}
			}
			m_peaks_pos.emplace_back(peak_pos);
		}
		//initial heights and slopes
		m_heights = {10,10};
		m_slopes = {10,10};

		generateAdLoadPF();
	}

	void TEST1::generateAdLoadPF() {
		size_t num = 2000;
		m_optima.reset(new Optima<>());
		for (size_t i = 0; i < num; i++) {
			VariableVector<Real> temp_var(m_number_variables);
			std::vector<Real> temp_obj(m_number_objectives);
			for (size_t j = 0; j < m_number_variables; j++) {
				temp_var[j] = m_peaks_pos[0][j] + (m_peaks_pos[1][j]-m_peaks_pos[0][j])*i/(num-1.);
			}
			evaluateObjective(temp_var.data(), temp_obj);
			dynamic_cast<Optima<>&>(*m_optima).appendVar(temp_var);
			m_optima->appendObj(temp_obj);
		}
		m_optima->setObjectiveGiven(true);
		m_optima->setVariableGiven(true);
	}

	void TEST1::evaluateObjective(Real* x, std::vector<Real>& obj) {
		std::vector<Real> temp_sol;
		for (size_t i = 0; i < m_number_variables; ++i) {
			temp_sol.push_back(x[i]);
		}
		std::vector<Real> temp_dist;
		for (size_t i = 0; i < m_number_objectives; ++i) {
			temp_dist.push_back(euclideanDistance(m_peaks_pos[i].begin(), m_peaks_pos[i].end(), temp_sol.begin()));
		}
		for (size_t i = 0; i < m_number_objectives; ++i) {
			//nonlinear
			//obj[i] = m_heights[i] * (m_slopes[i]*temp_dist[i])/(1+ m_slopes[i] * temp_dist[i]);
			//linear
			obj[i] = m_slopes[i] * temp_dist[i];
		}
	}
}