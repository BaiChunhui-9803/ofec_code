/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

#include "f10_composition2.h"
#include "f6_hybrid1.h"
#include "f7_hybrid2.h"
#include "f8_hybrid3.h"
#include "../../../../../../core/global.h"

namespace ofec {
	namespace cec2015 {
		F10_composition2::~F10_composition2() {
			for (auto &i : m_hybrid)
				if (i) delete i;
		}

		void F10_composition2::setFunction() {
			m_num_hybrid = 3;
			m_hybrid.resize(m_num_hybrid);
			m_height.resize(m_num_hybrid);
			m_converge_severity.resize(m_num_hybrid);
			m_f_bias.resize(m_num_hybrid);
			m_hy_bias.resize(m_num_hybrid);

			BasicFunctions f(m_num_hybrid);
			f[0] = &createFunction<F6_hybrid1>;
			f[1] = &createFunction<F7_hybrid2>;
			f[2] = &createFunction<F8_hybrid3>;

			for (size_t i = 0; i < m_num_function; ++i) {
				auto param = *m_param;
				param["problem name"] = m_name + "_part" + std::to_string(i + 1);
				m_param_fun[i] = std::make_shared<const ParameterMap>(param);
				m_hybrid[i] = dynamic_cast<Hybrid*>(f[i]());
				m_hybrid[i]->setIdRnd(m_random.get());
				m_hybrid[i]->setIdParam(m_param_fun[i]);
				m_hybrid[i]->initialize();
			}

			m_converge_severity[0] = 10;
			m_converge_severity[1] = 30;
			m_converge_severity[2] = 50;

			m_height[0] = 1;
			m_height[1] = 1;
			m_height[2] = 1;

			m_f_bias[0] = 0;
			m_f_bias[1] = 100;
			m_f_bias[2] = 200;

			m_hy_bias[0] = 600;
			m_hy_bias[1] = 700;
			m_hy_bias[2] = 800;
		}

		void F10_composition2::initialize_(Environment *env) {
			Continuous::initialize_();
			resizeObjective(1);
			m_optimize_mode[0] = OptimizeMode::kMinimize;

			auto& v = *m_param;;
			resizeVariable(std::get<int>(m_param->at("number of variables")));

//			resizeVariable(v.get<int>("number of variables"));
			setDomain(-5., 5.);

			for (auto &i : m_hybrid)
				if (i) delete i;
			for (auto id : m_param_fun)
				DEL_PARAM(id);

			setFunction();
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");  //data path
			for (size_t i = 0; i < m_hybrid.size(); ++i) {
				for (size_t j = 0; j < m_number_variables; ++j)
					s.variable()[j] = i->get_optima().variable(0)[j] + i->hybrid_translation()[j];
				i->get_optima().clear();
				i->get_optima().append(s.variable());
				i->evaluate_(s, caller::Problem, false, false);
				s.objective()[0] -= m_hy_bias[i];
				i->get_optima().append(s.objective());
			}
			// Set optimal solution
			for (size_t j = 0; j < m_number_variables; ++j)
				s.variable()[j] = m_hybrid[0]->get_optima().variable(0)[j];

			m_optima->append(s.variable());
			s.evaluate(false, caller::Problem);
			m_optima->append(s.objective());
			m_initialized = true;
		}

		int F10_composition2::evaluateObjective(Real *x, std::vector<Real> &obj) {
			std::vector<Real> x_(m_number_variables);
			std::copy(x, x + m_number_variables, x_.begin());
			std::vector<Real> weight(m_num_function, 0);

			set_weight(weight, x_);
			std::vector<Real> fit(m_num_function);
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
			std::vector<Real> hy_bias = { 600,700,800 };
			for (size_t i = 0; i < m_num_function; ++i) { // calculate objective value for each function
				s.variable().vect() = x_;
				for (size_t j = 0; j < m_number_variables; ++j)
					s.variable()[j] -= m_hybrid[i]->hybrid_translation()[j];
				m_hybrid[i]->evaluate_(s, caller::Problem, false, false);
				fit[i] = s.objective()[0] - hy_bias[i];

			}
			Real sumw = 0;
			for (size_t i = 0; i < m_num_function; ++i)
				sumw += weight[i];
			for (size_t i = 0; i < m_num_function; ++i)
				weight[i] /= sumw;

			Real temp = 0;
			for (size_t i = 0; i < m_num_function; ++i) {
				temp += weight[i] * (m_height[i] * fit[i] + m_f_bias[i]);
			}

			obj[0] = temp;
			obj[0] += 1000;
			return kNormalEval;
		}

		bool F10_composition2::loadTranslation(const std::string &path) {
			std::string s;
			std::stringstream ss;
			ss << m_number_variables << "Dim.txt";
			s = ss.str();
			s.insert(0, m_name + "_Shift_");
			s.insert(0, path);    // data path
			s.insert(0, g_working_directory);

			for (auto &i : m_hybrid)
				i->hybrid_translation().resize(m_number_variables);
			std::ifstream in(s);
			if (in.fail()) {
				setTranslation();
				std::ofstream out(s);
				for (size_t i = 0; i < m_num_function; ++i)
					for (size_t j = 0; j < m_number_variables; ++j)
						out << m_hybrid[i]->hybrid_translation()[j] << " ";

				out.close();
			}
			else {
				for (size_t i = 0; i < m_num_function; ++i)
					for (size_t j = 0; j < m_number_variables; ++j)
						in >> m_hybrid[i]->hybrid_translation()[j];
			}
			in.close();
			return true;
		}

		void F10_composition2::setTranslation() {
			for (int i = 0; i < m_num_function; i++)
				for (int j = 0; j < m_number_variables; j++)
					m_hybrid[i]->hybrid_translation()[j] = m_domain[j].limit.first + (m_domain[j].limit.second - m_domain[j].limit.first)*(1 - global::ms_global->m_uniform[caller::Problem]->next());
		}

		void F10_composition2::setWeight(std::vector<Real>& weight, const std::vector<Real>&x) {

			for (size_t i = 0; i < m_num_function; ++i) { // calculate weight for each function
				weight[i] = 0;
				for (size_t j = 0; j < m_number_variables; ++j) {
					//weight[i] += pow(x[j] - m_function[i]->translation()[j], 2);
					weight[i] += pow(x[j] - m_hybrid[i]->get_optima().variable(0)[j], 2);
				}

				if (fabs(weight[i])>1e-6) weight[i] = exp(-weight[i] / (2 * m_number_variables*m_converge_severity[i] * m_converge_severity[i])) / (Real)(pow(weight[i], 0.5));
				else {
					for (auto &m : weight) {
						m = 0;
					}
					weight[i] = 1;
					break;
				}
			}
		}

		Hybrid* F10_composition2::getHybrid(size_t num){
			return m_hybrid[num]; 
		}
	}
}


