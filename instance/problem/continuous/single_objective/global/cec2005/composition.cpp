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

#include <algorithm>
#include "composition.h"
#include "../../../../../../core/problem/solution.h"
#include "../../../../../../core/global.h"

namespace ofec {
	namespace cec2005 {
		size_t Composition::numFunctions() {
			return m_num_function;
		}

		void Composition::computeFmax() {  // calculate the estimate max value of funciton i
			VariableVector<> vars(m_number_variables);
			std::vector<Real> objs(m_number_objectives), cons(m_number_constraints);
			for (size_t i = 0; i < m_num_function; ++i) {
				for (size_t j = 0; j < m_number_variables; ++j) {
					vars[j] = m_domain[j].limit.second;
				}
				m_function[i]->evaluate(vars, objs, cons);
				m_fmax[i] = objs[0];
			}
		}

		void Composition::setWeight(std::vector<Real>& weight, const std::vector<Real>&x) const { //default CEC05
			for (size_t i = 0; i < m_num_function; ++i) { // calculate weight for each function
				weight[i] = 0;
				for (size_t j = 0; j < m_number_variables; ++j) {
					weight[i] += (x[j] - m_function[i]->translation()[j])*(x[j] - m_function[i]->translation()[j]);
				}
				weight[i] = exp(-weight[i] / (2 * m_number_variables*m_converge_severity[i] * m_converge_severity[i]));

			}
		}

		void Composition::addInputParameters() {
			m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
		}

		void Composition::initialize_(Environment *env) {
			Continuous::initialize_(env);
			m_number_objectives = 1;
			m_optimize_mode.resize(m_number_objectives);
			m_optimize_mode[0] = OptimizeMode::kMinimize;
			m_domain.resize(m_number_variables);
			setDomain(-5., 5.);
			m_height_normalize_severity = 2000;
			setFunction(env);
			loadRotation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005");
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2005/GOP_CEC2005");  //data path
			computeFmax();
		}

		void Composition::updateOptima(Environment *env) {
			m_optima.reset(new Optima<>());
			Solution<> sol(m_number_objectives, m_number_constraints, m_number_variables);
			sol.variable() = m_function[0]->translation();
			evaluate(sol.variableBase(), sol.objective(), sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(sol);
		}

		void Composition::evaluateObjective(Real *x, std::vector<Real> &obj) const {
			std::vector<Real> x_(m_number_variables);
			std::copy(x, x + m_number_variables, x_.begin());
			
			std::vector<Real> weight(m_num_function, 0);
			setWeight(weight, x_);

			std::vector<Real> fit(m_num_function);
			Solution<> s(m_number_objectives, m_number_constraints, m_number_variables);
			s.variable().vector() = x_;

			for (size_t i = 0; i < m_num_function; ++i) { // calculate objective value for each function
				m_function[i]->evaluate(s.variableBase(), s.objective(), s.constraint());
				fit[i] = s.objective()[0];
				if (fabs(m_fmax[i]) > 1e-6)
					fit[i] = m_height_normalize_severity*fit[i] / fabs(m_fmax[i]);
			}

			Real sumw = 0, wmax;
			wmax = *std::max_element(weight.begin(), weight.end());
			for (size_t i = 0; i < m_num_function; ++i) {
				if (weight[i] != wmax) {
					weight[i] = weight[i] * (1 - pow(wmax, 10));
				}
			}
			size_t same_wmax_num = 0;
			for (size_t i = 0; i < m_num_function; ++i) {
				if (weight[i] == wmax) ++same_wmax_num;
			}
			size_t i = m_num_function - 1;
			while (same_wmax_num > 1 && i >= 0) {
				if (wmax == weight[i]) {
					weight[i] = 0;
					--same_wmax_num;
				}
				--i;
			}

			for (size_t i = 0; i < m_num_function; ++i)
				sumw += weight[i];
			for (size_t i = 0; i < m_num_function; ++i)
				weight[i] /= sumw;

			Real temp = 0;
			for (size_t i = 0; i < m_num_function; ++i) {
				temp += weight[i] * (fit[i] + m_height[i]);
			}

			obj[0] = temp;
		}

		bool Composition::loadTranslation(const std::string &path) {
			std::string s;
			std::stringstream ss;
			ss << m_number_variables << "Dim.txt";
			s = ss.str();
			s.insert(0, "_F" + m_index_number + "_Shift_");
			s.insert(0, path);    // data path
			s.insert(0, g_working_directory);
			for (auto &i : m_function)
				i->translation().resize(m_number_variables);
			std::ifstream in(s);
			if (in.fail()) {
				setTranslation();
				std::ofstream out(s);
				out << "# The data is generated by OFEC randomly." << std::endl;
				for (size_t i = 0; i < m_num_function; ++i) {
                    for (size_t j = 0; j < m_number_variables; ++j)
                        out << m_function[i]->translation()[j] << " ";
                    out << std::endl;
                }
				
				out.close();
			}
			else {
				std::string row;
				std::getline(in, row); // Skip the first line of commments
				for (size_t i = 0; i < m_num_function; ++i) {
					std::getline(in, row);
					std::stringstream sstr_row(row);
					for (size_t j = 0; j < m_number_variables; ++j)
						sstr_row >> m_function[i]->translation()[j];
				}
			}
			in.close();
			for(auto &i: m_function)
				i->setTranlated(true);
			return true;
		}

		bool Composition::loadRotation(const std::string &path) {
			std::string s;
			std::stringstream ss;
			ss << m_number_variables << "Dim.txt";
			s = ss.str();
			s.insert(0, "_F" + m_index_number + "_RotM_");

			s.insert(0, path);// data path
			s.insert(0, g_working_directory);

			for (auto &i : m_function)
				i->rotation().resize(m_number_variables, m_number_variables);
			std::ifstream in(s);
			if (in.fail()) {
				setRotation();
				std::ofstream out(s);
				out << "# The data is generated by OFEC randomly." << std::endl;
				for (size_t i = 0; i < m_num_function; ++i)
					m_function[i]->rotation().print(out);
				out.close();
			}
			else {
				std::string first_line;
				std::getline(in, first_line);
				for (size_t i = 0; i < m_num_function; ++i) 
					m_function[i]->rotation().load(in);	
			}
			in.close();
			for (auto &i : m_function)
				i->setRotated(true);
			return true;
		}

		void Composition::setTranslation() {
			for (int i = 0; i < m_num_function; i++) 
				for (int j = 0; j < m_number_variables; j++) 
					m_function[i]->translation()[j] = m_domain[j].limit.first + (m_domain[j].limit.second - m_domain[j].limit.first)*(1 - m_random->uniform.next());
		}

		void Composition::setRotation() {
			for (auto i : m_function) 
				i->rotation().generateRotationClassical(&m_random->normal, i->conditionNumber());	
		}

		Function* Composition::getFunction(size_t num) {
			return m_function[num].get();
		}
		
	}
}

