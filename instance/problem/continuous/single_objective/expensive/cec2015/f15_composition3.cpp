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

#include "F15_expensive_composition3.h"
#include "../../global/classical/HGBat.h"
#include "../../global/classical/rastrigin.h"
#include "../../global/classical/schwefel.h"
#include "../../global/classical/weierstrass.h"
#include "../../global/classical/elliptic.h"



namespace ofec {
	namespace cec2015 {
		F15_expensive_composition3::F15_expensive_composition3(const ParameterMap &v) :
			F15_expensive_composition3((v.at("problem name")), (v.at("number of variables")), 1) {
			
		}
		F15_expensive_composition3::F15_expensive_composition3(const std::string &name, size_t size_var, size_t size_obj) :problem(name, size_var, size_obj), \
			composition_2015(name, size_var, size_obj) {
			
		}
		
		void F15_expensive_composition3::setFunction() {
			BasicFunctions f(5);
			f[0] = &createFunction<HGBat>;
			f[1] = &createFunction<rastrigin>;
			f[2] = &createFunction<schwefel>;
			f[3] = &createFunction<weierstrass>;
			f[4] = &createFunction<elliptic>;

			for (size_t i = 0; i < m_num_function; ++i) {
				m_function[i] = dynamic_cast<function*>(f[i]("", m_number_variables, m_number_objectives));
				m_function[i]->initialize();
				m_function[i]->setBias(0);
			}

			for (auto &i : m_function)
				i->setConditionNumber(2.);

			m_converge_severity[0] = 10;
			m_converge_severity[1] = 10;
			m_converge_severity[2] = 10;
			m_converge_severity[3] = 20;
			m_converge_severity[4] = 20;

			m_height[0] = 10;
			m_height[1] = 10;
			m_height[2] = 2.5;
			m_height[3] = 25;
			m_height[4] = 1e-6;

			m_f_bias[0] = 0;
			m_f_bias[1] = 100;
			m_f_bias[2] = 200;
			m_f_bias[3] = 300;
			m_f_bias[4] = 400;

			//setBias(1500);
		}
		void F15_expensive_composition3::initialize() {
			m_num_function = 5;
			m_function.resize(m_num_function);
			m_height.resize(m_num_function);
			m_converge_severity.resize(m_num_function);
			m_f_bias.resize(m_num_function);
			setDomain(-100., 100.);
			setInitialDomain(-100., 100.);
			m_variable_monitor = true;

			setFunction();
			loadTranslation("/instance/problem/continuous/expensive/cec2015");  //data path
			
			loadRotation("/instance/problem/continuous/expensive/cec2015");  //data path
			
			for (auto &i : m_function) {
				i->get_optima().clear();
				i->setGlobalOpt(i->translation().data());
			}
			// Set optimal solution
            Solution<VariableVector<Real>, Real> s(m_number_objectives, num_constraints(), m_number_variables);
			for (int i = 0; i < m_number_variables; ++i) {
				s.variable()[i] = m_function[0]->get_optima().variable(0)[i];
			}
			m_optima->append(s.variable());

            s.evaluate(false, caller::Problem);
            m_optima->append(s.objective());
			m_initialized = true;
		}

		int F15_expensive_composition3::evaluateObjective(Real *x, std::vector<Real> &obj) {
			composition_2015::evaluateObjective(x, obj);
			obj[0] += 1500;
			return kNormalEval;
		}

	}
}


