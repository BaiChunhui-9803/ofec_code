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
******************************************************************************************
*  Paper; A sequential niching memetic algorithm for Continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
/*******************************************************************************************
*  Paper: A sequential niching memetic algorithm for Continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
****************************************************************************************
*  LaTex:F(x)=\frac{1}{D}\sum^D_{i=1}{sin(10\log(x_i))}
*******************************************************************************************/

#include "vincent.h"

namespace ofec {
	void Vincent::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
		m_input_parameters.add("objective accuracy", new EnumeratedReal(
			m_objective_accuracy, { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 }, 1e-4));
	}

	void Vincent::initialize_(Environment *env) { // note
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		resizeVariable(m_number_variables);
		setDomain(0.25, 10.); // note
	}

	void Vincent::updateOptima(Environment *env) {
		if (m_number_variables > 5) {
			return;
		}
		m_optima.reset(new Optima<>());
		std::vector<size_t> divs = { 0 };
		std::vector<Real> opt_x(6);
		for (int i = 0; i < 6; i++)
			opt_x[i] = exp((4 * i - 7) * OFEC_PI / 20);
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		size_t cur_dim = 0;
		while (true) {
			if (divs[cur_dim] < 6)
				temp_sol.variable()[cur_dim] = opt_x[divs[cur_dim]];
			if (divs.back() == 6) {
				divs.pop_back();
				cur_dim--;
				if (divs.empty())
					break;
				divs.back()++;
			}
			else {
				if (divs.size() == m_number_variables) {
					evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
					dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
					divs.back()++;
				}
				else {
					divs.push_back(0);
					cur_dim++;
				}
			}
		}
		m_variable_niche_radius = 0.2;
	}

	void Vincent::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s = 0;
		for (int i = 0; i < m_number_variables; ++i) {
			s += sin(10 * log(x[i]));
		}
		s /= m_number_variables;
		obj[0] = s;  // note
	}
}