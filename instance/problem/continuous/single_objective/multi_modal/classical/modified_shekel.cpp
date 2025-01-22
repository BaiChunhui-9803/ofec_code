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
*  Paper; A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/

#include "modified_shekel.h"

namespace ofec {
	void ModifiedShekel::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 5, 2));
	}

	void ModifiedShekel::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		resizeVariable(m_number_variables);
		setDomain(0.0, 11.0);
	}

	void ModifiedShekel::updateOptima(Environment *env) {
		 //1 gopt+7 lopt
		m_optima.reset(new Optima<>());
		Real a[8][5] = {
			4,4,6.3,4,4,1,1,8.5,
			1,1,6,6,9.1,6,6,
			3.5,7.5,4,9,4,
			5,5,3,3,9,
			9.1,8.2,2,3,9,
			1.5,9.3,7.4,3,9,
			7.8,2.2,5.3,9,3 };
		Real c[8] = { 0.1,0.2,0.4,0.15,0.6,0.2,0.06,0.18 };
		std::copy(c, c + 8, m_c);
		for (size_t i = 0; i < 8; ++i)
			std::copy(a[i], a[i] + 5, m_a[i]);
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		for (size_t j = 0; j < m_number_variables; j++) {
			temp_sol.variable()[j] = m_a[6][j];
		}
		evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
		dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		for (size_t i = 0; i < 8; ++i) {
			if (i == 6) continue;
			for (size_t j = 0; j < m_number_variables; j++) {
				temp_sol.variable()[j] = m_a[i][j];
			}
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 0.2;
		m_objective_accuracy = 1.e-3;
	}

	void ModifiedShekel::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s = 0;
		for (int i = 0; i < 8; ++i) {
			Real b = 0;
			for (int j = 0; j < m_number_variables; ++j) {
				b += (x[j] - m_a[i][j])*(x[j] - m_a[i][j]);
			}
			s += pow(b + m_c[i], -1);
		}
		obj[0] = s;
	}
	
}