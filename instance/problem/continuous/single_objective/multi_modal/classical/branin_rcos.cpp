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
*  Paper: Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
*******************************************************************************************/

#include "branin_rcos.h"

namespace ofec {
	void BraninRCOS::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(2);
		m_domain.setRange(-5, 10, 0);
		m_domain.setRange(0, 15, 1);	
	}

	void BraninRCOS::updateOptima(Environment *env) {
		m_variable_niche_radius = 0.1;
		m_objective_accuracy = 1.e-5;
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<std::vector<Real>> var_data = { { -(Real)OFEC_PI, 12.275f },{ (Real)OFEC_PI, 2.275f },{9.42478f, 2.475f} };
		for (auto &i : var_data) {
			temp_sol.variable()[0] = i[0];
			temp_sol.variable()[1] = i[1];
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
	}

	void BraninRCOS::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s, a;
		a = x[1] - 5.1*x[0] * x[0] / (4 * OFEC_PI*OFEC_PI) + 5 * x[0] / OFEC_PI - 6;
		s = a*a + 10 * (1 - 1 / (8 * OFEC_PI))*cos(x[0]) + 10;
		obj[0] = s;
	}
	
}