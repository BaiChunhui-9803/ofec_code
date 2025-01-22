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

#include "max_global.h"

namespace ofec {
	void MaxGlobal::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		resizeVariable(1);
		setDomain(0, 1);
	}

	void MaxGlobal::updateOptima(Environment *env) {
		m_variable_niche_radius = 0.1;
		m_objective_accuracy = 1.e-5;
		//5 gopt
		m_optima.reset(new Optima<>());
		std::vector<Real> opt_x = { 0.5f, 0.1f, 0.3f, 0.7f, 0.9f };
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		for (Real x : opt_x) {
			temp_sol.variable()[0] = x;
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
	}

	void MaxGlobal::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		obj[0] = pow(sin(5 * OFEC_PI*x[0]), 6.);
	}
	
}