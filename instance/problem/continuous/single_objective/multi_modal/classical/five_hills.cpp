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

#include "five_hills.h"

namespace ofec {
	void FiveHills::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;

		resizeVariable(2);
		m_domain.setRange(-2.5, 3., 0);
		m_domain.setRange(-2, 2., 1);
	}

	void FiveHills::updateOptima(Environment *env) {
		 //1 gopt + 4 lopt
		m_optima.reset(new Optima<>());
		Real var_data[5][3] = {
			-8.24371e-014,-1.53082e-013,2.5,
			0.889286,-1.08335e-016,1.60084,
			-0.889286,-8.32809e-017,1.60084,
			-1.78391,1.09381e-016,0.699112,
			1.78391,1.09381e-016,0.699112 };
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		for (int i = 0; i < 5; ++i) {
			temp_sol.variable()[0] = var_data[i][0];
			temp_sol.variable()[1] = var_data[i][1];
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 0.2;
		m_objective_accuracy = 1.e-5;
	}

	void FiveHills::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s;
		s = sin(2.2*OFEC_PI*x[0] + 0.5*OFEC_PI)*(2 - fabs(x[1])) / 2 * (3 - fabs(x[0])) / 2 + sin(0.5*OFEC_PI*x[1] + 0.5*OFEC_PI)*(2 - fabs(x[1])) / 2 * (2 - fabs(x[0])) / 2;
		obj[0] = s;
	}
	
}