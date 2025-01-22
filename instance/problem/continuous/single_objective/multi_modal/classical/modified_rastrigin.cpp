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
/*******************************************************************************************
*  Paper: A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
****************************************************************************************/

#include "modified_rastrigin.h"

namespace ofec {
	void ModifiedRastrigin::addInputParameters() {
		m_input_parameters.add("objective accuracy", new EnumeratedReal(
			m_objective_accuracy, { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 }, 1e-4));
	}

	void ModifiedRastrigin::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(2);
		setDomain(0, 1); 
		m_k.resize(2);
		m_k[0] = 3; m_k[1] = 4;
	}

	void ModifiedRastrigin::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>());
		Real var_data[12][2] = {
			1 / 6.0, 1 / 8.0,
			1 / 6.0, 3 / 8.0,
			1 / 6.0, 5 / 8.0,
			1 / 6.0, 7 / 8.0,
			3 / 6.0, 1 / 8.0,
			3 / 6.0, 3 / 8.0,
			3 / 6.0, 5 / 8.0,
			3 / 6.0, 7 / 8.0,
			5 / 6.0, 1 / 8.0,
			5 / 6.0, 3 / 8.0,
			5 / 6.0, 5 / 8.0,
			5 / 6.0, 7 / 8.0 };
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		for (int i = 0; i < 12; ++i) {
			temp_sol.variable()[0] = var_data[i][0];
			temp_sol.variable()[1] = var_data[i][1];
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 1.e-2;
	}

	void ModifiedRastrigin::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s = 0;
		for (int i = 0; i < m_number_variables; ++i) {
			s += 10 + 9 * cos(2 * OFEC_PI * m_k[i] * x[i]);
		}
		obj[0] = s;  // note
	}
	
}