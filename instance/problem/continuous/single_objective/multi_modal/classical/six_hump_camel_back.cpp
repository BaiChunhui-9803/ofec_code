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

#include "six_hump_camel_back.h"

namespace ofec {
	void SixHumpCamelBack::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
		m_input_parameters.add("objective accuracy", new EnumeratedReal(
			m_objective_accuracy, { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 }, 1e-4));
	}

	void SixHumpCamelBack::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(2);
		m_domain.setRange(-1.9, 1.9, 0);
		m_domain.setRange(-1.1, 1.1, 1);
	}

	void SixHumpCamelBack::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<std::vector<Real>> var_data = {
			{ 0.089842008935272, -0.712656403019058},
			{-0.089842008935272,  0.712656403019058} };
		for (auto &var_datum : var_data) {
			temp_sol.variable().vector() = var_datum;
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 0.5;
	}

	void SixHumpCamelBack::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s = x[0] * x[0], t = x[1] * x[1];
		s = (4 - 2.1*s + pow(x[0], 4) / 3)*s + x[0] * x[1] + (-4 + 4 * t)*t;
		obj[0] = s;
	}
	
}