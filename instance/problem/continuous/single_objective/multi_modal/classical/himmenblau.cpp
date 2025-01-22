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
*  Paper: A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/

#include "himmenblau.h"

namespace ofec {
	void Himmenblau::addInputParameters() {
		m_input_parameters.add("objective accuracy", new EnumeratedReal(
			m_objective_accuracy, { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 }, 1e-4));
	}

	void Himmenblau::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(2);
		setDomain(-6, 6);
	}

	void Himmenblau::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<std::vector<Real>> var_data = {
			{ 3.0,2.0},
			{-2.805118094822989,  3.131312538494919},
			{-3.779310265963066, -3.283185984612214},
			{ 3.584428351760445, -1.848126540197251} };
		for (auto &var_datum : var_data) {
			temp_sol.variable().vector() = var_datum;
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
		m_variable_niche_radius = 0.01;
	}

	void Himmenblau::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s = 0;
		Real t0 = (x[0] * x[0] + x[1] - 11), t1 = (x[1] * x[1] + x[0] - 7);
		s = t0 * t0 + t1 * t1;
		obj[0] = s;
	}

}