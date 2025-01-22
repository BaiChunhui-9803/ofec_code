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
****************************************************************************************
*  LaTex:F_4(x)=\mathrm{e}^{-2\lg2*{{\frac{x-0.08}{0.854}}^2}}*{\sin{5\pi*{x^0.75-0.05}}}^6
*******************************************************************************************/

#include "beasley.h"

namespace ofec {
	void Beasley::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		resizeVariable(1);
		setDomain(0, 1.);
	}

	void Beasley::updateOptima(Environment *env) {
		m_variable_niche_radius = 0.01;
		m_objective_accuracy = 1.e-6;
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		std::vector<Real> var_data = { 0.08f, 0.25f, 0.45f, 0.68f, 0.93f };
		for (auto &i : var_data) {
			temp_sol.variable()[0] = i;
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
	}

	void Beasley::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s;
		s = exp(-2 * log(2.) * ((x[0] - 0.08) / 0.854) * ((x[0] - 0.08) / 0.854)) * pow(sin(5 * OFEC_PI * (pow(x[0], 0.75) - 0.05)), 6);
		obj[0] = s;
	}
}