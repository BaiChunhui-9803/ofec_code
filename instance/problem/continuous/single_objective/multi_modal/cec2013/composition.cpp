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

#include "composition.h"

namespace ofec::cec2013 {
	void Composition::initialize_(Environment *env) {
		Continuous::initialize_(env);
		m_number_objectives = 1;
		m_optimize_mode.resize(m_number_objectives);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		m_domain.resize(m_number_variables);
		setDomain(-5., 5.);
		m_height_normalize_severity = 2000;
		setFunction(env);
		loadRotation("/instance/problem/continuous/single_objective/multi_modal/cec2013/MMOP_CEC2013");
		computeFmax();
		loadTranslation("/instance/problem/continuous/single_objective/multi_modal/cec2013/MMOP_CEC2013");  //data path
	}

	void Composition::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		for (auto &i : m_function) {
			temp_sol.variable().vector() = i->translation();
			evaluate(temp_sol.variableBase(), temp_sol.objective(), temp_sol.constraint());
			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		}
	}

	void Composition::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		cec2005::Composition::evaluateObjective(x, obj);
		obj[0] = -obj[0];
	}
}


