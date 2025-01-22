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

#include "griewank.h"

namespace ofec {
	void Griewank::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void Griewank::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-600, 600);
		setOriginalGlobalOpt();
	}

	void Griewank::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void Griewank::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real result = 0;
		Real s1 = 0, s2 = 1;
		for (int i = 0; i < m_number_variables; i++) {
			s1 += x[i] * x[i] / 4000.;
			s2 *= cos(x[i] / sqrt((Real)(i + 1)));
		}
		result = s1 - s2 + 1. + m_bias;
		obj[0] = result;
	}
	
}