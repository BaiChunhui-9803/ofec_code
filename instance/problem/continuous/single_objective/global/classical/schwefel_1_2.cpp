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

#include "schwefel_1_2.h"

namespace ofec {
	void Schwefel_1_2::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void Schwefel_1_2::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-100, 100);
		setOriginalGlobalOpt();
	}

	void Schwefel_1_2::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void Schwefel_1_2::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real s1 = 0, s2 = 0;
		for (int i = 0; i < m_number_variables; i++) {
			for (int j = 0; j <= i; j++)
				s1 += x[j];
			s2 += s1*s1;
			s1 = 0;
		}
		obj[0] = s2 + m_bias;
	}
}