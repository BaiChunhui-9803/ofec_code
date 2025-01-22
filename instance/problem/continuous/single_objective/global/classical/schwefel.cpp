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

#include "schwefel.h"

namespace ofec {
	void Schwefel::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void Schwefel::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-500, 500);
		std::vector<Real> var(m_number_variables, 420.9687);
		setOriginalGlobalOpt(var.data());
	}

	void Schwefel::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1.0e-2;
	}

	void Schwefel::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real fitness = 0;
		for (int i = 0; i < m_number_variables; i++) {
			fitness += -x[i] * sin(sqrt(fabs(x[i])));
		}
		obj[0] = fitness + m_bias;
	}	
}