
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

#include "ackley.h"

namespace ofec {
	void Ackley::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void Ackley::initialize_(Environment *env) {
		Function::initialize_(env);
		resizeVariable(m_number_variables);
		setDomain(-32.768, 32.768);
		m_optimize_mode[0] = OptimizeMode::kMinimize;	
		setOriginalGlobalOpt();
	}

	void Ackley::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void Ackley::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real s1 = 0, s2 = 0;
		for (int i = 0; i < m_number_variables; i++) {
			s1 += x[i] * x[i];
			s2 += cos(2 * OFEC_PI*x[i]);
		}
		obj[0] = -20 * exp(-0.2*sqrt(s1 / m_number_variables)) - exp(s2 / m_number_variables) + 20 + OFEC_E + m_bias;
	}
	
}