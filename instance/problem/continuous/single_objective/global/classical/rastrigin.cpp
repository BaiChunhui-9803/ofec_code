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

#include "rastrigin.h"

namespace ofec {
	void Rastrigin::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 100, 2));
	}

	void Rastrigin::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-5.12, 5.12);
		setOriginalGlobalOpt();
	}

	void Rastrigin::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void Rastrigin::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real fit = 0;
		for (int i = 0; i < m_number_variables; i++)
			fit += x[i] * x[i] - 10.*cos(2 * OFEC_PI*x[i]) + 10.;
		obj[0] = fit + m_bias;
	}
	
}