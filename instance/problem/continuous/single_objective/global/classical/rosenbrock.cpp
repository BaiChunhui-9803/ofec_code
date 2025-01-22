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

#include "rosenbrock.h"

namespace ofec {
	void Rosenbrock::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void Rosenbrock::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-2.048, 2.048);
		std::vector<Real> v(m_number_variables, 1);
		setOriginalGlobalOpt(v.data());
	}

	void Rosenbrock::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void Rosenbrock::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real fitness = 0;
		for (int i = 0; i < m_number_variables - 1; i++) {
			fitness += 100 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]) + (x[i] - 1) * (x[i] - 1);
		}
		obj[0] = fitness + m_bias;
	}

}