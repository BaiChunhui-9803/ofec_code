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

#include "non_continuous_rastrigin.h"

namespace ofec {
	void NonContinuousRastrigin::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void NonContinuousRastrigin::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-5.12, 5.12);
		setOriginalGlobalOpt();
	}

	void NonContinuousRastrigin::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void NonContinuousRastrigin::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real fit = 0;
		std::vector<Real> y(m_number_variables);
		for (size_t i = 0; i < m_number_variables; ++i) {
			if (fabs(x[i]) < 0.5) y[i] = x[i];
			else {
				Real a, b = 2 * x[i];
				if (b <= 0 && b - (int)b < 0.5) a = (int)b;
				else if (b - (int)b < 0.5) a = (int)b;
				else a = (int)b + 1;
				y[i] = a / 2;
			}
		}
		for (size_t i = 0; i < m_number_variables; ++i) {
			fit = fit + y[i] * y[i] - 10.*cos(2 * OFEC_PI*y[i]) + 10.;
		}
		obj[0] = fit + m_bias;
	}
	
}
