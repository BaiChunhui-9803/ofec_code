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

#include "penalized_1.h"

namespace ofec {
	void Penalized_1::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void Penalized_1::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-50, 50);
		std::vector<Real> var(m_number_variables, -1);
		setOriginalGlobalOpt(var.data());
	}

	void Penalized_1::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void Penalized_1::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		std::vector<Real> y(m_number_variables);
		for (int i = 0; i < m_number_variables; i++) y[i] = (x[i] + 1) / 4. + 1;
		Real s = 0;
		for (int i = 0; i < m_number_variables - 1; i++)
			s += (y[i] - 1)*(y[i] - 1)*(1 + 10 * sin(OFEC_PI*y[i + 1])*sin(OFEC_PI*y[i + 1]));
		s += (y[m_number_variables - 1] - 1)*(y[m_number_variables - 1] - 1) + 10 * sin(OFEC_PI*y[0])*sin(OFEC_PI*y[0]);
		s = s*OFEC_PI / m_number_variables;
		for (int i = 0; i < m_number_variables; i++) {
			s += u(x[i], 10, 100, 4);
		}
		obj[0] = s + m_bias;
	}

	Real Penalized_1::u(Real x, Real a, Real k, Real m)const {
		if (x > a) return k*pow(x - a, m);
		else if (x < -a) return k*pow(-x - a, m);
		else return 0;
	}
}