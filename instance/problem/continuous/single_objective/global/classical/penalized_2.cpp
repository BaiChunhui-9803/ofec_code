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

#include "penalized_2.h"

namespace ofec {
	void Penalized_2::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void Penalized_2::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-50, 50);
		std::vector<Real> var(m_number_variables, 1);
		setOriginalGlobalOpt(var.data());
	}

	void Penalized_2::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void Penalized_2::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real s = 0;
		for (int i = 0; i < m_number_variables - 1; i++)
			s += (x[i] - 1)*(x[i] - 1)*(1 + sin(3 * OFEC_PI*x[i + 1])*sin(3 * OFEC_PI*x[i + 1]));
		s += (x[m_number_variables - 1] - 1)*(x[m_number_variables - 1] - 1)*(1 + sin(2 * OFEC_PI*x[m_number_variables - 1])*sin(2 * OFEC_PI*x[m_number_variables - 1])) + sin(3 * OFEC_PI*x[0])*sin(3 * OFEC_PI*x[0]);
		s = s*0.1;
		for (int i = 0; i < m_number_variables; i++)
			s += u(x[i], 5, 100, 4);
		obj[0] = s + m_bias;
	}

	Real Penalized_2::u(Real x, Real a, Real k, Real m)const {
		if (x > a) return k*pow(x - a, m);
		else if (x < -a) return k*pow(-x - a, m);
		else return 0;
	}
}
