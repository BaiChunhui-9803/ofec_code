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

#include "weierstrass.h"

namespace ofec {
	void Weierstrass::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void Weierstrass::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-0.5, 0.5);
		m_a = 0.5;
		m_b = 3;
		m_kmax = 20;
		setOriginalGlobalOpt();
	}

	void Weierstrass::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void Weierstrass::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real fit = 0, s = 0;
		for (int i = 0; i < m_number_variables; i++)
			for (int k = 0; k <= m_kmax; k++)
				fit += pow(m_a, k) * cos(2 * OFEC_PI * pow(m_b, k) * (x[i] + 0.5));

		for (int k = 0; k <= m_kmax; k++)
			s += pow(m_a, k) * cos(2 * OFEC_PI * pow(m_b, k) * 0.5);
		s = s * m_number_variables;
		obj[0] = fit - s + m_bias;
	}

}