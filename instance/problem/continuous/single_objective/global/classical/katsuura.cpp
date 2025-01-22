/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

#include "katsuura.h"

namespace ofec {
	void katsuura::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void katsuura::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-100., 100.);
		setOriginalGlobalOpt();
	}

	void katsuura::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void katsuura::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		if (m_translated)
			translate(x);
		if (m_scaled)
			scale(x);
		if (m_rotated)
			rotate(x);
		if (m_translated)
			translateOrigin(x);
		size_t i,j;
		Real temp, tmp1, tmp2, tmp3;
		obj[0] = 1.0;
		tmp3 = pow(1.0*m_number_variables, 1.2);
		for (i = 0; i<m_number_variables; i++) {
			temp = 0.0;
			for (j = 1; j <= 32; ++j) {
				tmp1 = pow(2.0, j);
				tmp2 = tmp1*x[i];
				temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
			}
			obj[0] *= pow(1.0 + (i + 1)*temp, 10.0 / tmp3);
		}
		tmp1 = 10.0 / m_number_variables / m_number_variables;
		obj[0] = obj[0] * tmp1 - tmp1;
		obj[0] += m_bias;
	}
}