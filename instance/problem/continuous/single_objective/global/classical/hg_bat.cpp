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

#include "hg_bat.h"

namespace ofec {
	void HGBat::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void HGBat::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-50, 50);
		setOriginalGlobalOpt();
	}

	void HGBat::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void HGBat::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real alpha, r2, sum_x;
		alpha = 1.0 / 4.0;
		r2 = 0.0;
		sum_x = 0.0;
		for (size_t i = 0; i<m_number_variables; ++i) {
			x[i] = x[i] - 1.0;//shift to orgin
			r2 += x[i] * x[i];
			sum_x += x[i];
		}
		obj[0] = pow(fabs(pow(r2, 2.0) - pow(sum_x, 2.0)), 2 * alpha) + (0.5*r2 + sum_x) / m_number_variables + 0.5;
		obj[0] += m_bias;
	}

}