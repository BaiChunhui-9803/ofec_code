
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

#include "bent_cigar.h"

namespace ofec {
	void BentCigar::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void BentCigar::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-100., 100.);
		setOriginalGlobalOpt();
	}

	void BentCigar::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void BentCigar::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		obj[0] = x[0] * x[0];
		for (size_t i = 1; i< m_number_variables; ++i)	{
			obj[0] += pow(10.0, 6.0) * x[i] * x[i];
		}
		obj[0] += m_bias;
	}

}