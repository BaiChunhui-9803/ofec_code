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

#include "modified_schwefel.h"

namespace ofec {
	void ModifiedSchwefel::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void ModifiedSchwefel::initialize_(Environment *env) {
		Function::initialize_(env);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-100., 100.);
		setOriginalGlobalOpt();
	}

	void ModifiedSchwefel::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>(m_original_optima));
		m_objective_accuracy = 1e-8;
	}

	void ModifiedSchwefel::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {
		Real tmp;
		obj[0] = 0.0;
		for (size_t i = 0; i < m_number_variables; i++){
			x[i] += 4.209687462275036e+002;
			if (x[i]>500) {
				obj[0] -= (500.0 - fmod(x[i], 500))*sin(pow(500.0 - fmod(x[i], 500), 0.5));
				tmp = (x[i] + 500.0) / 100;
				obj[0] += tmp*tmp / m_number_variables;
			}
			else if (x[i]<-500)	{
				obj[0] -= (-500.0 + fmod(fabs(x[i]), 500))*sin(pow(500.0 - fmod(fabs(x[i]), 500), 0.5));
				tmp = (x[i] + 500.0) / 100;
				obj[0] += tmp*tmp / m_number_variables;
			}
			else
				obj[0] -= x[i] * sin(pow(fabs(x[i]), 0.5));
		}
		obj[0] += 4.189828872724338e+002*m_number_variables;
		obj[0] += m_bias;
	}
}