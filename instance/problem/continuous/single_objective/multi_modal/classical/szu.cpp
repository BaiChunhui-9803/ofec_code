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
******************************************************************************************
*  Paper; A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
* //* 
//* F(vec3{X})=\sum_{i=1}^{D}{-x_i^2}
//*
*******************************************************************************************/

#include "szu.h"

namespace ofec {
	void Szu::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 2, 9, 2));
	}

	void Szu::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		resizeVariable(m_number_variables);
		setDomain(-5.0, 5.0);
	}

	void Szu::updateOptima(Environment *env) {
		m_optima.reset(new Optima<>());
		switch (m_number_variables) {
		case 2:	m_optima->appendObjective(-156.66); break;
		case 3:	m_optima->appendObjective(-235.00); break;
		case 4:	m_optima->appendObjective(-313.33); break;
		case 5:	m_optima->appendObjective(-391.66); break;
		case 6:	m_optima->appendObjective(-469.99); break;
		case 7:	m_optima->appendObjective(-548.33); break;
		case 8:	m_optima->appendObjective(-626.66); break;
		case 9:	m_optima->appendObjective(-704.99); break;
		}
		m_objective_accuracy = 0.1;
	}

	void Szu::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s = 0;
		size_t i;
		for (i = 0; i < m_number_variables; ++i) {
			s += pow(x[i], 4) - 16 * x[i] * x[i] + 5 * x[i];
		}
		obj[0] = s;
	}
	
}