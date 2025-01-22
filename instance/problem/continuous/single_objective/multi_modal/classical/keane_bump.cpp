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
*  Paper; Multimodal Optimization by Means of a Topological Species Conservation Algorithm
*		  IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, VOL.14,NO.6,DECEMBER 2010
*******************************************************************************************/

#include "keane_bump.h"

namespace ofec {
	void KeaneBump::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void KeaneBump::initialize_(Environment *env) { 
		Continuous::initialize_(env);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMaximize;
		resizeVariable(m_number_variables);
		setDomain(1e-8, 10);
	}

	void KeaneBump::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Real s, a = 0, b = 1, c = 0;
		int i;
		for (i = 0; i < m_number_variables; i++) {
			a += pow(cos(x[i]), 4);
			b *= cos(x[i]) * cos(x[i]);
			c += (i + 1) * x[i] * x[i];
		}
		s = abs((a - 2 * b) / sqrt(c));
		obj[0] = s;
	}
	
}