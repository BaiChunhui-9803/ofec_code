/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (ofec)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com 
* Language: C++
*************************************************************************
*  This file is part of ofec. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper: A sequential niching memetic algorithm for continuous multimodal
*		  Appled Mathematics and Computation 218(2012) 8242-8259
*******************************************************************************************/
// Created: 21 July 2011
// Last modified:

#include "fiba.h"

namespace ofec {
	void FIBA::addInputParameters() {
		m_input_parameters.add("case", new EnumeratedInt(m_case, { 1, 2 }, 1));
	}

	void FIBA::initialize_(Environment *env) {
		Function::initialize_(env);
		resizeVariable(2);
		setDomain(-4.0, 4.0);
		setCase(m_case);
	}

	void FIBA::updateOptima(Environment *env) {
		m_objective_accuracy = 1.e-6;
		m_optima.reset(new Optima<>(m_original_optima));
	}

	void FIBA::evaluateOriginalObj(Real *x, std::vector<Real> &obj) const {

		Real s = 0;
		Real t0 = x[0] * x[0] + x[1] * x[1];
		s = (t0) / (1 + t0) + m_kappa*(14 * t0 + pow(t0, 2)*m_chi*m_chi - 2 * sqrt(14.)*(pow(x[0], 3) - 3 * x[0] * x[1] * x[1])*m_chi) / (14 * (pow(1 + t0, 2)));
		obj[0] = s + m_bias;
	}

	void FIBA::setCase(int c) {
		m_case = c;
		if (m_case == 1) {
			m_variable_niche_radius = 0.08;
			std::vector<std::vector<Real>> var_data = {
				{ 0.45186f, 0.0f },
				{ -0.22593f, 0.39132f },
				{ -0.22593f, -0.39132f },
				{ 0.0f, 0.0f }
			};
			for (auto &i : var_data) {
				setOriginalGlobalOpt(i.data());
			}
		}
		else {
			m_variable_niche_radius = 0.5;
			std::vector<std::vector<Real>> var_data = {
				{ 0.0f, 0.0f },
				{ 1.2243f, 0.0f },
				{ -0.61215f, 1.0603f },
				{ -0.61215f, -1.0603f }
			};
			for (auto &i : var_data) {
				setOriginalGlobalOpt(i.data());
			}
		}
		setPara();
	}

	void FIBA::setPara() {
		if (m_case == 1) {
			m_kappa = -0.95;
			m_chi = -1.26;
		}
		else {
			m_kappa = 4;
			m_chi = 2;
		}
	}
}