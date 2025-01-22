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

#include "f11_composition3.h"
#include "../../global/classical/griewank.h"
#include "../../global/classical/weierstrass.h"
#include "../../global/classical/griewank_rosenbrock.h"

namespace ofec::cec2013 {
	void Composition3::addInputParameters() {
		m_input_parameters.add("objective accuracy", new EnumeratedReal(
			m_objective_accuracy, { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 }, 1e-4));
		m_index_number = "11";
	}

	void Composition3::updateOptima(Environment *env) {
		Composition::updateOptima(env);
		m_variable_niche_radius = 0.01;
	}

	void Composition3::setFunction(Environment *env) {
		m_num_function = 6;
		m_function.resize(m_num_function);
		m_param_fun.resize(m_num_function);
		m_fmax.resize(m_num_function);
		m_stretch_severity.resize(m_num_function);
		m_converge_severity.resize(m_num_function);
		m_height.resize(m_num_function);

		BasicFunctions f(3);
		f[0] = &createFunction<GriewankRosenbrock>;
		f[1] = &createFunction<Weierstrass>;
		f[2] = &createFunction<Griewank>;

		for (size_t i = 0; i < m_num_function; ++i) {
			m_param_fun[i].reset(new ParameterMap(m_archived_parameters));
			m_function[i].reset(dynamic_cast<Function*>(f[i / 2]()));
			m_function[i]->reset();
			m_function[i]->inputParameters().input(*m_param_fun[i]);
			m_function[i]->initialize(env);
			m_function[i]->setBias(0);
		}

		m_stretch_severity[0] = 1 / 4.; m_stretch_severity[1] = 1 / 10.;
		m_stretch_severity[2] = 2; m_stretch_severity[3] = 1;
		m_stretch_severity[4] = 2; m_stretch_severity[5] = 5;

		for (size_t i = 0; i<m_num_function; ++i) {
			m_function[i]->setScale(m_stretch_severity[i]);
		}

		for (size_t i = 0; i<m_num_function; ++i) {
			if (i == 0 || i == 1)	m_converge_severity[i] = 1;
			else m_converge_severity[i] = 2;
		}

		for (size_t i = 0; i<m_num_function; ++i) {
			m_height[i] = 0;
			m_function[i]->setConditionNumber(1.);
		}
	}
}



