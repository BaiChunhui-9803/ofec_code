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

#include "f10_composition2.h"
#include "../../global/classical/griewank.h"
#include "../../global/classical/weierstrass.h"
#include "../../global/classical/sphere.h"
#include "../../global/classical/rastrigin.h"

namespace ofec::cec2013 {
	void Composition2::addInputParameters() {
		m_input_parameters.add("objective accuracy", new EnumeratedReal(
			m_objective_accuracy, { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 }, 1e-4));
		m_index_number = "10";
	}

	void Composition2::updateOptima(Environment *env) {
		Composition::updateOptima(env);
		m_variable_niche_radius = 0.01;
	}

	void Composition2::setFunction(Environment *env) {
		m_num_function = 8;
		m_function.resize(m_num_function);
		m_param_fun.resize(m_num_function);
		m_height.resize(m_num_function);
		m_fmax.resize(m_num_function);
		m_converge_severity.resize(m_num_function);
		m_stretch_severity.resize(m_num_function);

		BasicFunctions f(4);
		f[0] = &createFunction<Rastrigin>;
		f[1] = &createFunction<Weierstrass>;
		f[2] = &createFunction<Griewank>;
		f[3] = &createFunction<Sphere>;

		for (size_t i = 0; i < m_num_function; ++i) {
			m_param_fun[i].reset(new ParameterMap(m_archived_parameters));
			m_function[i].reset(dynamic_cast<Function*>(f[i / 2]()));
			m_function[i]->reset();
			m_function[i]->inputParameters().input(*m_param_fun[i]);
			m_function[i]->initialize(env);
			m_function[i]->setBias(0);
		}

		m_stretch_severity[0] = 1;				m_stretch_severity[1] = 1;
		m_stretch_severity[2] = 10;				 m_stretch_severity[3] = 10;
		m_stretch_severity[4] = 1. / 10;			 m_stretch_severity[5] = 1. / 10;
		m_stretch_severity[6] = 1 / 7.;			m_stretch_severity[7] = 1 / 7.;

		for (int i = 0; i<m_num_function; i++) {
			m_function[i]->setScale(m_stretch_severity[i]);
		}

		for (int i = 0; i<m_num_function; i++) {
			m_converge_severity[i] = 1;
		}

		for (int i = 0; i<m_num_function; i++) {
			m_height[i] = 0;
		}
	}

	void Composition2::setRotation() {
		for (auto i : m_function)
			i->rotation().identify();
	}
}


