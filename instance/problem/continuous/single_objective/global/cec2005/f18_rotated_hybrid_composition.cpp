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

#include "f18_rotated_hybrid_composition.h"
#include "../classical/sphere.h"
#include "../classical/ackley.h"
#include "../classical/griewank.h"
#include "../classical/rastrigin.h"
#include "../classical/weierstrass.h"

namespace ofec::cec2005 {
	void RotatedHybridComposition2::addInputParameters() {
		m_index_number = "18";
	}

	bool RotatedHybridComposition2::loadTranslation(const std::string &path) {
		Composition::loadTranslation(path);
		for (size_t j = 0; j < m_number_variables; ++j)
			m_function[9]->translation()[j] = 0;
		return true;
	}

	void RotatedHybridComposition2::setFunction(Environment *env) {
		m_num_function = 10;
		m_function.resize(m_num_function);
		m_param_fun.resize(m_num_function);
		m_height.resize(m_num_function);
		m_fmax.resize(m_num_function);
		m_converge_severity.resize(m_num_function);
		m_stretch_severity.resize(m_num_function);

		BasicFunctions f(5);
		f[0] = &createFunction<Ackley>;
		f[1] = &createFunction<Rastrigin>;
		f[2] = &createFunction<Sphere>;
		f[3] = &createFunction<Weierstrass>;
		f[4] = &createFunction<Griewank>;

		for (size_t i = 0; i < m_num_function; ++i) {
			m_param_fun[i].reset(new ParameterMap(m_archived_parameters));
			m_function[i].reset(dynamic_cast<Function*>(f[i / 2]()));
			m_function[i]->reset();
			m_function[i]->inputParameters().input(*m_param_fun[i]);
			m_function[i]->initialize(env);
			m_function[i]->setBias(0);
		}

		m_function[0]->setConditionNumber(2); m_function[1]->setConditionNumber(3);
		m_function[2]->setConditionNumber(2); m_function[3]->setConditionNumber(3);
		m_function[4]->setConditionNumber(2); m_function[5]->setConditionNumber(3);
		m_function[6]->setConditionNumber(20); m_function[7]->setConditionNumber(30);
		m_function[8]->setConditionNumber(200); m_function[9]->setConditionNumber(300);

		for (int i = 0; i < m_num_function; i++) {
			m_height[i] = 100 * i;
		}

		m_function[0]->setDomain(-32, 32);   m_function[1]->setDomain(-32, 32);
		m_function[2]->setDomain(-5, 5);     m_function[3]->setDomain(-5, 5);
		m_function[4]->setDomain(-100, 100); m_function[5]->setDomain(-100, 100);
		m_function[6]->setDomain(-0.5, 0.5); m_function[7]->setDomain(-0.5, 0.5);
		m_function[8]->setDomain(-60, 60); m_function[9]->setDomain(-60, 60);

		m_stretch_severity[0] = 2.*5. / 32; m_stretch_severity[1] = 5. / 32;
		m_stretch_severity[2] = 2.;		m_stretch_severity[3] = 1.;
		m_stretch_severity[4] = 2 * 5. / 100;  m_stretch_severity[5] = 5. / 100;
		m_stretch_severity[6] = 20.;	m_stretch_severity[7] = 10.;
		m_stretch_severity[8] = 2. * 5 / 60;  m_stretch_severity[9] = 5. / 60;

		m_converge_severity[0] = 1.; m_converge_severity[1] = 2.;
		m_converge_severity[2] = 1.5;	m_converge_severity[3] = 1.5;
		m_converge_severity[4] = 1.;  m_converge_severity[5] = 1.;
		m_converge_severity[6] = 1.5;	m_converge_severity[7] = 1.5;
		m_converge_severity[8] = 2.;  m_converge_severity[9] = 2.;

		for (int i = 0; i < m_num_function; i++) {
			m_function[i]->setScale(m_stretch_severity[i]);
		}

		//setBias(10.);
	}

	void RotatedHybridComposition2::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Composition::evaluateObjective(x, obj);
		obj[0] += 10.;
	}

	void RotatedHybridComposition2::setTranslation() {
		for (size_t i = 0; i < m_num_function - 1; ++i)
			for (size_t j = 0; j < m_number_variables; ++j)
				m_function[i]->translation()[j] = m_domain[j].limit.first + (m_domain[j].limit.second - m_domain[j].limit.first) * (1 - m_random->uniform.next());
				
		for (size_t j = 0; j < m_number_variables; ++j) 
			m_function[m_num_function - 1]->translation()[j] = 0;
			
	}
	
	void RotatedHybridComposition2::updateOptima(Environment *env) {
		Composition::updateOptima(env);
		m_objective_accuracy = 1.0e-2;
	}
}
