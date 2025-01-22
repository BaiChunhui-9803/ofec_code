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

#include "f25_rotated_hybrid_no_bound.h"
#include "../classical/griewank_rosenbrock.h"
#include "../classical/ackley.h"
#include "../classical/rastrigin.h"
#include "../classical/weierstrass.h"
#include "../classical/griewank.h"
#include "../classical/non_continuous_scaffer_F6.h"
#include "../classical/non_continuous_rastrigin.h"
#include "../classical/rotated_scaffer_F6.h"
#include "../classical/elliptic.h"
#include "../classical/sphere_noisy.h"

namespace ofec::cec2005 {
	void RotatedHybridNoBound::addInputParameters() {
		m_index_number = "25";
	}

	void RotatedHybridNoBound::initialize_(Environment *env) {
		InitPopBounded::initialize_(env);
		Composition::initialize_(env);
		for (size_t i = 0; i < m_number_variables; i++)
			m_domain[i].limited = false;
		setDomainInitPop(2, 5);
	}

	void RotatedHybridNoBound::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Composition::evaluateObjective(x, obj);
		obj[0] += 260.;
	}

	void RotatedHybridNoBound::setTranslation() {	
		for (size_t j = 0; j < m_number_variables; ++j) 
			m_function[0]->translation()[j] = -5.0 + (7.0)*(m_random->uniform.next());

		for (size_t i = 1; i < m_num_function - 1; ++i)
			for (int j = 0; j < m_number_variables; ++j) 
				m_function[i]->translation()[j] = m_domain[j].limit.first + (m_domain[j].limit.second - m_domain[j].limit.first)*(1 - m_random->uniform.next());
		
		for (size_t j = 0; j < m_number_variables; ++j) 
			m_function[m_num_function - 1]->translation()[j] = 0;
	
	}

	void RotatedHybridNoBound::setFunction(Environment *env) {
		m_num_function = 10;
		m_function.resize(m_num_function);
		m_param_fun.resize(m_num_function);
		m_height.resize(m_num_function);
		m_fmax.resize(m_num_function);
		m_converge_severity.resize(m_num_function);
		m_stretch_severity.resize(m_num_function);

		BasicFunctions f(10);
		f[0] = &createFunction<Weierstrass>;
		f[1] = &createFunction<RotatedScafferF6>;
		f[2] = &createFunction<GriewankRosenbrock>;
		f[3] = &createFunction<Ackley>;
		f[4] = &createFunction<Rastrigin>;
		f[5] = &createFunction<Griewank>;
		f[6] = &createFunction<NonContinuousScafferF6>;
		f[7] = &createFunction<NonContinuousRastrigin>;
		f[8] = &createFunction<Elliptic>;
		f[9] = &createFunction<SphereNoisy>;

		for (size_t i = 0; i < m_num_function; ++i) {
			m_param_fun[i].reset(new ParameterMap(m_archived_parameters));
			m_function[i].reset(dynamic_cast<Function*>(f[i]()));
			m_function[i]->reset();
			m_function[i]->inputParameters().input(*m_param_fun[i]);
			m_function[i]->initialize(env);
			m_function[i]->setBias(0);
		}
		m_function[1]->loadRotation("/instance/problem/continuous/single_objective/global/classical/GOP_CEC2005_F25_part2");
		m_function[9]->setRandom(m_random);

		m_function[0]->setConditionNumber(100); m_function[1]->setConditionNumber(50);
		m_function[2]->setConditionNumber(30); m_function[3]->setConditionNumber(10);
		m_function[4]->setConditionNumber(5); m_function[5]->setConditionNumber(5);
		m_function[6]->setConditionNumber(4); m_function[7]->setConditionNumber(3);
		m_function[8]->setConditionNumber(2); m_function[9]->setConditionNumber(2);

		for (int i = 0; i < m_num_function; i++) {
			m_height[i] = 100 * i;
		}

		m_function[0]->setDomain(-0.5, 0.5);
		m_function[1]->setDomain(-100, 100);
		m_function[2]->setDomain(-5, 5);
		m_function[3]->setDomain(-32, 32);
		m_function[4]->setDomain(-5, 5);
		m_function[5]->setDomain(-5, 5);
		m_function[6]->setDomain(-100, 100);
		m_function[7]->setDomain(-5, 5);
		m_function[8]->setDomain(-100, 100);
		m_function[9]->setDomain(-100, 100);

		m_stretch_severity[0] = 10;
		m_stretch_severity[1] = 5.0 / 20;
		m_stretch_severity[2] = 1.0;
		m_stretch_severity[3] = 5.0 / 32;
		m_stretch_severity[4] = 1.0;
		m_stretch_severity[5] = 5.0 / 100;
		m_stretch_severity[6] = 5.0 / 50;
		m_stretch_severity[7] = 1.0;
		m_stretch_severity[8] = 5.0 / 100;
		m_stretch_severity[9] = 5.0 / 100;

		m_converge_severity[0] = 2.; ; m_converge_severity[1] = 2.;
		m_converge_severity[2] = 2.;	m_converge_severity[3] = 2.;
		m_converge_severity[4] = 2.;  m_converge_severity[5] = 2.;
		m_converge_severity[6] = 2;	m_converge_severity[7] = 2;
		m_converge_severity[8] = 2.;  m_converge_severity[9] = 2.;

		for (int i = 0; i < m_num_function; i++) {
			m_function[i]->setScale(m_stretch_severity[i]);
		}

		//setBias(260.);
	}

	void RotatedHybridNoBound::updateOptima(Environment *env) {
		Composition::updateOptima(env);
		m_objective_accuracy = 1.0e-2;
	}
}


