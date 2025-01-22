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

#include "f21_rotated_hybrid_composition.h"
#include "../classical/griewank_rosenbrock.h"
#include "../classical/griewank.h"
#include "../classical/rastrigin.h"
#include "../classical/weierstrass.h"
#include "../classical/rotated_scaffer_F6.h"

namespace ofec::cec2005 {
	void RotatedHybridComposition3::setFunction(Environment *env) {
		m_num_function = 10;
		m_function.resize(m_num_function);
		m_param_fun.resize(m_num_function);
		m_height.resize(m_num_function);
		m_fmax.resize(m_num_function);
		m_converge_severity.resize(m_num_function);
		m_stretch_severity.resize(m_num_function);

		BasicFunctions f(5);
		f[0] = &createFunction<RotatedScafferF6>;
		f[1] = &createFunction<Rastrigin>;
		f[2] = &createFunction<GriewankRosenbrock>;
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
		m_function[0]->loadRotation("/instance/problem/continuous/single_objective/global/classical/GOP_CEC2005_F21_part1");
		m_function[1]->loadRotation("/instance/problem/continuous/single_objective/global/classical/GOP_CEC2005_F21_part2");

		for (auto &i : m_function)
			i->setConditionNumber(2.);

		for (int i = 0; i < m_num_function; i++) {
			m_height[i] = 100 * i;
		}

		m_function[0]->setDomain(-100, 100);   m_function[1]->setDomain(-100, 100);
		m_function[2]->setDomain(-5, 5);     m_function[3]->setDomain(-5, 5);
		m_function[4]->setDomain(-5, 5); m_function[5]->setDomain(-5, 5);
		m_function[6]->setDomain(-0.5, 0.5); m_function[7]->setDomain(-0.5, 0.5);
		m_function[8]->setDomain(-200, 200); m_function[9]->setDomain(-200, 200);

		m_stretch_severity[0] = 5.*5. / 100; m_stretch_severity[1] = 5. / 100;
		m_stretch_severity[2] = 5.;		m_stretch_severity[3] = 1.;
		m_stretch_severity[4] = 5;  m_stretch_severity[5] = 1;
		m_stretch_severity[6] = 50.;	m_stretch_severity[7] = 10.;
		m_stretch_severity[8] = 5. * 5 / 200;  m_stretch_severity[9] = 5. / 200;

		m_converge_severity[0] = 1.; m_converge_severity[1] = 1.;
		m_converge_severity[2] = 1.;	m_converge_severity[3] = 1.;
		m_converge_severity[4] = 1.;  m_converge_severity[5] = 2.;
		m_converge_severity[6] = 2.;	m_converge_severity[7] = 2.;
		m_converge_severity[8] = 2.;  m_converge_severity[9] = 2.;

		for (int i = 0; i < m_num_function; i++) {
			m_function[i]->setScale(m_stretch_severity[i]);
		}
		//setBias(360.);
	}

	void RotatedHybridComposition3::addInputParameters() {
		m_index_number = "21";
	}

	void RotatedHybridComposition3::evaluateObjective(Real *x, std::vector<Real> &obj) const {
		Composition::evaluateObjective(x, obj);
		obj[0] += 360.;
	}

	void RotatedHybridComposition3::updateOptima(Environment *env) {
		Composition::updateOptima(env);
		m_objective_accuracy = 1.0e-2;
	}
}
