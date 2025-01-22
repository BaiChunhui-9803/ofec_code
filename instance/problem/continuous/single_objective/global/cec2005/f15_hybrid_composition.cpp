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

#include "f15_hybrid_composition.h"
#include "../classical/sphere.h"
#include "../classical/ackley.h"
#include "../classical/griewank.h"
#include "../classical/rastrigin.h"
#include "../classical/weierstrass.h"

namespace ofec {
	namespace cec2005 {
		void HybridComposition::addInputParameters() {
			m_index_number = "15";
		}

		void HybridComposition::setFunction(Environment *env) {
			m_num_function = 10;
			m_function.resize(m_num_function);
			m_param_fun.resize(m_num_function);
			m_height.resize(m_num_function);
			m_fmax.resize(m_num_function);
			m_converge_severity.resize(m_num_function);
			m_stretch_severity.resize(m_num_function);

			BasicFunctions f(5);
			f[0] = &createFunction<Rastrigin>;
			f[1] = &createFunction<Weierstrass>;
			f[2] = &createFunction<Griewank>;
			f[3] = &createFunction<Ackley>;
			f[4] = &createFunction<Sphere>;

			for (size_t i = 0; i < m_num_function; ++i) {
				m_param_fun[i].reset(new ParameterMap(m_archived_parameters));
				m_function[i].reset(dynamic_cast<Function*>(f[i / 2]()));
				m_function[i]->reset();
				m_function[i]->inputParameters().input(*m_param_fun[i]);
				m_function[i]->initialize(env);
				m_function[i]->setBias(0);
			}

			for (int i = 0; i < m_num_function; i++) {
				m_height[i] = 100 * i;
				m_converge_severity[i] = 1.;
			}

			m_function[0]->setDomain(-5, 5);     m_function[1]->setDomain(-5, 5);
			m_function[2]->setDomain(-0.5, 0.5); m_function[3]->setDomain(-0.5, 0.5);
			m_function[4]->setDomain(-60, 60);	 m_function[5]->setDomain(-60, 60);
			m_function[6]->setDomain(-32, 32);   m_function[7]->setDomain(-32, 32);
			m_function[8]->setDomain(-100, 100); m_function[9]->setDomain(-100, 100);

			m_stretch_severity[0] = 1.;		m_stretch_severity[1] = 1.;
			m_stretch_severity[2] = 10.;		m_stretch_severity[3] = 10.;
			m_stretch_severity[4] = 5. / 60;  m_stretch_severity[5] = 5. / 60;
			m_stretch_severity[6] = 5. / 32;	m_stretch_severity[7] = 5. / 32;
			m_stretch_severity[8] = 5. / 100;  m_stretch_severity[9] = 5. / 100;

			for (int i = 0; i < m_num_function; i++) {
				m_function[i]->setScale(m_stretch_severity[i]);
			}
		}

		void HybridComposition::evaluateObjective(Real *x, std::vector<Real> &obj) const {
			Composition::evaluateObjective(x, obj);
			obj[0] += 120.;   // add m_bias
		}

		void HybridComposition::setRotation() {
			for (auto i : m_function)
				i->rotation().identify();
		}
		
		void HybridComposition::updateOptima(Environment *env) {
			Composition::updateOptima(env);
			m_objective_accuracy = 1.0e-2;
		}
	}
}


