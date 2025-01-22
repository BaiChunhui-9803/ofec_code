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

#include "f9_composition1.h"
#include "../../global/classical/schwefel.h"
#include "../../global/classical/rastrigin.h"
#include "../../global/classical/hg_bat.h"
#include "../../../../../../core/problem/solution.h"

namespace ofec::cec2015 {
	void F9_composition1::initialize_(Environment *env) {
		Continuous::initialize_();
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;

		auto& v = *m_param;;
		resizeVariable(std::get<int>(m_param->at("number of variables")));
		setDomain(-100., 100.);

		setFunction();
		loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");  //data path
		loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");
		for (auto &i : m_function) {
			i->setGlobalOpt(i->translation().data());
		}

		// Set optimal solution
		auto new_optima = new Optima<>();
		new_optima->appendSolution(m_function[0]->optima()->solution(0));
		m_optima.reset(new_optima);
	}

	void F9_composition1::setFunction() {
		m_num_function = 3;
		m_function.resize(m_num_function);
		m_height.resize(m_num_function);
		m_converge_severity.resize(m_num_function);
		m_f_bias.resize(m_num_function);
		m_param_fun.resize(m_num_function);

		BasicFunctions f(3);
		f[0] = &createFunction<Schwefel>;
		f[1] = &createFunction<Rastrigin>;
		f[2] = &createFunction<HGBat>;

		for (size_t i = 0; i < m_num_function; ++i) {
			auto param = *m_param;
			param["problem name"] = m_name + "_part" + std::to_string(i + 1);
			m_param_fun[i] = std::make_shared<ParameterMap>(param);
			m_function[i].reset(dynamic_cast<Function*>(f[i]()));
			m_function[i]->setRandom(m_random);
			m_function[i]->getInputParameters().input(*m_param_fun[i]);
			m_function[i]->initialize();
			m_function[i]->setBias(0);
		}
		m_function[0]->setRotated(false);

		for (auto &i : m_function)
			i->setConditionNumber(2.);

		m_converge_severity[0] = 20;
		m_converge_severity[1] = 20;
		m_converge_severity[2] = 20;

		m_height[0] = 1;
		m_height[1] = 1;
		m_height[2] = 1;

		m_f_bias[0] = 0;
		m_f_bias[1] = 100;
		m_f_bias[2] = 200;
	}

	void F9_composition1::evaluateObjective(Real *x, std::vector<Real> &obj) {
		Composition::evaluateObjective(x, obj);
		obj[0] += 900;
	}
}
