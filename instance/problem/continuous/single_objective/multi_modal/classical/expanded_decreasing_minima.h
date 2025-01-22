/********* Begin Register Information **********
{
	"name": "Classic_decreasing_minima_expanded",
	"identifier": "ExpandedDecreasingMinima",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper: Problem Definitions and Evaluation Criteria for the CEC 2015
*  Competition on Single Objective Multi-Niche Optimization.
*******************************************************************************************/

#ifndef OFEC_EXPANDED_DECREASING_MINIMA_H
#define OFEC_EXPANDED_DECREASING_MINIMA_H

#include "../../function.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	class ExpandedDecreasingMinima : public Function, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(ExpandedDecreasingMinima)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateOriginalObj(Real *x, std::vector<Real> &obj) override;
	};
}
#endif // !OFEC_EXPANDED_DECREASING_MINIMA_H



