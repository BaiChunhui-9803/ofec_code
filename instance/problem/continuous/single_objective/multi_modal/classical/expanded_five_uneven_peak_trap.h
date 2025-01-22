/********* Begin Register Information **********
{
	"name": "Classic_five_uneven_peak_trap_expanded",
	"identifier": "ExpandedFiveUnevenPeakTrap",
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

#ifndef OFEC_EXPANDED_FIVE_UNEVEN_PEAK_TRAP_H
#define OFEC_EXPANDED_FIVE_UNEVEN_PEAK_TRAP_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {
	class ExpandedFiveUnevenPeakTrap : public Continuous, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(ExpandedFiveUnevenPeakTrap)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
	};
}
#endif // !OFEC_EXPANDED_FIVE_UNEVEN_PEAK_TRAP_H

