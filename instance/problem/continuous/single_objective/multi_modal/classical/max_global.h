/********* Begin Register Information **********
{
	"name": "Classic_max_global",
	"identifier": "MaxGlobal",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

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
// Created: 21 July 2011
// Last modified:
#ifndef OFEC_MAX_GLOBAL_H
#define OFEC_MAX_GLOBAL_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../single_objective/metrics_mmop.h"

namespace ofec {	
	class MaxGlobal : public Continuous, public MetricsMMOP {
		OFEC_CONCRETE_INSTANCE(MaxGlobal)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
	};
}

#endif // FMAX_GLOBAL1_H
