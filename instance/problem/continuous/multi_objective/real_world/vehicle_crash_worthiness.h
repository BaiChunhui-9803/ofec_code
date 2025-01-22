/********* Begin Register Information **********
{
	"name": "MOP_real_VCW",
	"identifier": "Vehicle_crash_worthiness",
	"problem tags": ["MOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 1 Apr. 2023
// Modified: 
// the landscape is diversed from the corner of the bound 
// two objs, a single multi-objective landscape structure in search space

#ifndef VEHICLE_CRASH_WORTHINESS_H
#define VEHICLE_CRASH_WORTHINESS_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"

namespace ofec {
	class Vehicle_crash_worthiness : public Continuous, public MetricsMOP  {
	protected:
		void initialize_() override;
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
		void generateAdLoadPF();

	private:

	};
}

#endif //VEHICLE_CRASH_WORTHINESS_H