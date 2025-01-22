/********* Begin Register Information **********
{
	"name": "MOP_real_DBD",
	"identifier": "Disc_brake_design",
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

#ifndef DISC_BRAKE_DESIGN_H
#define DISC_BRAKE_DESIGN_H

#include"../metrics_mop.h"
//#include"../oomop/components/mpb_class.h"

namespace ofec {
	class Disc_brake_design : public MetricsMOP {
	protected:
		void initialize_() override;
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
		void generateAdLoadPF();

	private:

	};
}

#endif //DISC_BRAKE_DESIGN_H