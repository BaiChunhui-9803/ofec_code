/********* Begin Register Information **********
{
	"name": "CEC2018_1",
	"identifier": "CEC2018_1",
	"problem tags": [ "DMOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Qingshan Tan
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

* see https://github.com/Changhe160/OFEC for more information
*************************************************************************/

/************************************************************************
Jiang S, Yang S, Yao X, et al. 
Benchmark Problems for CEC2018 Competition on Dynamic Multiobjective Optimisation[J]. 
Proc. CEC2018 Competition, 2018: 1-8.
************************************************************************/

// Created: 22 September 2020 by Qingshan Tan


#ifndef CEC2018_1_H
#define CEC2018_1_H

#include "../DMOPs.h"
#include"../metrics_dmop.h"

namespace ofec {
	class CEC2018_1 : public DMOPs, public MetricsDMOP {
	public:
		void initialize_() override;
		void generateAdLoadPF() override;
	private:
		int updateEvaluationTag(SolutionBase& s, Algorithm *alg) override;
		void evaluateObjective(Real* x, std::vector<Real>& obj);
	};
}

#endif //CEC2018_1_H


