/********* Begin Register Information **********
{
	"name": "JY01",
	"identifier": "JY1",
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
Shouyong Jiang and Shengxiang Yang.
Evolutionary dynamic multiobjective optimization: Benchmarks and algorithm comparisons.
IEEE transactions on cybernetics, 47(1):198-211, 2016
************************************************************************/

// Created: 5 August 2019 by Qingshan Tan
// Last modified at 5 August 2019 by Qingshan Tan


#ifndef JY1_H
#define JY1_H


#include "../DMOPs.h"
#include"../metrics_dmop.h"

namespace ofec {
	class JY1 : public DMOPs, public MetricsDMOP {
	public:
		void initialize_() override;
		void generateAdLoadPF() override;
	private:
		int updateEvaluationTag(SolutionBase& s, Algorithm *alg) override;
		void evaluateObjective(Real* x, std::vector<Real>& obj);
	};
}

#endif //JY1_H

