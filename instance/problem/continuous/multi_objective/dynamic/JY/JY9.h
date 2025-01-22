/********* Begin Register Information **********
{
	"name": "JY09",
	"identifier": "JY9",
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
Evolutionary dynamic multiobjective optimization: Benchmarks and algorithm comparisons.
IEEE transactions on cybernetics, 47(1):198-211, 2016
************************************************************************/

// Created: 5 August 2019 by Qingshan Tan
// Last modified at 5 August 2019 by Qingshan Tan


#ifndef JY9_H
#define JY9_H


#include "../DMOPs.h"
#include"../metrics_dmop.h"

namespace ofec {
	class JY9 : public DMOPs, public MetricsDMOP {
	public:
		void initialize_() override;
		void generateAdLoadPF() override;
		void set_sigma(int v) { m_sigma = v; }
		int get_sigma() { return m_sigma; }
	private:
		int updateEvaluationTag(SolutionBase& s, Algorithm *alg) override;
		void evaluateObjective(Real* x, std::vector<Real>& obj);
		int m_sigma = 0;
	};
}

#endif //JY9_H

