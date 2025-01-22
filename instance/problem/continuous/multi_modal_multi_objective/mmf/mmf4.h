/********* Begin Register Information **********
{
	"name": "MMOP_MMF4",
	"identifier": "MMF4",
	"problem tags": [ "MMOP","MOP", "ConOP" ]
}
*********** End Register Information **********/

/*********************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
**********************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@gmail.com
* Language: C++
**********************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

**********************************************************************************
  Reference:
  Yue C, Qu B, Liang J.
  A multiobjective particle swarm optimizer using ring topology for solving multimodal
  multiobjective problems[J].
  IEEE Transactions on Evolutionary Computation, 2017, 22(5): 805-817.
**********************************************************************************/


/*
   GLT1: 2 objs, PS is a sin curve but is continuous, dims diff, PF is two discontinuous lines
   GLT2: 2 objs, PS is a sin curve, dims diff, PF is convex
   GLT3: 2 objs, PS is a sin curve, dims diff, PF is extremely convex
   GLT4: 2 objs, PS is sin curve, dims diff, PS and PF are all discontinuous
   GLT5: 3 objs, PS is a sin nonlinear surface, dims diff, PF is convex
   GLT6: 3 objs, PS is a sin nonlinear surface, dims diff, PS and PF are all discontinuous
*/

#ifndef OFEC_MMF4_H
#define OFEC_MMF4_H

#include "../metrics_mmop.h"

namespace ofec {
	class MMF4 : public MetricsMMOP {
	protected:
		void initialize_() override;
		void updateOptima() override;
		void loadParetoFront(size_t sample_num);
		void sampleParetoSets(size_t sample_num);
		void sampleParetoFront(size_t sample_num);
		std::vector<Real> createVar(const std::vector<Real>& s) const override;

		void evaluateObjective(Real* x, std::vector<Real>& obj) override;

		/*void loadParetoFront();
		void generateOptimalSolution();*/
	};
}

#endif

