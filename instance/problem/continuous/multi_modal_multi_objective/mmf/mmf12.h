/********* Begin Register Information **********
{
	"name": "MMOP_MMF12",
	"identifier": "MMF12",
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
   MMF12: 2 objs, 2 vars, one global PS, other are local PS. PS is a sin curve but is continuous, dims diff, PF is two discontinuous lines
*/

#ifndef OFEC_MMF12_H
#define OFEC_MMF12_H

#include "../metrics_mmop.h"

namespace ofec {
	class MMF12 : public MetricsMMOP {
	public:
		size_t m_num_ps;
		size_t m_num_pf;
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

