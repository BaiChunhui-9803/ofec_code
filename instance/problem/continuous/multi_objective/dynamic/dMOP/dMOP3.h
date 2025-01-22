/********* Begin Register Information **********
{
	"name": "dMOP3",
	"identifier": "dMOP3",
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
Chi-Keong Goh and Kay Chen Tan (2009).
A Competitive-Cooperative Coevolutionary Paradigm for Dynamic Multiobjective Optimization
IEEE Transactions on evolutionary computation, 13(1), 103-127.
************************************************************************/

// Created: 5 August 2019 by Qingshan Tan
// Last modified at 5 August 2019 by Qingshan Tan


#ifndef dMOP3_H
#define dMOP3_H


#include "../DMOPs.h"
#include"../metrics_dmop.h"

namespace ofec {
	class dMOP3 : public DMOPs, public MetricsDMOP {
	public:
		void initialize_() override;
		void generateAdLoadPF() override;
		void set_r(size_t v) { m_r = v; }
		size_t get_r() { return m_r; }
	private:
		int updateEvaluationTag(SolutionBase& s, Algorithm *alg) override;
		void evaluateObjective(Real* x, std::vector<Real>& obj);
		size_t m_r;//the random dimension index
	};
}

#endif //dMOP3_H

