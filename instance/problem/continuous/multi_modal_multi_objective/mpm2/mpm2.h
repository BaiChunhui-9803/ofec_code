/********* Begin Register Information **********
[
	{ "name":"MPM2", "identifier":"MPM2", "problem tags":["ConOP","GOP","MOP"] }
]
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yaqi Ti
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// source code refered to https://github.com/jakobbossek/smoof/blob/master/inst/mpm2.py

#ifndef OFEC_MPM2_H
#define OFEC_MPM2_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"
#include"mpm.h"

namespace ofec {
	class MPM2 : public Continuous, public MetricsMOP {
	public:
		std::vector<std::unique_ptr<MultiplePeaksModel>> mpm_objs;

		void initialize_() override;
		void evaluateObjective(Real* x, std::vector<Real>& obj);
		void compute_opt();
		//void compute_Fopt();

	};
}


#endif // !OFEC_MPM2_H
