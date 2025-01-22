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
  Gu F, Liu H L, Tan K C.
  A multiobjective evolutionary algorithm using dynamic weight design method[J].
  International Journal of Innovative Computing, Information and Control,2012,8(5(B)):3677-3688.
**********************************************************************************/
// recommend n=10

/*
   GLT1: 2 objs, PS is a sin curve but is continuous, dims diff, PF is two discontinuous lines
   GLT2: 2 objs, PS is a sin curve, dims diff, PF is convex
   GLT3: 2 objs, PS is a sin curve, dims diff, PF is extremely convex
   GLT4: 2 objs, PS is sin curve, dims diff, PS and PF are all discontinuous 
   GLT5: 3 objs, PS is a sin nonlinear surface, dims diff, PF is convex
   GLT6: 3 objs, PS is a sin nonlinear surface, dims diff, PS and PF are all discontinuous
*/

#ifndef OFEC_GLT_H
#define OFEC_GLT_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../problem/continuous/multi_objective/metrics_mop.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/nondominated_sorting/filter_sort.h"

namespace ofec {
	class GLT : public Continuous, public MetricsMOP  {
	protected:
		void initialize_() override;
		/*void updateOptima() override;
		void loadParetoFront(size_t sample_num);
		void sampleParetoSets(size_t sample_num);
		void sampleParetoFront(size_t sample_num);*/
		std::vector<Real> createVar(const std::vector<Real>& s) const override;
		
		size_t m_num_reference_points;
		/*void loadParetoFront();
		void generateOptimalSolution();*/
	};
}

#endif

