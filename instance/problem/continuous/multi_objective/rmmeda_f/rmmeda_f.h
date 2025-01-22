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
  Zhang Q, Zhou A, Jin Y. 
  RM-MEDA: A regularity model-based multiobjective estimation of distribution algorithm[J]. 
  IEEE Transactions on Evolutionary Computation, 2008, 12(1): 41-63.
**********************************************************************************/
// added by qingshan tan in Jan 20th, 2024
// recommend n=30

/*
   RMMEDA_F1: 2 objs, PS is a line, dims equal, PF is convex
   RMMEDA_F2: 2 objs, PS is a line, dims equal, PF is concave
   RMMEDA_F3: 2 objs, PS is a line, dims equal, PF is concave
   RMMEDA_F4: 3 objs, PS is a 2-D linear surface, dims equal, PF is a sphere, concave
   RMMEDA_F5: 2 objs, PS is a sqrt nonlinear curve, dims equal, PF is convex
   RMMEDA_F6: 2 objs, PS is a sqrt nonlinear curve, dims equal, PF is concave
   RMMEDA_F7: 2 objs, PS is a sqrt nonlinear curve, dims equal, PF is concave, has bias
   RMMEDA_F8: 3 objs, PS is a 2-D sqrt nonlinear surface, dims equal, PF is a sphere, concave
   RMMEDA_F9: 2 objs, PS is a sqrt nonlinear curve, dims equal, PF is convex, has bias, g is multi-modal
   RMMEDA_F10: 2 objs, PS is a sqrt nonlinear curve, dims equal, PF is convex, has bias, g is multi-modal
*/

#ifndef OFEC_RMMEDA_F_H
#define OFEC_RMMEDA_F_H

#include "../../../../problem/continuous/multi_objective/metrics_mop.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/nondominated_sorting/filter_sort.h"

namespace ofec {
	class RMMEDA_F : public Continuous, public MetricsMOP {
	
	protected:
		void initialize_() override;
		/*void updateOptima() override;
		void loadParetoFront(size_t sample_num);
		void sampleParetoSets(size_t sample_num);
		void sampleParetoFront(size_t sample_num);*/
		std::vector<Real> createVar(const std::vector<Real>& s) const override;
		size_t m_num_reference_points;
		//void loadParetoSets();
		//void generateOptimalSolution();
	};
}

#endif