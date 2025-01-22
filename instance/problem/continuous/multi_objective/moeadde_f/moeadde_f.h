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
  H. Li and Q. Zhang, Multiobjective optimization problems with complicated
  Pareto sets, MOEA/D and NSGA-II, IEEE Transactions on Evolutionary
  Computation, 2009, 13(2): 284-302.
**********************************************************************************/
// added by qingshan tan in Jan 20th, 2024
// recommend n=10

/* 
   variables grouping
   MOEADDE_F1: 2 objs, PS is power nonlinear curve, dims diff, PF is convex
   MOEADDE_F2: 2 objs, PS is a 2-period sin nonlinear curve, dims diff, PF is convex
   MOEADDE_F3: 2 objs, PS is a sin nonlinear curve, dims diff, PF is convex
   MOEADDE_F4: 3 objs, PS is a sin nonlinear curve, dims diff, PF is convex
   MOEADDE_F5: 2 objs, PS is a complicated sin nonlinear curve, dims diff, PF is convex
   MOEADDE_F6: 3 objs, PS is a sin nonlinear curve, dims diff, PF is concave
   MOEADDE_F7: 2 objs, PS is a power nonlinear curve, dims diff, PF is convex, has local Pareto sets
   MOEADDE_F8: 2 objs, PS is a power nonlinear curve, dims diff, PF is convex, has local Pareto sets
   MOEADDE_F9: 2 objs, PS is a 3-period sin nonlinear curve, dims diff, PF is concave
*/

#ifndef OFEC_MOEADDE_F_H
#define OFEC_MOEADDE_F_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../problem/continuous/multi_objective/metrics_mop.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/nondominated_sorting/filter_sort.h"

namespace ofec {
	class MOEADDE_F : public Continuous,MetricsMOP {
	protected:
		void initialize_() override;
		//void updateOptima() override;
		//void loadParetoFront(size_t sample_num);
		//void sampleParetoSets(size_t sample_num);
		//void sampleParetoFront(size_t sample_num);
		std::vector<Real> createVar(const std::vector<Real>& s) const override;
		size_t m_num_reference_points;
		//void loadParetoSets();
		//void generateOptimalSolution();
	};
}

#endif