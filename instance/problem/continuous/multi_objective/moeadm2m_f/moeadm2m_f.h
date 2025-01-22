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
  H. Liu, F. Gu, and Q. Zhang, 
  Decomposition of a multiobjective optimization problem into a number of simple
  multiobjective subproblems,
  IEEE Transactions on Evolutionary Computation, 2014, 18(3): 450-455.
**********************************************************************************/
// added by qingshan tan in Jan 20th, 2024
// recommend n=10

/*
   MOEADM2M_F1: 2 objs, PS is a sin nonlinear curve about x1, dims equal, PF is convex
   MOEADM2M_F2: 2 objs, PS is a sin nonlinear curve about x1, dims equal, PF is concave
   MOEADM2M_F3: 2 objs, PS is a sin nonlinear curve about x1, dims equal, PF is concave
   MOEADM2M_F4: 2 objs, PS is a sin nonlinear curve about x1, dims equal, PS and PF is continuous
   MOEADM2M_F5: 2 objs, PS is a sin nonlinear curve about x1, dims equal, PF is convex
   MOEADM2M_F6: 3 objs, PS is a 2-D sin nonlinear manifold, dims equal, PF is a surface
   MOEADM2M_F7: 3 objs, PS is a 2-D sin nonlinear manifold, dims equal, PF is a sphere
*/

#ifndef OFEC_MOEADM2M_F_H
#define OFEC_MOEADM2M_F_H

#include "../../../../problem/continuous/multi_objective/metrics_mop.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/nondominated_sorting/filter_sort.h"

namespace ofec {
	class MOEADM2M_F : public Continuous, public MetricsMOP {
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