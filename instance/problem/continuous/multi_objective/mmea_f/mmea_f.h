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
  Zhou A, Zhang Q, Jin Y. Approximating the set of Pareto-optimal solutions in 
  both the decision and objective spaces by an estimation of distribution algorithm[J]. 
  IEEE transactions on evolutionary computation, 2009, 13(5): 1167-1189.
**********************************************************************************/
// added by qingshan tan in Jan 20th, 2024
// recommend n=20

/*
  MMEA_F3: 2 objs, PS is a 2-D sin nonlinear surface, PF is convex
  MMEA_F4: 2 objs, PS is a 2-D sin nonlinear surface, PF is concave
  MMEA_F5: 2 objs, PS is a 2-D sin nonlinear surface, PF is mix and discontinuous
  MMEA_F6: 2 objs, PS is a 3-D continuous sin nonlinear manifold, PF is concave
  MMEA_F7: 3 objs, PS is a 3-D continuous sin nonlinear manifold, PF is concave
*/

#ifndef OFEC_MMEA_F_H
#define OFEC_MMEA_F_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../problem/continuous/multi_objective/metrics_mop.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/nondominated_sorting/filter_sort.h"

namespace ofec {
	class MMEA_F : public Continuous,public MetricsMOP {
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