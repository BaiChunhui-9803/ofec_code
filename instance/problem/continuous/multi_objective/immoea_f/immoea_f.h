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
  R. Cheng, Y. Jin, K. Narukawa, and B. Sendhoff, A multiobjective
  evolutionary algorithm using Gaussian process-based inverse modeling,
  IEEE Transactions on Evolutionary Computation, 2015, 19(6): 838-856.
**********************************************************************************/
// added by qingshan tan in Jan 20th, 2024
// recommend n=30

/*
   IMMOEA_F1: 2 objs, PS is a line, dims diff, PF is convex
   IMMOEA_F2: 2 objs, PS is a line, dims diff, PF is concave
   IMMOEA_F3: 2 objs, PS is a line, dims diff, PF is concave
   IMMOEA_F4: 3 objs, PS is a 2-D linear surface, dims diff, PF is a sphere, concave
   IMMOEA_F5: 2 objs, PS is a power nonlinear curve, dims diff, PF is convex
   IMMOEA_F6: 2 objs, PS is a power nonlinear curve, dims diff, PF is concave
   IMMOEA_F7: 2 objs, PS is a power nonlinear curve, dims diff, PF is concave, has bias
   IMMOEA_F8: 3 objs, PS is a 2-D power nonlinear surface, dims diff, PF is a sphere, concave
   IMMOEA_F9: 2 objs, PS is a power nonlinear curve, dims diff, PF is convex, has bias, g is multi-modal
   IMMOEA_F10: 2 objs, PS is a power nonlinear curve, dims diff, PF is convex, has bias, g is multi-modal
*/

#ifndef OFEC_IMMOEA_F_H
#define OFEC_IMMOEA_F_H

#include "../../../../problem/continuous/multi_objective/metrics_mop.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/nondominated_sorting/filter_sort.h"

namespace ofec {
	class IMMOEA_F : public Continuous, public MetricsMOP {
	protected:
		Real m_alpha;
		Real m_beta;
		size_t m_num_reference_points;
	protected:
		void initialize_() override;
		/*void updateOptima() override;
		void loadParetoFront(size_t sample_num);
		void sampleParetoSets(size_t sample_num);
		void sampleParetoFront(size_t sample_num);*/
		std::vector<Real> createVar(const std::vector<Real>& s) const override;
		//void loadParetoSets();
		//void generateOptimalSolution();
	};
}

#endif