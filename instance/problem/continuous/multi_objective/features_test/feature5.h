/********* Begin Register Information **********
{
	"name": "MOP_feature5",
	"identifier": "Feature5",
	"problem tags": ["MOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 1 Apr. 2023
// Modified: 
// a local Pareto region which is derive from the boundary and extend to the middle of space 
// two objs, a single multi-objective landscape structure in search space

#ifndef FEATURE5_H
#define FEATURE5_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"
#include "../oomop/components/mpb_class.h"

namespace ofec {
	class Feature5 : public Continuous, public MetricsMOP  {
	protected:
		void initialize_() override;
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
		void generateAdLoadPF();

	private:
		int m_peak_num;
		Real m_max_height;
		std::vector<Real> m_peak_width;
		std::vector<Real> m_peak_angle;
		std::vector<Real> m_peak_norm;
		std::vector<std::shared_ptr<Mpb_class>> m_obj_mpb;
		std::vector<std::vector<std::vector<Real>>> m_peaks_pos;

		std::vector<std::vector<std::vector<Real>>> m_total_peaks_pos;
	};
}

#endif //FEATURE5_H