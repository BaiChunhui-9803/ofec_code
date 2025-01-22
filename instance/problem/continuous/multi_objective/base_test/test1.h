/********* Begin Register Information **********
{
	"name": "MOP_TEST1",
	"identifier": "TEST1",
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
// Created: 9 Feb. 2023
// two obj, unimodal function

#ifndef TEST1_H
#define TEST1_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"

namespace ofec {
	class TEST1 : public Continuous, public MetricsMOP  {
	protected:
		void initialize_() override;
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
		void generateAdLoadPF();

	private:
		//peaks' position is on a circle
		std::vector<std::vector<Real>> m_peaks_pos;
		std::vector<Real> m_center_pos;
		Real m_radius;
		std::vector<std::vector<Real>> m_angles;//单个子目标在圆上的角度
		bool m_rand_flag;
		std::vector<Real> m_heights;
		std::vector<Real> m_slopes;
	};
}

#endif //TEST1_H