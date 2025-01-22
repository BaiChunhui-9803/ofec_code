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

*************************************************************************
   Reference:
   Deb, K., Thiele, L., Laumanns, M., & Zitzler, E.
   Scalable multi-objective optimization test problems.
   In Evolutionary Computation, 2002. CEC'02. Proceedings of the 2002 Congress on IEEE.
*******************************************************************************/
// Created: 31 December 2014
// Modified: 29 Mar 2018 by Junchen Wang (wangjunchen@cug.edu.cn)


#ifndef OFEC_ZDT_H
#define OFEC_ZDT_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../problem/continuous/multi_objective/metrics_mop.h"

namespace ofec {
	class ZDT : public Continuous, public MetricsMOP  {
	protected:
		void initialize_() override;
		size_t m_num_reference_points;
		/*void updateOptima() override;
		void loadParetoFront(size_t sample_num);
		void sampleParetoSets(size_t sample_num);
		void sampleParetoFront(size_t sample_num);*/
		//void generateAdLoadPF();
	};
}

#endif //OFEC_ZDT_H