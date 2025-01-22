/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
*
*********************************************************************************/

#ifndef METRIC_UNCERTAINTY_H
#define METRIC_UNCERTAINTY_H

#include<vector>
#include"../../../core/problem/problem.h"

namespace ofec {
#define CAST_MetricComUncertianty(pro) dynamic_cast<MetricComUncertianty*>(pro)

	class MetricComUncertianty : virtual public Problem {

	protected:
		Real m_ratio = 0.5;
		int m_metrix_idRnd = -1;

	public:
	    virtual const std::vector<int>& getNeighborSize() = 0;
		virtual Real getEffectiveValue(Real cur_obj, const std::vector<std::vector<Real>>& nei_objs) {
			Real eff_val(cur_obj* m_ratio);
			Real ratio = (1.0- m_ratio)/ nei_objs.size();
			
			Real sum_obj(0);
			for (int idx(0); idx < nei_objs.size(); ++idx) {
				sum_obj = 0;
				for (auto& it : nei_objs[idx]) {
					sum_obj += it;
				}
				eff_val += sum_obj * ratio;
			}
			return eff_val;
		}







	protected:
	};
}

#endif 