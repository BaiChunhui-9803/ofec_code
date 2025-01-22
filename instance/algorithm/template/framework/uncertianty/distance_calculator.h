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
* class Algorithm is an abstract for all algorithms.
*
*********************************************************************************/
#ifndef OFEC_DISTANCE_CALCULATOR_BASE_H
#define OFEC_DISTANCE_CALCULATOR_BASE_H

#include <vector>
#include "../../../../../utility/parameter/param_map.h"
#include "distance_calculator_base.h"

namespace ofec {

	template<typename TSolution>
	class DistanceCalculatorBase : public DistanceCalculatorBaseB<TSolution>{
	public:


		using SolutionType = typename TSolution;
	protected:


		std::unique_ptr<SolutionType> m_center;

	protected:


		//virtual void updateRadius(const std::vector<SolutionType*>& indis) {
		//	m_avg_radius = 0;
		//	m_max_radius = 0;
		//	m_min_radius = std::numeric_limits<double>::max();
		//	double cur_dis(0);
		//	for (auto& it : indis) {
		//		cur_dis = disToPop(*it);
		//		m_avg_radius += cur_dis;
		//		m_max_radius = std::max(m_max_radius, cur_dis);
		//		m_min_radius = std::min(m_min_radius, cur_dis);
		//	}
		//	m_avg_radius /= indis.size();


		//}

	public:





		virtual int updateMemory(
			const std::vector<SolutionType*>& indis
		) override{
			if (indis.empty()) {
				m_center.reset(nullptr);
				m_radius_threadhold = 0;
				return 0;
			}
			int center_idx(0);
			for (int idx(0); idx < indis.size(); ++idx) {
				if (indis[idx]->fitness() > indis[center_idx]->fitness()) {
					center_idx = idx;
				}
			}
			m_center.reset(new SolutionType(*indis[center_idx]));
			updateRadius(indis);

			updateRadiusThreadhold();
			//m_radius_threadhold = m_avg_radius* m_threadhold_ratio;
			return 0;
		}
		virtual double disBetweenInds(
			const SolutionType& sol1,
			const SolutionType& sol2)const override {
			return sol1.variableDistance(sol2, m_problem.get());
		}
		virtual double disToPop(const SolutionType& sol) const override {
			return m_center->variableDistance(sol, m_problem.get());
		}
		//virtual double innerDisToPop(const SolutionType& sol)const {
		//	return m_center->variableDistance(sol, m_problem.get());
		//}

	};


}

#endif