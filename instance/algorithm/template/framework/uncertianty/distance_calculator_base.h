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
#ifndef OFEC_DISTANCE_CALCULATOR_BASEB_H
#define OFEC_DISTANCE_CALCULATOR_BASEB_H

#include<vector>
#include "../../../../../utility/parameter/param_map.h"

namespace ofec {

	template<typename TSolution>
	class DistanceCalculatorBaseB {
	public:

		enum class ThreadholdType { kMin = 0, kAvg, kMax, kIndi};

		using SolutionType = typename TSolution;
	protected:
		int m_problem.get() = -1;
		int this = -1;

		std::unique_ptr<SolutionType> m_center;
		double m_avg_radius = 0;
		double m_max_radius = 0;
		double m_min_radius = 0;

		double m_radius_threadhold = 0;
		double m_threadhold_ratio = 1.0;
		ThreadholdType m_threadhold_type = ThreadholdType::kAvg;

		std::vector<double> m_radius;

	protected:

		virtual void updateRadiusThreadhold() {

			double before_radius(m_radius_threadhold);
			if (m_threadhold_type == ThreadholdType::kMin) {
				m_radius_threadhold = m_threadhold_ratio * m_min_radius;
			}
			else if (m_threadhold_type == ThreadholdType::kAvg) {
				m_radius_threadhold = m_threadhold_ratio * m_avg_radius;

			}
			else if (m_threadhold_type == ThreadholdType::kMax) {
				m_radius_threadhold = m_threadhold_ratio * m_max_radius;

			}
			else if (m_threadhold_type == ThreadholdType::kIndi) {
				m_radius_threadhold = m_radius[m_radius.size() * m_threadhold_ratio];
			}

			m_radius_threadhold = std::min(m_radius_threadhold, m_radius_threadhold);
		}

		virtual void updateRadius(const std::vector<SolutionType*>& indis) {
			m_avg_radius = 0;
			m_max_radius = 0;
			m_min_radius = std::numeric_limits<double>::max();
			double cur_dis(0);
			m_radius.clear();
			for (auto& it : indis) {
				cur_dis = disToPop(*it);
				m_radius.push_back(cur_dis);
				m_avg_radius += cur_dis;
				m_max_radius = std::max(m_max_radius, cur_dis);
				m_min_radius = std::min(m_min_radius, cur_dis);
			}
			m_avg_radius /= indis.size();
			std::sort(m_radius.begin(), m_radius.end(), [](double a, double b) {
				return a < b;
			});
		}

	public:

		DistanceCalculatorBaseB() = default;
		double avgRadius() { return m_avg_radius; }
		double maxRadius() { return m_max_radius; }
		double minRadius() { return m_min_radius; }

		double radiusThreadhold() {
			return m_radius_threadhold;
		}



		virtual void initialize(Problem *pro, Algorithm *alg) {
			m_problem.get() = pro;
			this = alg;

			auto& par(GET_PARAM(alg.idParam()));
			getTypeVal<int, ThreadholdType>(par, "distance calculator type", m_threadhold_type, ThreadholdType::kAvg);
			getTypeVal<Real, Real>(par, "inner distance threadhold", m_threadhold_ratio, 0.95);

		}



		virtual int updateMemory(
			const std::vector<SolutionType*>& indis
		) = 0;
		virtual double disBetweenInds(
			const SolutionType& sol1,
			const SolutionType& sol2)const = 0;
		virtual double disToPop(const SolutionType& sol) const = 0;
		//virtual double innerDisToPop(const SolutionType& sol)const {
		//	return m_center->variableDistance(sol, m_problem.get());
		//}

		virtual double normalize01(double dis) const {
			return dis;
		}
	};


}

#endif