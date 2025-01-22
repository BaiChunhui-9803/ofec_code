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
#ifndef OFEC_DISTANCE_CALCULATOR_WEIGHT_H
#define OFEC_DISTANCE_CALCULATOR_WEIGHT_H


#include "../../template/framework/uncertianty/distance_calculator.h"
#include "fitness_weight_mapper.h"

namespace ofec {

	template<typename TSolution>
	class DistanceCalculatorWeight :public DistanceCalculatorBase<TSolution> {
	public:
		using SolutionType = typename TSolution;
	protected:

		std::vector<SolutionType> m_indis;
		std::vector<double> m_weight;
		FitnessWeightMapper m_fitness_weight;
		//double m_sum_weight = 0;

		double m_weight_scale = 0.98;

	protected:

	public:



		DistanceCalculatorWeight() = default;
		virtual void initialize(Problem *pro, Algorithm *alg) override {
			DistanceCalculatorBase::initialize(pro, alg);
			m_fitness_weight.reset();
		}
		virtual int updateMemory(
			const std::vector<SolutionType*>& indis
		) override {
			m_indis.resize(indis.size());
			for (int idx(0); idx < m_indis.size(); ++idx) {
				m_indis[idx] = *indis[idx];
			}

			m_weight.resize(indis.size());
			for (int idx(0); idx < indis.size(); ++idx) {
				m_weight[idx] = indis[idx]->fitness();
			}
			m_fitness_weight.reset();
			m_fitness_weight.update(m_weight);
			for (int idx(0); idx < indis.size(); ++idx) {
				m_weight[idx] = m_fitness_weight.getPopFitness(m_weight[idx]);
			}

			double sum_weight = 0;
			for (auto& it : m_weight) {
				sum_weight += it;
			}
			for (auto& it : m_weight) {
				it /= sum_weight;
			}

			updateRadius(indis);
			updateRadiusThreadhold();
			//m_radius_threadhold = m_max_radius*m_weight_scale;
			return 0;
		}
		virtual double disBetweenInds(
			const SolutionType& sol1,
			const SolutionType& sol2)const override {
			return sol1.variableDistance(sol2, m_problem.get());
		}
		virtual double disToPop(const SolutionType& sol) const override {
			
			double dis(0);
			for (int idx(0); idx < m_indis.size(); ++idx) {
				dis+= sol.variableDistance(m_indis[idx], m_problem.get())*m_weight[idx];
			}
			//dis /= m_indis.size();
			return dis;
		}
		//virtual double innerDisToPop(const SolutionType& sol)const override {
		//	return disToPop(sol);
		//}

		virtual double normalize01(double dis) const override {
			return dis;
		}
	};
}

#endif