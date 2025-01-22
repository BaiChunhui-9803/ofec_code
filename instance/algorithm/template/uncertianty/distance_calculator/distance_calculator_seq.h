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
#ifndef OFEC_DISTANCE_CALCULATOR_MAT_H
#define OFEC_DISTANCE_CALCULATOR_MAT_H

#include "../../../../utility/functional.h"
#include "../../template/framework/uncertianty/distance_calculator.h"
#include "fitness_weight_mapper.h"

// for test


namespace ofec {
	template<typename TSolution>
	class DistanceCalculatorSeqMat :public DistanceCalculatorBase<TSolution> {
	public:



		using SolutionType = typename TSolution;
		std::vector<std::vector<double>> m_density_matrix;
		//	std::vector<std::vector<double>> m_temp_matrix;
		//	double m_sum_fitness = 0;
		std::vector<int> m_mat_size;
		double m_num_var = 0;

		std::vector<double> m_weight;
		FitnessWeightMapper m_fitness_weight;

		double m_Emax_radius = 0;
		double m_Emin_radius = 0;

	//	double m_thread_ratio = 0.88;

	

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

		virtual void initialize(Problem *pro, Algorithm *alg) override {
			DistanceCalculatorBase::initialize(pro, alg);
			m_num_var = pro->numberVariables();
			auto& matSize(GET_ASeq(alg)->interpreter().getMatrixSize());
			UTILITY::assignVVector<double>(m_density_matrix, matSize, 0);
			m_mat_size = matSize;
		}


		virtual int updateMemory(
			const std::vector<SolutionType*>& indis
		)override {

			std::vector<double> fitness;
			m_weight.resize(indis.size());
			for (int idx(0); idx < indis.size(); ++idx) {
				m_weight[idx] = indis[idx]->fitness();
			}
			fitness = m_weight;
			m_fitness_weight.reset();
			m_fitness_weight.update(m_weight);
			for (int idx(0); idx < indis.size(); ++idx) {
				m_weight[idx] = m_fitness_weight.getPopFitness(m_weight[idx]);
			}


			UTILITY::assignVVector<double>(m_density_matrix, 0);
			for (int idx(0); idx < indis.size(); ++idx) {
				for (auto& edge : indis[idx]->edges()) {
					m_density_matrix[edge.first][edge.second] += m_weight[idx];
				}
			}

			double max_val(0);
			for (auto& it : m_density_matrix) {
				for (auto& it2 : it) {
					max_val = std::max(it2, max_val);
				}
			}
			m_Emax_radius = m_num_var * (max_val);
			m_Emin_radius = 0;

			std::vector<double> radius(indis.size());
			for (int idx(0); idx < indis.size(); ++idx) {
				radius[idx] = disToPop(*indis[idx]);
			}

			updateRadius(indis);
			
			updateRadiusThreadhold();
			//m_radius_threadhold = m_max_radius * m_thread_ratio;
			return 0;

		}
		double disBetweenInds(const SolutionType& sol1, const SolutionType& sol2)const override {

			static std::vector<std::vector<bool>> temp_mat;

			UTILITY::assignVVector<bool>(temp_mat, m_mat_size, false);

			for (auto& edge : sol1.edges()) {
				temp_mat[edge.first][edge.second] = true;
			}
			double totalFit(0);
			for (auto& edge : sol2.edges()) {
				if (temp_mat[edge.first][edge.second]) {
					totalFit += temp_mat[edge.first][edge.second];
				}
			}
			return  m_Emax_radius - totalFit;
		}
		virtual double disToPop(const SolutionType& sol) const override {
			double dis(0);
			for (auto& edge : sol.edges()) {
				dis += m_density_matrix[edge.first][edge.second];
			}
			return m_Emax_radius - dis;
		}
		//double innerDisToPop(const TIndi& sol)const {
		//	//auto & inter(GET_AS
		//	double dis(0);
		//	double sumFit(m_sum_fitness);
		//	for (auto& edge : sol.edges()) {
		//		dis += m_density_matrix[edge.first][edge.second];
		//	}
		//	dis = dis / (sumFit * m_num_var);
		//	return - dis;
		//}

		virtual double normalize01(double dis)const override {
			return mapReal<double>(dis, 0, m_Emax_radius - m_Emin_radius, 0, 1);
		}
	};
}

#endif