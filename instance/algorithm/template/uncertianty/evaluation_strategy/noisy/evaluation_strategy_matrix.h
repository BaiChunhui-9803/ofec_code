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
#ifndef OFEC_EVALUATION_STRATEGY_MATRIX_H
#define OFEC_EVALUATION_STRATEGY_MATRIX_H

#include<vector>
#include"../../template/combination/sequence/sequence_algorithm.h"
#include"../../../../utility/function/custom_function.h"

//#include"../../combination/sequence/SP/sp_interpreter.h"
namespace ofec {

	template<typename TSolution>
	class EvaluationStrategyMatrix :public EvaluationStrategyBase<TSolution> {
	public:
		using SolutionType = typename TSolution;
	protected:
		long long m_curTimes = 0;
		long long m_init_updateTimes = -100;
		double m_alpha = 0.5;
		double m_beta = 0.8;
		//std::vector<std::vector<bool>> m_edge_visited;
		std::vector<std::vector<double>> m_weight;
		std::vector<std::vector<long long>> m_updatedTimes;
		std::vector<std::vector<std::pair<double,int>>> m_temp;
		
	
	public:
		virtual void initialize(Problem *pro, Algorithm *alg)override {
			EvaluationStrategyBase::initialize(pro, alg);
			m_curTimes = 0;
			auto & matSize(GET_ASeq(alg)->interpreter().getMatrixSize());
			//UTILITY::assignVVector<bool>(m_edge_visited, matSize, false);
			UTILITY::assignVVector<double>(m_weight, matSize, 0);
			UTILITY::assignVVector<long long>(m_updatedTimes, matSize, m_init_updateTimes);
			
		}
		virtual void updateSolutions(
			const std::vector<const SolutionType*>& indis)override {
			++m_curTimes;
			auto& matSize(GET_ASeq(this).interpreter().getMatrixSize());
			UTILITY::assignVVector<std::pair<double, int>>(m_temp, matSize, std::pair<double, int>(0, 0));
			for (auto& it : indis) {
				for (auto& edge : it->edges()) {
					m_temp[edge.first][edge.second].first += it->objective()[0];
					++m_temp[edge.first][edge.second].second;
				}
			}
			
			for (int from(0); from < matSize.size(); ++from) {
				for (int to(0); to < matSize[from]; ++to) {
					if (m_temp[from][to].second > 0) {
						m_temp[from][to].first /= m_temp[from][to].second;
						double ratio = std::pow<double>(m_alpha, m_curTimes - m_updatedTimes[from][to]);
						m_updatedTimes[from][to] = m_curTimes;
						m_weight[from][to] = m_weight[from][to] * m_alpha + (1.0 - m_alpha) * m_temp[from][to].first;
						m_temp[from][to].first = 1.0;
					}
					else {
						m_temp[from][to].first = std::pow<double>(m_alpha, m_curTimes - m_updatedTimes[from][to]);
					}
				}
			}
		};



		virtual void updateFitness(SolutionType& curSol) override {
			double fitness(0);
			curSol.setFitness(curSol.objective()[0]);
			for (auto& edge : curSol.edges()) {
				double ratio = m_temp[edge.first][edge.second].first;
				fitness += m_weight[edge.first][edge.second] * ratio + (1.0 - ratio) * curSol.fitness();
			}
			fitness /= curSol.edges().size();
			curSol.setFitness(curSol.fitness() * (1.0 - m_beta) + m_beta * fitness);
			objectiveToFitness(curSol, curSol.fitness());
		}
		int evaluate(SolutionType& curSol, int popIdx = 0){
			return curSol.evaluate(m_problem.get(), this, true);
		}
	};
}

#endif