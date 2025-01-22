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
#ifndef OFEC_EVALUATION_STRATEGY_HLS_H
#define OFEC_EVALUATION_STRATEGY_HLS_H



#include<vector>
#include"../../../../instance/problem/combination/selection_problem/selection_problem.h"
#include"../../template/framework/dynamic_noisy/evaluation_strategy.h"
#include"../../combination/sequence/SP/sp_interpreter.h"
namespace ofec {

	template<typename TSolution>
	class EvaluationStrategyHLS : EvaluationStrategyBase<TSolution> {
	public:
		using SolutionType = typename TSolution;
	protected:
		int m_hls_neighbor_size = 100;
	public:
		virtual int evaluate(SolutionType& curSol, int popIdx = 0){
			std::vector<std::unique_ptr<SolutionBase>> samples(m_hls_neighbor_size);
			GET_DYN_SP(m_problem.get())->generate_HLS_samples(curSol, samples, m_random.get());
			int rf(0);
			double obj(0);
			for (int idx(0); idx < m_hls_neighbor_size; ++idx) {
				rf|= samples[idx]->evaluate(m_problem.get(), this, true);
				obj += samples[idx]->objective()[0];
			}
			obj /= m_hls_neighbor_size;
			rf|=curSol.evaluate(m_problem.get(), this, true);
			curSol.setFitness(obj * 0.5 + curSol.objective()[0] * 5);
			curSol.setFitness(objectiveToFitness(curSol, curSol.fitness()));
			return rf;
		}

		virtual int calFitness(
			std::vector<SolutionType*>& indis,
			std::vector<SolutionType>& curSol,
			int popIdx = 0
		) override{
			int rf(0);
			for (auto& it : curSol) {
				rf |= evaluate(it);
				it.setFitness(objectiveToFitness(it.objective()[0]));
			}

			if (m_problem->hasTag(ProblemTag::kDOP)) {
				for (auto& it : indis) {
					rf |= evaluate(*it);
					it->setFitness(objectiveToFitness(it.objective()[0]));
				}
			}
			return rf;
		}

		//virtual void update(Problem *pro, Algorithm *alg, Random *rnd, const SolutionType& otherSol, SolutionType& curSol)const {}

	};
}

#endif