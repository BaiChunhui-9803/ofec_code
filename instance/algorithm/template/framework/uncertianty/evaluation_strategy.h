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
#ifndef OFEC_EVALUATION_STRATEGY_H
#define OFEC_EVALUATION_STRATEGY_H

#include<vector>

namespace ofec {

	template<typename TSolution> 
	class EvaluationStrategyBase {
	public:
		using SolutionType = typename TSolution;
	protected:
		int m_problem.get() = -1;
		int this = -1;
		int m_random.get() = -1;

		virtual double objectiveToFitness(double obj) {
			if (m_problem->optimizeMode(0) == OptimizeMode::kMinimize)
				return -obj;
			else
				return obj;
		}
	public:
		virtual void initialize(Problem *pro, Algorithm *alg) {
			m_problem.get() = pro;
			this = alg;
			m_random.get() = GET_ALG(this).idRandom();
		}

		virtual void clearMemory(int idPop) {}

		virtual int evaluate(SolutionType& curSol,int popIdx=0){
			int rf(0);
			rf = curSol.evaluate(m_problem.get(), this, true);
			curSol.setFitness(objectiveToFitness(curSol.objective()[0]));
			return rf;
		}

		virtual int calFitness(
			std::vector<SolutionType*>& indis,
			int popIdx = 0
		) {
			int rf(0);
			for (auto& it : indis) {
				rf |= evaluate(*it,popIdx);
			}
			return rf;
		}
		virtual int calFitness(
			std::vector<SolutionType*>& indis,
			std::vector<SolutionType>& curSol,
			int popIdx = 0
		) {
			int rf(0);
			for (auto& it : curSol) {
				rf |= evaluate(it, popIdx);
			}
			if (m_problem->hasTag(ProblemTag::kDOP)) {
				for (auto& it : indis) {
					rf |= evaluate(*it, popIdx);
				}
			}
			return rf;
		}

		virtual int calFitness(SolutionType& curSol,
			int pop_idx = 0){
			return evaluate(curSol, pop_idx);
		}

		virtual int updateMemory(
			const std::vector<SolutionType*>& indis,
			int popIdx = 0
		) {
			return 0;
		}
	
	};
}

#endif