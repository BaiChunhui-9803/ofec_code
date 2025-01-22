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
#ifndef OFEC_EVALUATION_STRATEGY_SAMPLES_H
#define OFEC_EVALUATION_STRATEGY_SAMPLES_H

#include<vector>
#include "evaluation_hls.h"
#include"../../combination/sequence/SP/sp_interpreter.h"
namespace ofec {

	template<typename TSolution>
	class EvaluationStrategySamples: EvaluationStrategyHLS<TSolution>{
	public:
		using SolutionType = typename TSolution;
	protected:
		int m_samples_size = 100;
	public:
		virtual int evaluate(SolutionType& curSol, int popIdx = 0) override{
			int rf(0);
			double obj(0);
			for (int idx(0); idx < m_samples_size; ++idx) {
				rf |= curSol.evaluate(pro, alg, true);
				obj += curSol.objective()[0];
			}
			curSol.objective()[0] = obj / m_samples_size;
			return rf;
		}
		//virtual void update(Problem *pro, Algorithm *alg, Random *rnd, const SolutionType& otherSol, SolutionType& curSol)const {}

	};
}

#endif