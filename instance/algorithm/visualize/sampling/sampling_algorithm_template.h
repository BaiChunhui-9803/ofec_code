/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*-------------------------------------------------------------------------------
*
*  Created by Diao Yiya on 2023 11 13
*  Modified by Diao Yiya on 2024 07 10 
*
*-------------------------------------------------------------------------------
*************************************************************************/

#ifndef OFEC_SAMPLING_ALGORITHM_TEMPLATE_H
#define OFEC_SAMPLING_ALGORITHM_TEMPLATE_H

#include "../../../../core/algorithm/algorithm.h"
#include "../../../../core/environment/environment.h"
#include "../../../../core/parameter/parameter_variant.h"


namespace ofec {

#define CAST_SAMPLING_ALG(alg) dynamic_cast<SamplingAlgorithm*>(alg)

	class SamplingAlgorithm: virtual public Algorithm{
	public:
		
		struct SolutionInfo {
			int m_runId = 0;
			int m_iter = 0;
			int m_popId = 0;
			int m_indiId = 0;

			virtual ~SolutionInfo() = default;

			// for output to file 
			virtual void toParameterVariants(ParameterVariantStream& out)const {
				out << m_runId << m_iter << m_popId << m_indiId;
			}
			virtual void fromParameterVariants(ParameterVariantStream& in) {
				in >> m_runId >> m_iter >> m_popId >> m_indiId;
			}

		};


	protected:
		std::vector<std::shared_ptr<SolutionBase>> m_sols;
		std::vector<std::shared_ptr<SolutionInfo>> m_solInfos;

		int m_iteration = 0;
	protected:

		virtual void evolve() {}
		virtual void getCurrentSolutions(std::vector<const SolutionBase*>& sols, std::vector<std::shared_ptr<SolutionInfo>>& solInfos, Environment* env)const = 0;

		virtual void recordEachGeneration(Environment* env) {
			std::vector<const SolutionBase*> sols;
			std::vector<std::shared_ptr<SolutionInfo>> solInfos;
			getCurrentSolutions(sols, solInfos,env);
		//	std::vector<std::shared_ptr<SolutionBase>> totalSols;
			for (auto& it : sols) {
				std::shared_ptr<ofec::SolutionBase> cursol(env->problem()->createSolution(*it));
				m_sols.emplace_back(cursol);
			}
			for (auto& it : solInfos) {
				m_solInfos.emplace_back(it);
				
			}
			//m_sols.emplace_back(totalSols);
		}

	public:
		virtual void samplingDuringRun(Environment* env) {
			while (!terminating()) {
				evolve();
				recordEachGeneration(env);
			}
		}
		double getRandom() {
			return m_random->uniform.next();
		}
		const std::vector<std::shared_ptr<SolutionBase>>& getSols() {
			return m_sols;
		}
		const std::vector<std::shared_ptr<SolutionInfo>>& getSolInfos() {
			return m_solInfos;
		}


	//	SamplingAlgorithm() = default;
	};
}

#endif