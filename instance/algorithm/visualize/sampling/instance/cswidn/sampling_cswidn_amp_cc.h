/********* Begin Register Information **********
{
	"name": "AMP-GL-SaDE-sampling",
	"identifier": "AMP_GL_SaDE_Sampling",
	"tags": [ "contamination source identification for water distribution network" ]
}
*********** End Register Information **********/

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
*  Created by Diao Yiya on 2024-07-10
*
*************************************************************************/

#ifndef OFEC_SAMPLING_ALGORITHM_INSTANCE_H
#define OFEC_SAMPLING_ALGORITHM_INSTANCE_H


#include "../../sampling_algorithm_template.h"
#include "../../../../realworld/csiwdn/amp_cc/amp_cc_instance.h"
#include "../../../../../../core/environment/environment.h"

namespace ofec {



	class AMP_GL_SaDE_Sampling : virtual public SamplingAlgorithm, virtual public AMP_GL_SaDE {
		OFEC_CONCRETE_INSTANCE(AMP_GL_SaDE_Sampling)
	public:

		struct SolInfo : public SolutionInfo {
			std::pair<int, int> m_source_idx;
			int m_current_phase = 0;
			int m_cc_popId = 0;


			// for output to file 
			virtual void toParameterVariants(ParameterVariantStream& out)const {
				SolutionInfo::toParameterVariants(out);
				out << m_source_idx.first << m_source_idx.second<< m_current_phase << m_cc_popId;
			}
			virtual void fromParameterVariants(ParameterVariantStream& in) {
				SolutionInfo::fromParameterVariants(in);
				in >> m_source_idx.first >> m_source_idx.second >> m_current_phase >> m_cc_popId;
			}
		};

	protected:

		virtual void initialize_(Environment* env) {
			SamplingAlgorithm::initialize_(env);
			AMP_GL_SaDE::initialize_(env);
			SamplingAlgorithm::m_iteration = 0;
			recordEachGeneration(env);

		}

		virtual void evolve(Environment* env) {
			AMP_GL_SaDE::evolve(env);
		}
		
		//virtual void evolve() override {
		//	m_env->evolve();
		//}
		virtual void getCurrentSolutions(std::vector<const SolutionBase*>& sols, std::vector<std::shared_ptr<SolutionInfo>>& solInfos, Environment* env)const override {
			auto pro = env->problem();

			SolInfo curInfo;
			curInfo.m_iter = SamplingAlgorithm::m_iteration;
			curInfo.m_current_phase = CAST_CSIWDN(pro)->phase();
			curInfo.m_source_idx = m_source_idx;
			for (size_t popId(0); popId < m_subpop.size(); ++popId) {
				auto& curPop = m_subpop.at(popId);
				curInfo.m_popId = popId;

				{
					auto& subpop = curPop.first_pop();
					curInfo.m_cc_popId = 1;
					for (size_t solId(0); solId < subpop.size(); ++solId) {
						curInfo.m_indiId = solId;
						solInfos.emplace_back(new SolInfo(curInfo));
						sols.push_back(&subpop.at(solId));
					}
				}
				{
					auto& subpop = curPop.second_pop();
					curInfo.m_cc_popId = 2;
					for (size_t solId(0); solId < subpop.size(); ++solId) {
						curInfo.m_indiId = solId;
						solInfos.emplace_back(new SolInfo(curInfo));
						sols.push_back(&subpop.at(solId));
					}
				}
			}

			//sols.resize(curpop.size());
			//sols.clear();
			//for (auto& it : curpop) {
			//	it.updateSolBase();
			//	sols.push_back(&it);
			//}
		}
		void addInputParameters() {}
	public:
		//virtual void samplingDuringRun(Environment* env)override {

		//	initializePop(env);
		//	recordEachGeneration(env);
		//	int tag = kNormalEval;
		//	while (!terminating() && tag != kTerminate) {
		//		tag= AMP_GL_SaDE::evolve(env);
		//		recordEachGeneration(env);
		//	}
		//}
	};
}

#endif