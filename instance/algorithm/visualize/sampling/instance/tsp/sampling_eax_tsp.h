/********* Begin Register Information **********
{
	"name": "EAX-TSP-sampling",
	"identifier": "EAX_TSP_Sampling",
	"tags": [ "travelling salesman problem" ]
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
*  Created by Diao Yiya on 2024-04-02
*
*************************************************************************/

#ifndef OFEC_SAMPLING_EAX_TSP_INSTANCE_H
#define OFEC_SAMPLING_EAX_TSP_INSTANCE_H


#include "../../sampling_algorithm_template.h"
#include "../../../../combination/eax_tsp/eax_tsp_origin/eax_tsp_alg.h"
#include "../../../../../../core/environment/environment.h"

namespace ofec {
	class EAX_TSP_Sampling : virtual public SamplingAlgorithm, virtual public EAX_TSP {
		OFEC_CONCRETE_INSTANCE(EAX_TSP_Sampling)

	protected:
		

		virtual void evolve(Environment* env) override{
			EAX_TSP::evolve(env);
		}
		virtual void getCurrentSolutions(std::vector<const SolutionBase*>& sols, std::vector<std::shared_ptr<SolutionInfo>>& solInfos, Environment* env)const override {
			auto& curpop = m_env->getCurPop();
			//sols.resize(curpop.size());
			sols.clear();
			SolutionInfo curInfo;
			curInfo.m_iter = m_iteration;
			curInfo.m_popId = 0;
			curInfo.m_indiId = 0;
			for (auto& it : curpop) {
				it.updateSolBase();
				sols.push_back(&it);
				solInfos.emplace_back(new SolutionInfo(curInfo));
				curInfo.m_indiId++;
			}
		}
	public:

		void addInputParameters() {}
		void initialize_(Environment* env) override {
			SamplingAlgorithm::initialize_(env);
			EAX_TSP::initialize_(env);
		}


	};
}

#endif