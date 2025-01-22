/********* Begin Register Information **********
{
	"name": "CMAES-sampling",
	"identifier": "CMEAS_Sampling",
	"tags": [ "continuous", "single-objective" ]
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

#ifndef OFEC_SAMPLING_CMAES_INSTANCE_H
#define OFEC_SAMPLING_CMAES_INSTANCE_H


#include "../../sampling_algorithm_template.h"
#include "../../../../continuous/single_objective/global/cma_es/cma_es.h"
#include "../../../../continuous/single_objective/global/cma_es/cma_es_pop.h"


namespace ofec {
	class CMEAS_Sampling : virtual public SamplingAlgorithm, virtual public CMA_ES {
		OFEC_CONCRETE_INSTANCE(CMEAS_Sampling)

	protected:
		PopCMA_ES m_pop;

		virtual void evolve(Environment* env) override {
			CMA_ES::evolve(env);
		}
		virtual void getCurrentSolutions(std::vector<const SolutionBase*>& sols, std::vector<std::shared_ptr<SolutionInfo>>& solInfos, Environment* env)const override {
			sols.clear();
			SolutionInfo curInfo;
			curInfo.m_iter = m_iteration - 1;
			curInfo.m_popId = 0;
			curInfo.m_indiId = 0;
			for (auto& it : m_pop) {
				sols.push_back(it.get());
				solInfos.emplace_back(new SolutionInfo(curInfo));
				curInfo.m_indiId++;
			}
		}
	public:

		void addInputParameters() {}
		void initialize_(Environment* env) override {
			SamplingAlgorithm::initialize_(env);
			CMA_ES::initialize_(env);
			//m_pop.resize(m_pop_size, env);
			//m_pop.initialize(env, m_random.get());
			//m_pop.reproduce(env, m_random.get());
			//m_pop.evaluate(env);
			recordEachGeneration(env);
			++m_iteration;
		}


		//virtual void samplingDuringRun(Environment* env)override {
		//	while (!terminating()) {
		//		m_iteration++;
		//		m_pop.evolve(env, m_random.get());
		//		recordEachGeneration(env);
		//		if (m_pop.conditionNumber() > 1e14 || m_pop.equalFitness()) {
		//			terminate();
		//		}

		//		datumUpdated(env);
		//	}
		//}


	};
}

#endif