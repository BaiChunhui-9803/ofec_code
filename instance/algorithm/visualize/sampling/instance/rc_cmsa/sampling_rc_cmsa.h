/********* Begin Register Information **********
{
	"name": "RS-CMSA-sampling",
	"identifier": "RS_CMSA_Sampling",
	"tags": [ "continuous", "single-objective" ],
	"dependency on libraries": [ "Eigen" ]
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

#ifndef OFEC_SAMPLING_RS_CMSA_INSTANCE_H
#define OFEC_SAMPLING_RS_CMSA_INSTANCE_H


#include "../../sampling_algorithm_template.h"
#include "../../../../continuous/single_objective/multi_modal/rs_cmsa/rs_cmsa.h"


namespace ofec {
	class RS_CMSA_Sampling : virtual public SamplingAlgorithm, virtual public RS_CMSA {
		OFEC_CONCRETE_INSTANCE(RS_CMSA_Sampling)

	protected:

		virtual void evolve(Environment* env) override{
			RS_CMSA::evolve(env);
		}

		virtual void getCurrentSolutions(std::vector<const SolutionBase*>& sols, std::vector<std::shared_ptr<SolutionInfo>>& solInfos, Environment* env)const override {
			sols.clear();


			SolutionInfo curInfo;
			curInfo.m_iter = m_iteration;
			curInfo.m_popId = 0;
			curInfo.m_indiId = 0;
			
			for (int idx(0); idx < m_subpops.size(); ++idx) {
				curInfo.m_popId = idx;
				auto& curpop = m_subpops[idx];

				for (auto& it : curpop) {
					sols.push_back(it.get());
					solInfos.emplace_back(new SolutionInfo(curInfo));
					curInfo.m_indiId++;
				}
			}

		}
	public:

		void addInputParameters() {}
		void initialize_(Environment* env) override {
			SamplingAlgorithm::initialize_(env);
			RS_CMSA::initialize_(env);
		}



	};
}

#endif