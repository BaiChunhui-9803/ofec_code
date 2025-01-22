/********* Begin Register Information **********
{
	"name": "HillVallEA-sampling",
	"identifier": "HillVallEA_Sampling",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*-------------------------------------------------------------------------------
*
*  Created by Diao Yiya on 2024-10-03
*
*************************************************************************/

#ifndef OFEC_SAMPLING_HillVallEA_INSTANCE_H
#define OFEC_SAMPLING_HillVallEA_INSTANCE_H


#include "../../../../sampling_algorithm_multipop.h"
#include "../../../../../../continuous/single_objective/multi_modal/hill_vallea/hill_vallea.h"
#include "../../../../../../../../datum/algorithm/hill_vallea/hill_vallea_popinfo.h"


namespace ofec {

#define CAST_HillValleaSampling(alg) dynamic_cast<HillVallEA_Sampling*>(alg)

	class HillVallEA_Sampling : public SamplingAlgorithm<HillVallEA> {
		OFEC_CONCRETE_INSTANCE(HillVallEA_Sampling)


		std::vector<HillValleaPopInfo> m_popInfos;
	protected:

	public:
		const std::vector< HillValleaPopInfo>& getPopInfos() const{
			return m_popInfos;
		}

		virtual void initialize_(Environment* env) override {
			SamplingAlgorithm<HillVallEA>::initialize_(env);
			m_popInfos.clear();
		}
		

		virtual void handleDatumUpdated(Environment* env) override {
			SamplingAlgorithm<HillVallEA>::handleDatumUpdated(env);
			if (g_hillvallue_popinfo.updateFlag()) {
				g_hillvallue_popinfo.m_iter = m_iteration;
				m_popInfos.push_back(g_hillvallue_popinfo);
			}
		}

	};

}
#endif

