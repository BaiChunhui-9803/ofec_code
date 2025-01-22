/********* Begin Register Information **********
{
	"name": "AMP_CC_CSWIDN-sampling",
	"identifier": "AMP_CC_CSWIDN_Data_Sampling",
	"tags": [ "contamination source identification for water distribution network" ]
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
*  Created by Diao Yiya on 2024-09-16
*
*************************************************************************/

#ifndef OFEC_SAMPLING_AMP_CC_CSWIDN_INSTANCE_H
#define OFEC_SAMPLING_AMP_CC_CSWIDN_INSTANCE_H


#include "../../../sampling_algorithm_multipop.h"
#include "../../../../../../algorithm/realworld/csiwdn/amp_cc/amp_cc_instance.h"

namespace ofec {

	namespace CSWIDN_samping {
		struct SolutionInfo : public SamplingData::SolutionInfo {
			std::pair<int, int> m_source_idx;
			int m_current_phase = 0;
			int m_cc_popId = 0;


			// for output to file 
			virtual void toParameterVariants(ParameterVariantStream& out)const {
				SamplingData::SolutionInfo::toParameterVariants(out);
				out << m_source_idx.first << m_source_idx.second << m_current_phase << m_cc_popId;
			}
			virtual void fromParameterVariants(ParameterVariantStream& in) {
				SamplingData::SolutionInfo::fromParameterVariants(in);
				in >> m_source_idx.first >> m_source_idx.second >> m_current_phase >> m_cc_popId;
			}
		};
	}


	class AMP_CC_CSWIDN_Data_Sampling : public  SamplingAlgorithm<AMP_GL_SaDE> {
		OFEC_CONCRETE_INSTANCE(AMP_CC_CSWIDN_Data_Sampling)

			
	protected:
		virtual void handleDatumUpdated(Environment* env) override {
			if (g_multi_pop.updateFlag()) {
				auto pro = env->problem();

				++m_iteration;
				CSWIDN_samping::SolutionInfo curInfo;
				curInfo.m_iter = m_iteration;
				curInfo.m_current_phase = CAST_CSIWDN(pro)->phase();
				curInfo.m_source_idx = CAST_CSIWDN(pro)->sourceIdx();
				curInfo.m_popId = 0;
				curInfo.m_indiId = 0;
				auto& pops = g_multi_pop.pops;
				for (int idpop(0); idpop < pops.size(); ++idpop) {
					curInfo.m_popId = idpop;
					for (int idIndi(0); idIndi < pops[idpop].size(); ++idIndi) {
						curInfo.m_indiId = idIndi;
						m_sols.emplace_back(pro->createSolution(*pops[idpop][idIndi]));
						m_solInfos.emplace_back(new SolutionInfo(curInfo));
					}

				}
			}
		}

	};
}


#endif