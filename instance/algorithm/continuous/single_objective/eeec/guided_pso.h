/********* Begin Register Information **********
{
    "description": "",
    "identifier": "GuidedPSO",
    "name": "Guided-PSO",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_GUIDED_PSO_H
#define OFEC_GUIDED_PSO_H

#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class GuidedPSO : public Algorithm {
		OFEC_CONCRETE_INSTANCE(GuidedPSO)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_GUIDED_PSO_H
