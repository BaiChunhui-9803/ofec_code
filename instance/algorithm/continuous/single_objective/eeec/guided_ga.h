/********* Begin Register Information **********
{
    "description": "",
    "identifier": "GuidedGA",
    "name": "Guided-GA",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_GUIDED_GA_H
#define OFEC_GUIDED_GA_H

#include "sbx_ga.h"

namespace ofec {
	class GuidedGA : public SBX_GA {
		OFEC_CONCRETE_INSTANCE(GuidedGA)
	protected:
		void addInputParameters();
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_GUIDED_GA_H
