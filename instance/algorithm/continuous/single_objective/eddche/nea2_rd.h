/********* Begin Register Information **********
{
    "description": "NEA2 that doubles #samples upon restart",
    "identifier": "NEA2_RD",
    "name": "NEA2-RD",
    "tags": [
        "single-objective",
        "continuous"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_NEA2_RD_H
#define OFEC_NEA2_RD_H

#include "../multi_modal/nea2/nea2.h"

namespace ofec {
	class NEA2_RD : public NEA2 {
		OFEC_CONCRETE_INSTANCE(NEA2_RD)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_NEA2_RD_H
