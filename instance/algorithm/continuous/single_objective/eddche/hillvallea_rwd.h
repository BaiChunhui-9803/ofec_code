/********* Begin Register Information **********
{
    "description": "HillVallEA that restarts without doubling #samples",
    "identifier": "HillVallEA_RWD",
    "name": "HillVallEA-RWD",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_HILLVALLEA_RWD_H
#define OFEC_HILLVALLEA_RWD_H

#include "../multi_modal/hill_vallea/hill_vallea.h"

namespace ofec {
	class HillVallEA_RWD : public HillVallEA {
		OFEC_CONCRETE_INSTANCE(HillVallEA_RWD)
	protected:
		void addInputParameters();
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_HILLVALLEA_RWD_H
