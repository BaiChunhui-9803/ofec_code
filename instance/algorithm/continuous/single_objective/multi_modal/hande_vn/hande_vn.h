/********* Begin Register Information **********
{
    "description": "History archive assisted niching differential evolution with variable\nneighborhood",
    "identifier": "HANDE_VN",
    "name": "HANDE/VN",
    "tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/


#ifndef OFEC_HANDE_VN_H
#define OFEC_HANDE_VN_H

#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class HANDE_VN : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(HANDE_VN)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_HANDE_VN_H
