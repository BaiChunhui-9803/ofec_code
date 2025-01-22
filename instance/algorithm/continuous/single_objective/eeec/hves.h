/********* Begin Register Information **********
{
    "description": "",
    "identifier": "HVES",
    "name": "HVES",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_HVES_H
#define OFEC_HVES_H

#include "../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class HVES : public Algorithm {
		OFEC_CONCRETE_INSTANCE(HVES)
	protected:
		size_t m_num_samples, m_num_gradations;

		void addInputParameters();
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_HVES_H
