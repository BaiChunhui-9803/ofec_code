/********* Begin Register Information **********
{
    "identifier": "RestartingDE",
    "name": "Restarting-DE",
	"tags": [ "single-objective", "continuous" ]
}
*********** End Register Information **********/


#ifndef OFEC_RESTARING_DE_H
#define OFEC_RESTARING_DE_H

#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class RestartingDE : public Algorithm {
		OFEC_CONCRETE_INSTANCE(RestartingDE)
	protected:
		void addInputParameters();
		void run_(Environment *env) override;

        size_t m_pop_size;
	};
}

#endif // ! OFEC_RESTARING_DE_H
