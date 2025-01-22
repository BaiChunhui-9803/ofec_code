/********* Begin Register Information **********
{
	"name": "HBE",
	"identifier": "HBE",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_HBE_H
#define OFEC_HBE_H

#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class HBE : public Algorithm {
		OFEC_CONCRETE_INSTANCE(HBE)
	protected:
		size_t m_num_increase_each_iter;

		void addInputParameters();
		void run_(Environment *env) override;
	};
}

#endif // !OFEC_HBE_H
