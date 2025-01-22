/********* Begin Register Information **********
{
	"name": "Imbalance-GA",
	"identifier": "ImbalanceGA",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_IMBALANCE_GA_H
#define OFEC_IMBALANCE_GA_H

#include "sbx_ga.h"

namespace ofec {
	class ImbalanceGA : public SBX_GA {
		OFEC_CONCRETE_INSTANCE(ImbalanceGA)
	protected:
		void addInputParameters();
		void run_(Environment *env) override;

		Real m_shrink_rate;
		bool m_quarter_change;
	};
}

#endif // !OFEC_IMBALANCE_GA_H
