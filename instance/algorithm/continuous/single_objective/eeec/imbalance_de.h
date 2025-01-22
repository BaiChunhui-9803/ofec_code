/********* Begin Register Information **********
{
	"identifier": "ImbalanceDE",
	"name": "Imbalance-DE",
	"tags": [ "single-objective", "continuous" ]
}
*********** End Register Information **********/


#ifndef OFEC_IMBALANCE_DE_H
#define OFEC_IMBALANCE_DE_H

#include "canonical_de.h"

namespace ofec {
	class ImbalanceDE : public CanonicalDE {
		OFEC_CONCRETE_INSTANCE(ImbalanceDE)
	protected:
		void addInputParameters();
		void run_(Environment *env) override;

		Real m_shrink_rate;
		bool m_quarter_change;
	};
}

#endif // ! OFEC_IMBALANCE_DE_H
