/********* Begin Register Information **********
{
    "description": "",
    "identifier": "ImbalancePSO",
    "name": "Imbalance-PSO",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_IMBALANCE_PSO_H
#define OFEC_IMBALANCE_PSO_H

#include "spso11.h"

namespace ofec {
	class ImbalancePSO : public SPSO11 {
		OFEC_CONCRETE_INSTANCE(ImbalancePSO)
	protected:
		void addInputParameters();
		void run_(Environment *env) override;
	private:
		Real m_shrink_rate;
		bool m_quarter_change;
	};
}

#endif // ! OFEC_IMBALANCE_PSO_H
