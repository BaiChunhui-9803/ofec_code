/********* Begin Register Information **********
{
    "description": "kNPSO with restart-and-double strategy",
    "identifier": "kNPSO_RD",
    "name": "kNPSO-RD",
    "tags": [
        "single-objective",
        "continuous"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_KNPSO_RD_H
#define OFEC_KNPSO_RD_H

#include "../eeec/knpso.h"

namespace ofec {
	class kNPSO_RD : public kNPSO {
		OFEC_CONCRETE_INSTANCE(kNPSO_RD)
	protected:
		size_t m_max_stagant_iters;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_KNPSO_RD_H
