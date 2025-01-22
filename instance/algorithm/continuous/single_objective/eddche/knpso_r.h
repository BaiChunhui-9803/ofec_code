/********* Begin Register Information **********
{
    "description": "kNPSO with restart strategy",
    "identifier": "kNPSO_R",
    "name": "kNPSO-R",
    "tags": [
        "single-objective",
        "continuous"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_KNPSO_R_H
#define OFEC_KNPSO_R_H

#include "../eeec/knpso.h"

namespace ofec {
	class kNPSO_R : public kNPSO {
		OFEC_CONCRETE_INSTANCE(kNPSO_R)
	protected:
		size_t m_max_stagant_iters;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_KNPSO_R_H
