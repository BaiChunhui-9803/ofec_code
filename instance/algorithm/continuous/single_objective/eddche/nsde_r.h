/********* Begin Register Information **********
{
    "description": "NSDE with restart strategy",
    "identifier": "NSDE_R",
    "name": "NSDE-R",
    "tags": [
        "single-objective",
        "continuous"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_NSDE_R_H
#define OFEC_NSDE_R_H

#include "../multi_modal/nmde/nsde.h"

namespace ofec {
	class NSDE_R : public NSDE {
		OFEC_CONCRETE_INSTANCE(NSDE_R)
	protected:
		size_t m_max_stagant_iters;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_NSDE_R_H
